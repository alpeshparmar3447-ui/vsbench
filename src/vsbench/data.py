"""Target loading, descriptors, seeded sampling and property matching.

The loader is not tied to LIT-PCBA. A target is any directory holding ``actives.smi``,
``inactives.smi`` and at least one ``<id>_protein.mol2`` / ``<id>_ligand.mol2`` pair, so a
user's own actives and confirmed inactives work unchanged; see README, "Running your own
target". ``validate_target_dir`` checks such a directory before any compute is spent on it.
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
from scipy.spatial import cKDTree

RDLogger.DisableLog("rdApp.*")

DESCRIPTOR_COLUMNS = ["mw", "clogp", "hbd", "hba", "rotb", "tpsa", "charge"]


@dataclass(frozen=True)
class Target:
    name: str
    root: Path
    actives: pd.DataFrame
    inactives: pd.DataFrame
    templates: tuple[str, ...]

    def protein_mol2(self, template: str) -> Path:
        return self.root / f"{template}_protein.mol2"

    def ligand_mol2(self, template: str) -> Path:
        return self.root / f"{template}_ligand.mol2"


def read_smi(path: Path) -> pd.DataFrame:
    """Read a whitespace-separated ``SMILES ID`` file into columns ``smiles`` and ``id``."""
    rows = []
    for n, line in enumerate(Path(path).read_text().splitlines(), 1):
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) != 2:
            raise ValueError(f"{path}:{n}: expected 'SMILES ID', got {len(parts)} fields")
        rows.append((parts[0], parts[1]))
    df = pd.DataFrame(rows, columns=["smiles", "id"])
    if df["id"].duplicated().any():
        raise ValueError(f"{path}: duplicate compound ids")
    return df


def load_target(data_root: Path, name: str) -> Target:
    root = Path(data_root) / name
    if not root.is_dir():
        raise FileNotFoundError(root)
    templates = tuple(sorted(p.name.split("_")[0] for p in root.glob("*_protein.mol2")))
    for t in templates:
        if not (root / f"{t}_ligand.mol2").exists():
            raise FileNotFoundError(f"missing ligand file for template {t}")
    return Target(name, root, read_smi(root / "actives.smi"), read_smi(root / "inactives.smi"),
                  templates)


@dataclass(frozen=True)
class TargetCheck:
    """What a candidate target directory contains, and what is wrong with it."""
    name: str
    n_actives: int = 0
    n_inactives: int = 0
    templates: tuple[str, ...] = ()
    n_unparsable: int = 0
    problems: tuple[str, ...] = ()

    @property
    def ok(self) -> bool:
        return not self.problems


def validate_target_dir(target_dir: Path, n_actives: int = 40, per_active_matched: int = 8,
                        n_random: int = 300) -> TargetCheck:
    """Check one target directory against the layout the pipeline needs, before any docking.

    Every file is opened and parsed rather than merely checked for existence: a structure file
    that exists but carries no readable residue, or a SMILES file whose compounds RDKit cannot
    build, fails only once a stage is well under way, by which time GPU time has been spent.
    """
    target_dir = Path(target_dir)
    name = target_dir.name
    if not target_dir.is_dir():
        return TargetCheck(name, problems=(f"{target_dir} is not a directory",))
    problems: list[str] = []
    counts = {}
    for stem in ("actives", "inactives"):
        path = target_dir / f"{stem}.smi"
        if not path.exists():
            problems.append(f"missing {stem}.smi")
            counts[stem] = None
            continue
        try:
            counts[stem] = len(read_smi(path))
        except ValueError as exc:
            problems.append(f"{stem}.smi unreadable: {exc}")
            counts[stem] = None

    proteins = sorted(p.name.split("_")[0] for p in target_dir.glob("*_protein.mol2"))
    templates = tuple(t for t in proteins if (target_dir / f"{t}_ligand.mol2").exists())
    for t in proteins:
        if t not in templates:
            problems.append(f"template {t}: {t}_protein.mol2 has no matching {t}_ligand.mol2")
    if not templates:
        problems.append("no <pdb>_protein.mol2 / <pdb>_ligand.mol2 pair found")

    # A receptor whose residue identity cannot be read is the failure that produced this
    # project's own Correction 1: pocket detection and interaction analysis both return
    # nothing, with no error from either.
    from .structure_qc import read_mol2_residues
    for t in templates:
        try:
            residues = read_mol2_residues(target_dir / f"{t}_protein.mol2")
        except Exception as exc:                                        # unreadable mol2
            problems.append(f"template {t}: protein mol2 unreadable: {type(exc).__name__}: {exc}")
            continue
        if not residues:
            problems.append(f"template {t}: no residue could be read from the protein mol2")
        elif all(r.name.upper() == "UNK" for r in residues):
            problems.append(f"template {t}: every residue reads as UNK, so chain and residue "
                            f"identity are lost")

    n_act_file, n_inact_file = counts.get("actives"), counts.get("inactives")
    n_unparsable = 0
    if n_act_file is not None:
        used = min(n_actives, n_act_file)
        if used < 10:
            problems.append(f"only {n_act_file} actives: too few for a bootstrap interval")
        desc = compute_descriptors(read_smi(target_dir / "actives.smi")["smiles"])
        n_unparsable = int(desc.isna().any(axis=1).sum())
        if n_unparsable:
            problems.append(f"{n_unparsable} of {n_act_file} actives cannot be built by RDKit")
        if n_inact_file is not None:
            need = n_random + used * per_active_matched
            if n_inact_file < need:
                problems.append(f"{n_inact_file} inactives, but {need} are needed for "
                                f"{n_random} random plus {used}x{per_active_matched} matched")
    return TargetCheck(name, n_act_file or 0, n_inact_file or 0, templates, n_unparsable,
                       tuple(problems))


def validate_data_root(data_root: Path, targets=None, **kwargs) -> list[TargetCheck]:
    """Validate every configured target, or every subdirectory if none is named."""
    root = Path(data_root)
    if not root.is_dir():
        raise FileNotFoundError(f"data root {root} does not exist")
    names = list(targets) if targets else sorted(p.name for p in root.iterdir() if p.is_dir())
    return [validate_target_dir(root / n, **kwargs) for n in names]


def _describe_one(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (np.nan,) * len(DESCRIPTOR_COLUMNS)
    return (
        Descriptors.MolWt(mol),
        Crippen.MolLogP(mol),
        Lipinski.NumHDonors(mol),
        Lipinski.NumHAcceptors(mol),
        rdMolDescriptors.CalcNumRotatableBonds(mol),
        rdMolDescriptors.CalcTPSA(mol),
        sum(a.GetFormalCharge() for a in mol.GetAtoms()),
    )


def compute_descriptors(smiles, n_jobs: int = 1) -> pd.DataFrame:
    """Six 1-D properties plus net charge for each SMILES (NaN if unparsable)."""
    smiles = list(smiles)
    if n_jobs > 1 and len(smiles) > 2000:
        with ProcessPoolExecutor(n_jobs) as ex:
            rows = list(ex.map(_describe_one, smiles, chunksize=500))
    else:
        rows = [_describe_one(s) for s in smiles]
    return pd.DataFrame(rows, columns=DESCRIPTOR_COLUMNS)


def sample_rows(df: pd.DataFrame, n: int, seed: int) -> pd.DataFrame:
    """Seeded sample without replacement; returns everything if n >= len(df)."""
    if n < 0:
        raise ValueError("n must be non-negative")
    if n >= len(df):
        return df.copy().reset_index(drop=True)
    rng = np.random.default_rng(seed)
    return df.iloc[np.sort(rng.choice(len(df), n, replace=False))].reset_index(drop=True)


def match_inactives(active_desc: pd.DataFrame, pool_desc: pd.DataFrame, per_active: int,
                    seed: int, return_owners: bool = False):
    """Pick pool rows whose properties resemble the actives, without replacement.

    Features are standardised on the pool; each active in turn (random order,
    repeated ``per_active`` times) takes its nearest unused pool compound. The
    result is a set of experimentally inactive compounds that no one-dimensional
    property can separate from the actives. Returns positional pool indices.
    """
    if per_active < 1:
        raise ValueError("per_active must be at least 1")
    a = active_desc[DESCRIPTOR_COLUMNS].to_numpy(float)
    p = pool_desc[DESCRIPTOR_COLUMNS].to_numpy(float)
    if np.isnan(a).any() or np.isnan(p).any():
        raise ValueError("descriptors contain NaN; drop unparsable structures first")
    need = len(a) * per_active
    if need > len(p):
        raise ValueError(f"pool has {len(p)} compounds but {need} are needed")
    mu, sd = p.mean(axis=0), p.std(axis=0)
    sd[sd == 0] = 1.0
    tree = cKDTree((p - mu) / sd)
    queries = (a - mu) / sd
    used = np.zeros(len(p), dtype=bool)
    rng = np.random.default_rng(seed)
    chosen, owners = [], []
    for _ in range(per_active):
        for i in rng.permutation(len(a)):
            k = min(64, len(p))
            while True:
                _, idx = tree.query(queries[i], k=k)
                idx = np.atleast_1d(idx)
                free = idx[~used[idx]]
                if len(free):
                    used[free[0]] = True
                    chosen.append(int(free[0]))
                    owners.append(int(i))
                    break
                if k >= len(p):
                    raise RuntimeError("no unused pool compound left")  # pragma: no cover
                k = min(len(p), k * 4)
    if return_owners:
        return np.asarray(chosen), np.asarray(owners)
    return np.asarray(chosen)


def standardized_mean_difference(a: pd.DataFrame, b: pd.DataFrame) -> pd.Series:
    """Per-property (mean_a - mean_b) / pooled SD. |SMD| < 0.1 is conventionally balanced."""
    out = {}
    for c in DESCRIPTOR_COLUMNS:
        va, vb = a[c].to_numpy(float), b[c].to_numpy(float)
        pooled = np.sqrt((va.var(ddof=1) + vb.var(ddof=1)) / 2.0)
        scale = max(1.0, abs(va.mean()), abs(vb.mean()))
        # A constant property has zero variance up to floating-point noise; dividing
        # by that noise would report a spurious, enormous difference.
        out[c] = 0.0 if pooled <= 1e-9 * scale else (va.mean() - vb.mean()) / pooled
    return pd.Series(out)


@dataclass(frozen=True)
class NegativeSets:
    """Actives plus two separately drawn sets of experimentally inactive compounds.

    ``random`` is a uniform draw from all inactives. ``matched`` is drawn from a
    disjoint candidate pool so that, compound for compound, its 1-D properties
    resemble the actives'. Descriptor columns are included in every frame.
    """
    actives: pd.DataFrame
    random: pd.DataFrame
    matched: pd.DataFrame
    n_unparsable: int
    decoy: pd.DataFrame | None = None       # DUD-E-protocol decoys (never tested against the target)

    def balance(self) -> pd.DataFrame:
        cols = {"random": standardized_mean_difference(self.actives, self.random),
                "matched": standardized_mean_difference(self.actives, self.matched)}
        if self.decoy is not None:
            cols["decoy"] = standardized_mean_difference(self.actives, self.decoy)
        return pd.DataFrame(cols)


def _with_descriptors(df: pd.DataFrame, n_jobs: int) -> pd.DataFrame:
    out = pd.concat([df.reset_index(drop=True), compute_descriptors(df["smiles"], n_jobs)], axis=1)
    return out


def draw_negative_sets(actives: pd.DataFrame, inactives: pd.DataFrame, n_random: int,
                       per_active: int, seed: int, pool_size: int = 20000,
                       n_jobs: int = 1, decoy_pool: pd.DataFrame | None = None,
                       decoy_per_active: int = 0) -> NegativeSets:
    """Draw the random and the property-matched negative sets (see NegativeSets)."""
    if n_random < 1 or per_active < 1:
        raise ValueError("n_random and per_active must be at least 1")
    perm = np.random.default_rng(seed).permutation(len(inactives))
    rand_rows = inactives.iloc[perm[:n_random]]
    pool_rows = inactives.iloc[perm[n_random:n_random + pool_size]]
    act = _with_descriptors(actives, n_jobs)
    rnd = _with_descriptors(rand_rows, n_jobs)
    pool = _with_descriptors(pool_rows, n_jobs)
    bad = int(act[DESCRIPTOR_COLUMNS].isna().any(axis=1).sum()
              + rnd[DESCRIPTOR_COLUMNS].isna().any(axis=1).sum()
              + pool[DESCRIPTOR_COLUMNS].isna().any(axis=1).sum())
    act, rnd, pool = (d.dropna(subset=DESCRIPTOR_COLUMNS).reset_index(drop=True) for d in (act, rnd, pool))
    need = len(act) * per_active
    if len(pool) < need:
        raise ValueError(f"candidate pool has {len(pool)} compounds; {need} needed for matching")
    idx = match_inactives(act, pool, per_active, seed)
    decoy = None
    if decoy_per_active and decoy_pool is not None:
        decoy = draw_dude_decoys(act, list(act["smiles"]), decoy_pool, decoy_per_active, seed, n_jobs=n_jobs)
    return NegativeSets(act, rnd, pool.iloc[idx].reset_index(drop=True), bad, decoy)


# ------------------------------------------------------------------ DUD-E-protocol decoys
def build_decoy_pool(data_root: Path, exclude_target: str, per_target: int = 4000, seed: int = 0) -> pd.DataFrame:
    """Compounds never tested against ``exclude_target``: inactives of the other targets under
    ``data_root``, minus anything that appears in this target's actives or inactives. They play
    the role ZINC compounds play in DUD-E: presumed, not demonstrated, to be inactive here."""
    root = Path(data_root)
    own = set()
    for name in ("actives.smi", "inactives.smi"):
        own |= set(read_smi(root / exclude_target / name)["id"])
    parts = []
    for i, d in enumerate(sorted(p for p in root.iterdir() if p.is_dir() and p.name != exclude_target)):
        df = read_smi(d / "inactives.smi")
        parts.append(sample_rows(df, per_target, seed + i))
    pool = pd.concat(parts, ignore_index=True).drop_duplicates("id")
    pool = pool[~pool["id"].isin(own)].drop_duplicates("smiles").reset_index(drop=True)
    return pool


def draw_dude_decoys(active_desc: pd.DataFrame, active_smiles, pool: pd.DataFrame, per_active: int, seed: int,
                     candidates: int = 100, drop_most_similar: float = 0.25, n_jobs: int = 1) -> pd.DataFrame:
    """Property-matched, topologically dissimilar decoys, after Mysinger et al. (DUD-E, 2012).

    Each active takes its ``candidates`` nearest unused pool compounds in standardised property space;
    the ``drop_most_similar`` fraction most similar (ECFP4 Tanimoto) to ANY active is discarded and the
    nearest remaining compound is taken, repeated ``per_active`` times.
    """
    from .similarity import max_similarity, morgan_fps
    pool = pool.reset_index(drop=True)
    pdesc = compute_descriptors(pool["smiles"], n_jobs)
    ok = ~pdesc.isna().any(axis=1)
    pool, pdesc = pool[ok].reset_index(drop=True), pdesc[ok].reset_index(drop=True)
    a = active_desc[DESCRIPTOR_COLUMNS].to_numpy(float)
    p = pdesc[DESCRIPTOR_COLUMNS].to_numpy(float)
    if len(a) * per_active > len(p):
        raise ValueError("decoy pool too small")
    mu, sd = p.mean(axis=0), p.std(axis=0)
    sd[sd == 0] = 1.0
    tree = cKDTree((p - mu) / sd)
    q = (a - mu) / sd
    act_fps = morgan_fps(list(active_smiles))
    pool_fps = morgan_fps(list(pool["smiles"]))
    sim_all = max_similarity(pool_fps, act_fps)                       # similarity of every pool compound to the actives
    used = np.zeros(len(p), dtype=bool)
    rng = np.random.default_rng(seed)
    chosen, owners = [], []
    for _ in range(per_active):
        for i in rng.permutation(len(a)):
            k = min(len(p), max(candidates * 4, 400))
            while True:
                _, idx = tree.query(q[i], k=k)
                idx = np.atleast_1d(idx)
                free = idx[~used[idx]][:candidates]
                if len(free) >= 4:
                    keep = int(np.ceil(len(free) * (1.0 - drop_most_similar)))
                    ranked = free[np.argsort(sim_all[free], kind="mergesort")][:keep]   # least similar first
                    nearest = min(ranked, key=lambda j: np.linalg.norm(p[j] - a[i]))
                    used[nearest] = True
                    chosen.append(int(nearest))
                    break
                if k >= len(p):
                    raise RuntimeError("no decoy candidates left")
                k = min(len(p), k * 4)
    out = pool.iloc[chosen].reset_index(drop=True)
    return pd.concat([out, pdesc.iloc[chosen].reset_index(drop=True)], axis=1)
