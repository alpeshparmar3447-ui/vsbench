"""Receptor-level assessment of every candidate template: structure quality, binding-site
detection with druggability scoring, and the metal-ion audit.

These three studies describe the *receptors*, not the screen, so they run over every template
of a target rather than over the one template redocking selected, and their tables are written
once per data set rather than once per run tier (see :func:`receptor_dir`).

Nothing here is tied to a particular benchmark. A target is any directory the loader accepts,
a template is any protein/ligand pair inside it, and the one service that needs a public
accession, DoGSiteScorer, is used only for templates whose identifier actually looks like a PDB
code. Configure the whole group through the ``receptors`` key of a run config::

    "receptors": {"fpocket": true, "dogsite": false, "p2rank": false,
                  "controls": true, "restore_metals": true, "n_jobs": 2}
"""
from __future__ import annotations

import re
import traceback
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from . import data, metals as K, sites as S, structure_qc as Q
from .pipeline import Config

PDB_CODE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")

DEFAULTS = {
    "fpocket": True,        # local, no network
    "dogsite": False,       # public web service, public PDB codes only, minutes per structure
    "p2rank": "auto",       # local, needs Java: "auto" enables it when both are found
    "controls": True,       # the four fixed druggability control structures
    "restore_metals": True,  # compare each template against its deposition
    "offline": False,       # True forbids every outbound request (see options())
    "n_jobs": 2,            # concurrent templates; keep low when a web service is involved
}


def options(cfg: Config) -> dict:
    """Receptor-study options: the defaults, the config on top, then local resolution.

    Two keys are resolved rather than taken literally. ``p2rank: "auto"`` becomes True only
    when both the P2Rank launcher and a Java runtime are actually found, so a config can ask
    for the local detector without knowing whether this machine has one. ``offline: true``
    forces off everything that would reach the network, which at present is DoGSiteScorer and
    the deposition download, so that a run provably makes no outbound request.
    """
    opt = {**DEFAULTS, **(getattr(cfg, "receptors", None) or {})}

    if opt.get("p2rank") == "auto" or (opt.get("p2rank") and not opt.get("p2rank_bin")):
        from . import tools
        prank, java = tools.find_p2rank(), tools.find_java()
        if prank and java:
            opt["p2rank"] = True
            opt.setdefault("p2rank_bin", str(prank))
            opt.setdefault("java_bin", str(java))
        elif opt.get("p2rank") == "auto":
            opt["p2rank"] = False
        # an explicit p2rank:true with nothing found is left on, so the run reports the failure

    if opt.get("offline"):
        opt["dogsite"] = False
    return opt


def receptor_dir(cfg: Config) -> Path:
    """Where the template-level tables live.

    They describe the data set and not the run, so by default they sit one level above the
    run's own results directory and are shared by every tier run against the same data.
    ``results_dir`` itself is used when the two are the same place.
    """
    explicit = (getattr(cfg, "receptors", None) or {}).get("dir")
    out = Path(explicit) if explicit else Path(cfg.results_dir).parent
    out.mkdir(parents=True, exist_ok=True)
    return out


def pdb_cache(cfg: Config) -> Path:
    """Where deposited structures are cached.

    The config may name it. Otherwise an existing ``pdb_cache`` directory at or above the data
    root is reused, because the depositions belong to the data set and not to any one run, and
    only when none exists is one created beside the data root.
    """
    explicit = (getattr(cfg, "receptors", None) or {}).get("pdb_cache")
    if explicit:
        out = Path(explicit)
    else:
        root = Path(cfg.data_root).resolve()
        found = next((p / "pdb_cache" for p in (root, *root.parents)
                      if (p / "pdb_cache").is_dir()), None)
        out = found or root.parent / "pdb_cache"
    out.mkdir(parents=True, exist_ok=True)
    return out


def mol2_heavy_xyz(path) -> np.ndarray:
    """Heavy-atom coordinates of a mol2 file, as an (n, 3) array."""
    xyz, in_atoms = [], False
    for line in Path(path).read_text().splitlines():
        if line.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            continue
        if line.startswith("@<TRIPOS>BOND"):
            break
        if in_atoms and line.strip():
            t = line.split()
            if len(t) >= 6 and not t[5].startswith("H"):
                xyz.append((float(t[2]), float(t[3]), float(t[4])))
    return np.array(xyz) if xyz else np.empty((0, 3))


def _templates(cfg: Config, targets) -> list[tuple[str, str]]:
    root = Path(cfg.data_root)
    out = []
    for t in targets:
        tg = data.load_target(root, t)
        out.extend((t, tm) for tm in tg.templates)
    return out


def all_targets(cfg: Config) -> list[str]:
    """Every target directory under the data root, not only those the config screens."""
    root = Path(cfg.data_root)
    return sorted(p.name for p in root.iterdir()
                  if p.is_dir() and (p / "actives.smi").exists())


def _cached_deposition(template: str, cache: Path) -> Path | None:
    for ext in ("pdb", "cif"):
        p = cache / f"{template.upper()}.{ext}"
        if p.exists():
            return p
    return None


# --------------------------------------------------------------- structure quality
def qc_study(cfg: Config, targets, log=print) -> pd.DataFrame:
    """Backbone geometry, chain breaks, engineered residues and resolution per template."""
    cache = pdb_cache(cfg)
    rows = []
    for target, tm in _templates(cfg, targets):
        tg = data.load_target(Path(cfg.data_root), target)
        q = Q.qc_template(tm, tg.protein_mol2(tm), mol2_heavy_xyz(tg.ligand_mol2(tm)),
                          _cached_deposition(tm, cache)).as_dict()
        q["target"] = target
        rows.append(q)
    df = pd.DataFrame(rows)
    if df.empty:
        log("no templates to assess")
        return df
    df["rama_outlier_template"] = False
    for target, g in df.groupby("target"):
        bad = Q.flag_rama_outliers(dict(zip(g.template, g.rama_outside)))
        df.loc[df.template.isin(bad) & (df.target == target), "rama_outlier_template"] = True
    out = receptor_dir(cfg) / "structure_qc.csv"
    df.to_csv(out, index=False)
    log(f"structure QC: {len(df)} templates -> {out}")
    return df


# --------------------------------------------------- binding site and druggability
def _site_rows(cfg: Config, opt: dict, work: Path, dogsite_cache: Path, args) -> list[dict]:
    """Pocket rows for one template. Never raises: a template that cannot be read, or a
    detector that fails, becomes a row carrying the error, because one bad template must not
    discard the detections already made for every other one."""
    target, tm = args
    try:
        tg = data.load_target(Path(cfg.data_root), target)
        lig = mol2_heavy_xyz(tg.ligand_mol2(tm))
        native = S.native_site(Q.read_mol2_residues(tg.protein_mol2(tm)), lig)
        pdb = work / f"{tm}.pdb"
        S.write_pdb_from_mol2(tg.protein_mol2(tm), pdb)
        base = {"target": target, "n_native_residues": len(native),
                "ligand_buried_fraction": S.ligand_buried_fraction(pdb, lig)}
    except Exception as e:
        return [{"template": tm, "target": target, "method": "setup",
                 "error": f"{type(e).__name__}: {e}"[:200]}]
    out = []

    def attempt(method, call):
        try:
            out.append({"template": tm, "method": method,
                        **S.site_recovery(call(), native), **base})
        except Exception as e:                                   # one detector failing is data
            out.append({"template": tm, "method": method,
                        "error": f"{type(e).__name__}: {e}"[:200], **base})

    if opt["fpocket"]:
        attempt("fpocket", lambda: S.run_fpocket(pdb, work / "fpocket" / tm))
    if opt["p2rank"] and opt.get("p2rank_bin") and opt.get("java_bin"):
        attempt("p2rank", lambda: S.run_p2rank(pdb, work / "p2rank" / tm,
                                               Path(opt["p2rank_bin"]), Path(opt["java_bin"])))
    # DoGSiteScorer resolves a public accession, so it can only see templates named by one
    if opt["dogsite"] and PDB_CODE.match(tm):
        attempt("dogsite", lambda: S.run_dogsite(tm, dogsite_cache))
    return out


def _control_rows(cfg: Config, opt: dict, work: Path, dogsite_cache: Path, item) -> list[dict]:
    pdb_id, (role, resname, desc) = item
    try:
        text = K.fetch_deposition(pdb_id, pdb_cache(cfg),
                                  offline=opt.get("offline", False)).read_text()
        native, lig = S.control_native_site(text, resname)
        heavy_pdb = S.heavy_protein_pdb(text, work / f"{pdb_id}.pdb")
        base = {"target": "CONTROL", "control_role": role, "description": desc,
                "n_native_residues": len(native),
                "ligand_buried_fraction": S.ligand_buried_fraction(heavy_pdb, lig)}
        rows = []

        def attempt(method, call):
            try:
                rows.append({"template": pdb_id, "method": method,
                             **S.site_recovery(call(), native), **base})
            except Exception as e:
                rows.append({"template": pdb_id, "method": method,
                             "error": f"{type(e).__name__}: {e}"[:200], **base})

        # Every enabled detector must face the controls, or its druggability score is reported
        # without ever having been held to the rule that decides whether it means anything.
        if opt["fpocket"]:
            attempt("fpocket", lambda: S.run_fpocket(heavy_pdb, work / "fpocket" / pdb_id))
        if opt["p2rank"] and opt.get("p2rank_bin") and opt.get("java_bin"):
            attempt("p2rank", lambda: S.run_p2rank(heavy_pdb, work / "p2rank" / pdb_id,
                                                   Path(opt["p2rank_bin"]), Path(opt["java_bin"])))
        if opt["dogsite"]:
            attempt("dogsite", lambda: S.run_dogsite(pdb_id, dogsite_cache))
        return rows
    except Exception as e:
        return [{"template": pdb_id, "method": "error", "target": "CONTROL",
                 "error": f"{type(e).__name__}: {e}"[:200]}]


def sites_study(cfg: Config, targets, log=print, controls=None) -> pd.DataFrame:
    """Detect pockets in every template, score their druggability, and run the controls.

    Writes ``sites.csv`` (one row per template and detector: whether the detector's top pocket,
    top three, or any pocket recovered the site the co-crystal ligand occupies, plus that
    pocket's druggability score) and, when controls are enabled, ``site_controls.csv``.

    A druggability score is called informative only where the detector ranks the negative
    control below every positive, which is checked when the controls table is read.

    One column is not reproducible to the decimal: fpocket estimates ``site_pocket_volume`` by
    Monte Carlo and returns values differing by a few parts in a thousand between runs on
    byte-identical input. Pocket count, ranking, recovery and druggability score are stable, so
    every column a conclusion rests on is reproducible; the volume is reported as an estimate.
    """
    opt = options(cfg)
    out_dir = receptor_dir(cfg)
    # Beside the tables, not inside the run: pocket output and the DoGSite response cache
    # describe the templates, so every tier run against this data set reuses them.
    work = out_dir / "work" / "sites"
    work.mkdir(parents=True, exist_ok=True)
    # The DoGSite responses are keyed by accession, so they outlive any one data set; a config
    # may point several data sets at one cache rather than re-querying the service.
    dogsite_cache = Path(opt.get("dogsite_cache") or out_dir / "work" / "dogsite")
    jobs = _templates(cfg, targets)
    if not jobs:
        log("no templates to assess")
        return pd.DataFrame()
    detectors = [k for k in ("fpocket", "p2rank", "dogsite") if opt[k]]
    log(f"sites: {len(jobs)} templates × {', '.join(detectors) or 'no detector enabled'}")

    rows: list[dict] = []
    dest = out_dir / "sites.csv"
    with ThreadPoolExecutor(max(1, int(opt["n_jobs"]))) as ex:
        for i, r in enumerate(ex.map(lambda a: _site_rows(cfg, opt, work, dogsite_cache, a), jobs), 1):
            rows += r
            pd.DataFrame(rows).to_csv(dest, index=False)       # resumable: written as it goes
            if i % 10 == 0:
                log(f"  {i}/{len(jobs)} templates")
    log(f"sites -> {dest}")

    if opt["controls"]:
        items = (controls or S.CONTROLS).items()
        crows = [r for item in items for r in _control_rows(cfg, opt, work, dogsite_cache, item)]
        cdest = out_dir / "site_controls.csv"
        pd.DataFrame(crows).to_csv(cdest, index=False)
        log(f"druggability controls -> {cdest}")
    return pd.DataFrame(rows)


def druggability_verdict(controls: pd.DataFrame) -> dict:
    """Per detector: is its druggability score usable on this data set?

    The rule fixed in advance is that a score counts as informative only where the detector
    scores the negative control below every positive control. Returns the verdict and the
    scores behind it, so the caller never has to restate the rule.
    """
    out = {}
    if controls is None or controls.empty or "site_pocket_druggability" not in controls:
        return out
    for method, g in controls.groupby("method"):
        pos = g[g.control_role == "positive"]["site_pocket_druggability"].dropna()
        neg = g[g.control_role == "negative"]["site_pocket_druggability"].dropna()
        if pos.empty or neg.empty:
            continue
        out[method] = {"informative": bool(neg.max() < pos.min()),
                       "negative_control": float(neg.max()),
                       "lowest_positive": float(pos.min()),
                       "n_positive": int(len(pos))}
    return out


# --------------------------------------------------------------------- metal audit
#: Functional ions this far from the ligand are still reported, so that a structural ion
#: outside the pocket is visible rather than silently absent.
WIDE_ION_CUTOFF = 30.0


def _metal_rows(cfg: Config, opt: dict, args) -> list[dict]:
    """One row per ion worth retaining, or a single row when the template has none."""
    target, tm = args
    tg = data.load_target(Path(cfg.data_root), target)
    try:
        protein = tg.protein_mol2(tm)
        lig = mol2_heavy_xyz(tg.ligand_mol2(tm))
        ca = K.read_ca_coordinates(protein)
        in_file = K.read_structure_ions(protein)
        base = dict(target=target, template=tm, n_ca_mol2=len(ca),
                    ions_in_benchmark_file=len(in_file))

        if not (opt["restore_metals"] and PDB_CODE.match(tm)):
            # No accession to audit against: report what the file itself holds and say so.
            return [{**base, "frame_overlap": np.nan, "ions_in_deposition": np.nan,
                     "relevant_ions": len(K.ions_near(in_file, lig, K.DEFAULT_POCKET_CUTOFF)),
                     "note": "no deposition checked"}]

        dep = K.fetch_deposition(tm, pdb_cache(cfg), offline=opt.get("offline", False))
        dep_ions = K.read_structure_ions(dep)
        base.update(frame_overlap=round(K.frame_overlap(ca, K.read_ca_coordinates(dep)), 3),
                    ions_in_deposition=len(dep_ions))
        wide = K.ions_near([i for i in dep_ions if i.element in K.RETAIN_ELEMENTS],
                           lig, WIDE_ION_CUTOFF)
        base["nearest_functional_ion_within_30A"] = (
            f"{wide[0][0].element}@{wide[0][1]:.1f}" if wide else "")

        relevant = K.relevant_ions(dep_ions, lig)
        if not relevant:
            return [{**base, "relevant_ions": 0}]
        atoms = (list(K.read_pdb_protein_atoms(dep.read_text(errors="ignore")))
                 if dep.suffix == ".pdb" else None)
        rows = []
        for ion, dist in relevant:
            donors = len(K.coordination(ion.xyz, atoms)) if atoms is not None else -1
            nearest_ca = (float(np.min(np.linalg.norm(ca - np.array(ion.xyz), axis=1)))
                          if len(ca) else np.nan)
            rows.append({**base, "relevant_ions": len(relevant), "element": ion.element,
                         "ion": ion.key(), "dist_to_ligand_A": round(dist, 2),
                         "protein_donors_le_2.8A": donors,
                         "dist_to_nearest_CA_in_benchmark_A": round(nearest_ca, 1)})
        return rows
    except Exception as e:
        return [dict(target=target, template=tm, error=f"{type(e).__name__}: {e}"[:200])]


def metal_study(cfg: Config, targets, log=print) -> pd.DataFrame:
    """Audit every template for metal ions against its original deposition.

    Answers three things per template: do the benchmark file and the deposition share a
    coordinate frame, which ions does the benchmark file itself contain, and which functional
    ions sit near the co-crystal ligand in the deposition. A template whose identifier is not
    a PDB code is reported from its own file alone, with the omission stated in a ``note``.
    """
    opt = options(cfg)
    jobs = _templates(cfg, targets)
    if not jobs:
        log("no templates to assess")
        return pd.DataFrame()
    with ThreadPoolExecutor(max(1, int(opt.get("fetch_jobs", 8)))) as ex:
        rows = [r for group in ex.map(lambda a: _metal_rows(cfg, opt, a), jobs) for r in group]
    df = pd.DataFrame(rows)
    out = receptor_dir(cfg) / "metal_audit.csv"
    df.to_csv(out, index=False)
    with_ion = int(df.get("relevant_ions", pd.Series(dtype=float)).fillna(0).gt(0).sum())
    log(f"metal audit: {len(jobs)} templates, {with_ion} rows with a pocket ion -> {out}")
    return df
