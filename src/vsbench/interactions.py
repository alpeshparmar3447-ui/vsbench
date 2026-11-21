"""Protein-ligand interaction analysis (plan section 4.11), with PLIP.

An interaction is the pair (type, chain, residue number). The reference is the co-crystal ligand; a
docked pose is scored by how much of that reference it reproduces. The same protein file serves
reference and poses so that any protonation choices cancel.
"""
from __future__ import annotations

import contextlib
import io
from pathlib import Path

from rdkit import Chem

LIG_RES, LIG_CHAIN = "LIG", "Z"


def ligand_pdb_lines(mol: Chem.Mol, conf_id: int = -1) -> list[str]:
    """HETATM lines (plus CONECT) for a ligand, with a fixed residue name and chain."""
    block = Chem.MolToPDBBlock(mol, confId=conf_id)
    out = []
    for l in block.splitlines():
        if l.startswith(("ATOM", "HETATM")):
            l = "HETATM" + l[6:17] + f"{LIG_RES:>3s} {LIG_CHAIN}   1" + l[26:]
            out.append(l)
        elif l.startswith("CONECT"):
            out.append(l)
    return out


def write_complex(protein_pdb: Path, ligand: Chem.Mol, out: Path) -> Path:
    prot = [l for l in Path(protein_pdb).read_text().splitlines() if l.startswith(("ATOM", "TER"))]
    Path(out).write_text("\n".join(prot + ligand_pdb_lines(ligand)) + "\nEND\n")
    return Path(out)


def interaction_set(complex_pdb: Path) -> frozenset:
    """Residue-level interactions of the LIG ligand in a complex file."""
    from plip.basic import config as plip_config
    from plip.structure.preparation import PDBComplex
    plip_config.VERBOSE = False
    plip_config.QUIET = True
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        c = PDBComplex()
        # PLIP writes a repaired copy of the complex through tempfile.mktemp and never removes
        # it, so the default of the system temp directory leaves one file per pose behind: a
        # screen of this size left 10,000 of them, 2.6 GB, on a tmpfs that other work shares.
        # Directing them beside the complex puts them in the caller's own temporary directory,
        # which is removed with it.
        c.output_path = str(Path(complex_pdb).parent)
        c.load_pdb(str(complex_pdb))
        c.analyze()
    lig = next((l for l in c.ligands if l.hetid == LIG_RES), None)
    if lig is None:
        return frozenset()
    inter = c.interaction_sets[":".join([lig.hetid, lig.chain, str(lig.position)])]
    kinds = {
        "hbond": list(inter.hbonds_ldon) + list(inter.hbonds_pdon),
        "hydrophobic": list(inter.hydrophobic_contacts),
        "pi_stack": list(inter.pistacking),
        "pi_cation": list(inter.pication_laro) + list(inter.pication_paro),
        "salt_bridge": list(inter.saltbridge_lneg) + list(inter.saltbridge_pneg),
        "halogen": list(inter.halogen_bonds),
    }
    out = set()
    for kind, items in kinds.items():
        for it in items:
            out.add((kind, getattr(it, "reschain", "A"), int(it.resnr)))
    for it in inter.metal_complexes:
        out.add(("metal", getattr(it, "reschain", "A"), int(getattr(it, "resnr", 0))))
    return frozenset(out)


def recovery(reference: frozenset, pose: frozenset) -> dict:
    """How much of the reference interaction set the pose reproduces."""
    inter = reference & pose
    union = reference | pose
    return {"n_reference": len(reference), "n_pose": len(pose), "n_shared": len(inter),
            "recall": len(inter) / len(reference) if reference else float("nan"),
            "precision": len(inter) / len(pose) if pose else float("nan"),
            "tanimoto": len(inter) / len(union) if union else float("nan")}


def residue_level(interactions: frozenset) -> frozenset:
    """The set of contacted residues regardless of interaction type."""
    return frozenset((c, r) for _, c, r in interactions)


# ------------------------------------------------------------------ batch scoring
def _interaction_job(args):
    """Worker: interaction set of one pose given its pre-built ligand lines (picklable)."""
    import tempfile
    protein_lines, lig_lines, name = args
    with tempfile.TemporaryDirectory(prefix="plip_") as d:
        path = Path(d) / "c.pdb"
        path.write_text("\n".join(protein_lines + lig_lines) + "\nEND\n")
        try:
            return name, interaction_set(path), ""
        except Exception as exc:                                    # PLIP/OpenBabel failure on one pose
            return name, None, f"{type(exc).__name__}: {exc}"[:200]


def protein_lines(protein_pdb: Path) -> list[str]:
    return [l for l in Path(protein_pdb).read_text().splitlines() if l.startswith(("ATOM", "TER"))]


def reference_set(protein_pdb: Path, ligand: Chem.Mol) -> frozenset:
    """Interactions of the co-crystal ligand (the reference the poses are compared with)."""
    _, s, err = _interaction_job((protein_lines(protein_pdb), ligand_pdb_lines(ligand), "reference"))
    if s is None:
        raise RuntimeError(f"PLIP failed on the reference complex: {err}")
    if not s:
        # A deposited ligand sits in its own pocket by construction, so an empty reference means the
        # complex was built wrongly, not that the site is featureless. Recovery against it is undefined
        # and every pose would be scored NaN, which reads like a result rather than a broken input.
        raise RuntimeError(f"PLIP found no interaction for the co-crystal ligand in {protein_pdb}: "
                           f"the receptor is missing residue identity or the ligand is misplaced")
    return s


def score_poses(protein_pdb: Path, reference: frozenset, poses: dict, n_jobs: int = 8):
    """Recovery of the reference interactions for each pose (name -> Mol). A pose PLIP cannot analyse, or
    a compound with no pose, is recorded with ok=False; the caller decides how to count it."""
    import pandas as pd
    from joblib import Parallel, delayed
    prot = protein_lines(protein_pdb)
    names = list(poses)
    jobs = [(prot, ligand_pdb_lines(poses[n]), n) for n in names]
    res = Parallel(n_jobs=n_jobs, backend="loky")(delayed(_interaction_job)(j) for j in jobs)
    ref_res = residue_level(reference)
    rows = []
    for name, s, err in res:
        if s is None:
            rows.append({"name": name, "ok": False, "error": err})
            continue
        r = recovery(reference, s)
        rr = residue_level(s)
        rows.append({"name": name, "ok": True, "error": "", "n_pose": r["n_pose"], "n_shared": r["n_shared"],
                     "recall": r["recall"], "precision": r["precision"], "tanimoto": r["tanimoto"],
                     "residue_recall": len(ref_res & rr) / len(ref_res) if ref_res else float("nan")})
    return pd.DataFrame(rows)
