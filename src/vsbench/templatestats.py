"""What the redocking gate is actually measuring: template difficulty from the co-crystal ligand.

A template fails in-place redocking for a reason, and the reason is usually visible in its deposited ligand
before any docking runs. This module joins the redocking RMSDs already on disk to properties of each
template's co-crystal ligand and asks which properties predict failure.

The motivating case: nine of eleven MTORC1 templates in LIT-PCBA carry a 29-membered macrocycle of 914 to
980 Da. No rigid-receptor docking engine reproduces such a pose, so a validation gate applied to that target
is testing the structures, not the engines. Per PROTOCOL.md Amendment 3 this is reported as a property of the
benchmark's structures.

Correlations use Spearman's rho with a permutation p-value, so no distributional assumption is made and no
new dependency is needed. A ligand whose mol2 record cannot be sanitised keeps its heavy-atom count, is
counted in the per-target totals, and is left out of the correlations.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors as rd

from .gridstrat import wilson

RDLogger.DisableLog("rdApp.*")

MACROCYCLE_MIN_RING = 12        # IUPAC-style threshold; rapamycin's macrolactone is 29-membered
VALIDATED_RMSD = 2.0
PROPERTIES = ("heavy_atoms", "mol_weight", "rotatable_bonds", "n_rings", "max_ring_size")


def ligand_properties(mol2: Path) -> dict:
    """Properties of one deposited ligand. ``sanitized`` is False when only the atom count is trustworthy."""
    mol2 = Path(mol2)
    out = {"ligand_file": mol2.name, "sanitized": False,
           **{p: np.nan for p in PROPERTIES}, "macrocycle": None}
    if not mol2.exists():
        return out
    m = Chem.MolFromMol2File(str(mol2), sanitize=True)
    if m is None:
        raw = Chem.MolFromMol2File(str(mol2), sanitize=False)
        if raw is not None:
            out["heavy_atoms"] = float(raw.GetNumHeavyAtoms())
        return out
    rings = m.GetRingInfo().AtomRings()
    largest = max((len(r) for r in rings), default=0)
    out.update(sanitized=True, heavy_atoms=float(m.GetNumHeavyAtoms()),
               mol_weight=float(Descriptors.MolWt(m)),
               rotatable_bonds=float(rd.CalcNumRotatableBonds(m)),
               n_rings=float(rd.CalcNumRings(m)), max_ring_size=float(largest),
               macrocycle=bool(largest >= MACROCYCLE_MIN_RING))
    return out


def template_table(data_root: Path, results_dir: Path, targets) -> pd.DataFrame:
    """One row per template: its redocking RMSDs and its co-crystal ligand's properties."""
    data_root, results_dir = Path(data_root), Path(results_dir)
    rows = []
    for target in targets:
        rc = results_dir / target / "redock.csv"
        if not rc.exists():
            continue
        for _, r in pd.read_csv(rc).iterrows():
            pdb = str(r["template"])
            ad, gn = float(r.get("ad_lowest_rmsd", np.nan)), float(r.get("gn_top_rmsd", np.nan))
            worst = np.nanmax([ad, gn]) if np.isfinite([ad, gn]).any() else np.nan
            rows.append({"target": target, "template": pdb, "ran": bool(r.get("ok", False)),
                         "ad_rmsd": ad, "gn_rmsd": gn, "worst_rmsd": worst,
                         "validated": bool(np.isfinite(worst) and worst <= VALIDATED_RMSD),
                         **ligand_properties(data_root / target / f"{pdb}_ligand.mol2")})
    return pd.DataFrame(rows)


def spearman(x, y, n_perm: int = 2000, seed: int = 0) -> tuple[float, float, int]:
    """Spearman's rho with a two-sided permutation p-value. Pairs with a missing value are dropped."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    n = len(x)
    if n < 3 or len(np.unique(x)) < 2 or len(np.unique(y)) < 2:
        return float("nan"), float("nan"), n
    rx = pd.Series(x).rank().to_numpy()
    ry = pd.Series(y).rank().to_numpy()
    rho = float(np.corrcoef(rx, ry)[0, 1])
    rng = np.random.default_rng(seed)
    null = np.array([np.corrcoef(rx, rng.permutation(ry))[0, 1] for _ in range(n_perm)])
    p = float((np.sum(np.abs(null) >= abs(rho)) + 1) / (n_perm + 1))     # add-one: never reports p = 0
    return rho, p, n


def correlations(table: pd.DataFrame, n_perm: int = 2000, seed: int = 0) -> pd.DataFrame:
    """Spearman's rho between each ligand property and the worst-engine redocking RMSD."""
    use = table[table["sanitized"] & np.isfinite(table["worst_rmsd"])]
    rows = []
    for p in PROPERTIES:
        rho, pval, n = spearman(use[p], use["worst_rmsd"], n_perm=n_perm, seed=seed)
        rows.append({"property": p, "spearman_rho": rho, "p_permutation": pval, "n_templates": n})
    return pd.DataFrame(rows)


def by_macrocycle(table: pd.DataFrame) -> pd.DataFrame:
    """Validation rate split by the macrocycle flag, with Wilson intervals."""
    rows = []
    for flag, label in ((False, "non-macrocyclic"), (True, "macrocyclic")):
        sub = table[table["macrocycle"] == flag]
        n, k = len(sub), int(sub["validated"].sum())
        lo, hi = wilson(k, n) if n else (np.nan, np.nan)
        rows.append({"group": label, "n_templates": n, "validated": k,
                     "rate": k / n if n else np.nan, "ci_low": lo, "ci_high": hi})
    unknown = table[table["macrocycle"].isna()]
    rows.append({"group": "ligand not sanitizable", "n_templates": len(unknown),
                 "validated": int(unknown["validated"].sum()),
                 "rate": float(unknown["validated"].mean()) if len(unknown) else np.nan,
                 "ci_low": np.nan, "ci_high": np.nan})
    return pd.DataFrame(rows)


def per_target(table: pd.DataFrame) -> pd.DataFrame:
    """Per target: template counts, how many validated, and the median ligand size and flexibility."""
    rows = []
    for target, g in table.groupby("target", sort=False):
        rows.append({"target": target, "templates": len(g), "ran": int(g["ran"].sum()),
                     "validated": int(g["validated"].sum()),
                     "macrocyclic": int((g["macrocycle"] == True).sum()),          # noqa: E712 - NA-safe
                     "median_heavy_atoms": float(g["heavy_atoms"].median()),
                     "median_mol_weight": float(g["mol_weight"].median()),
                     "median_rotatable_bonds": float(g["rotatable_bonds"].median()),
                     "best_rmsd": float(g["worst_rmsd"].min()),
                     "has_non_macrocyclic_template": bool((g["macrocycle"] == False).any())})  # noqa: E712
    return pd.DataFrame(rows)
