"""Quality-gate ablation (plan section 5.4): dock identical ligands with and without each safeguard."""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

from . import engines, gates
from .ligprep import Ligand

CONDITIONS = ("baseline", "charges_x4_uncapped", "charges_x4_capped", "torsions_frozen")


def _variant(l: Ligand, condition: str, factor: float, cap: float) -> Ligand:
    if condition == "baseline":
        return l
    if condition == "charges_x4_uncapped":
        return replace(l, pdbqt=gates.scale_charges(l.pdbqt, factor))
    if condition == "charges_x4_capped":
        return replace(l, pdbqt=gates.cap_charges(gates.scale_charges(l.pdbqt, factor), cap))
    if condition == "torsions_frozen":
        return replace(l, pdbqt=gates.freeze_torsions(l.pdbqt))
    raise ValueError(condition)


def run_ablation(rec, ligands, workdir: Path, nrun: int = 20, seed: int = 42, factor: float = 4.0,
                 cap: float = gates.CHARGE_CAP, n_boot: int = 1000):
    """Returns (per-ligand energies, summary by condition). The ligands must have prepared PDBQT."""
    ok = [l for l in ligands if l.ok]
    per = {}
    for cond in CONDITIONS:
        variants = [_variant(l, cond, factor, cap) for l in ok]
        df, _ = engines.run_autodock_gpu(rec, variants, Path(workdir) / cond, nrun=nrun, seed=seed)
        per[cond] = df.set_index("name")["ad_energy"]
    wide = pd.DataFrame(per).loc[[l.name for l in ok]]
    rows = []
    rng = np.random.default_rng(0)
    for cond in CONDITIONS:
        e = wide[cond].to_numpy(float)
        flagged = int(sum(
            gates.charge_gate(_variant(l, cond, factor, cap).pdbqt)["charge_flag"] for l in ok))
        frozen = int(sum(gates.torsion_gate(l.smiles, gates.torsion_count(_variant(l, cond, factor, cap).pdbqt))["frozen_flexible"]
                         for l in ok))
        d = e - wide["baseline"].to_numpy(float)
        ok_d = np.isfinite(d)
        boots = [np.nanmean(d[rng.integers(0, ok_d.sum(), ok_d.sum())]) for _ in range(n_boot)] if ok_d.sum() > 1 else [np.nan]
        rows.append({"condition": cond, "n": int(np.isfinite(e).sum()), "mean_dG": float(np.nanmean(e)),
                     "median_dG": float(np.nanmedian(e)), "min_dG": float(np.nanmin(e)), "max_dG": float(np.nanmax(e)),
                     "non_binding_dG_ge_0": int(np.nansum(e >= 0)), "dG_below_minus20": int(np.nansum(e < -20)),
                     "n_charge_flagged": flagged, "n_torsion_frozen_flagged": frozen,
                     "mean_shift_vs_baseline": float(np.nanmean(d)),
                     "shift_ci_low": float(np.nanquantile(boots, 0.025)), "shift_ci_high": float(np.nanquantile(boots, 0.975))})
    return wide, pd.DataFrame(rows)
