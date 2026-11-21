"""Ligand quality gates and their ablation (plan sections 4.8 and 5.4).

Two failure modes of docking-ready ligands are checked:
* a flexible ligand written with every torsion frozen (it is then docked as a rigid body);
* implausibly large partial charges (they can swamp the electrostatic term).
Gates only flag. Ablation helpers build the inputs needed to measure what each gate is worth:
charges scaled up, charges capped, torsions frozen.
"""
from __future__ import annotations

import re

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

CHARGE_CAP = 1.0          # e: |q| above this on any atom flags the ligand
MIN_ROTB_FOR_FREEZE = 3   # a ligand with this many strict rotatable bonds should not have zero torsions
_ATOM = ("ATOM", "HETATM")


def pdbqt_charges(pdbqt: str) -> list[float]:
    return [float(l[70:76]) for l in pdbqt.splitlines() if l.startswith(_ATOM)]


def strict_rotatable_bonds(smiles: str) -> int:
    m = Chem.MolFromSmiles(smiles)
    return rdMolDescriptors.CalcNumRotatableBonds(m, strict=True) if m is not None else 0


def torsion_gate(smiles: str, n_torsions: int) -> dict:
    """Flag a ligand whose torsion tree is empty although it has several rotatable bonds.

    Amide and similar bonds are legitimately inactive in the torsion tree, so the gate compares with
    zero, not with the rotatable-bond count: a deficit is expected, a total absence is a defect.
    """
    rotb = strict_rotatable_bonds(smiles)
    return {"rotb_strict": rotb, "n_torsions": n_torsions,
            "frozen_flexible": bool(n_torsions == 0 and rotb >= MIN_ROTB_FOR_FREEZE)}


def charge_gate(pdbqt: str, cap: float = CHARGE_CAP) -> dict:
    q = np.abs(pdbqt_charges(pdbqt))
    return {"max_abs_charge": float(q.max()) if len(q) else 0.0, "charge_flag": bool(len(q) and q.max() > cap),
            "net_charge": float(sum(pdbqt_charges(pdbqt)))}


def apply_gates(ligands) -> list[dict]:
    """Gate results for prepared ligands (objects with .name, .smiles, .ok, .pdbqt, .n_torsions)."""
    out = []
    for l in ligands:
        row = {"name": l.name, "frozen_flexible": False, "rotb_strict": 0, "max_abs_charge": float("nan"),
               "charge_flag": False}
        if l.ok:
            row.update(torsion_gate(l.smiles, l.n_torsions))
            row.update({k: v for k, v in charge_gate(l.pdbqt).items() if k in ("max_abs_charge", "charge_flag")})
        out.append(row)
    return out


# ---------------------------------------------------------------------------- ablation inputs
def _rewrite_charges(pdbqt: str, fn) -> str:
    out = []
    for l in pdbqt.splitlines():
        if l.startswith(_ATOM):
            l = l[:70] + f"{fn(float(l[70:76])):+6.3f}" + l[76:]
        out.append(l)
    return "\n".join(out) + "\n"


def scale_charges(pdbqt: str, factor: float) -> str:
    """Multiply every partial charge (the artificial stress used to trigger the charge gate)."""
    return _rewrite_charges(pdbqt, lambda q: max(-9.999, min(9.999, q * factor)))


def cap_charges(pdbqt: str, cap: float = CHARGE_CAP) -> str:
    """Clamp every partial charge to [-cap, +cap]."""
    return _rewrite_charges(pdbqt, lambda q: max(-cap, min(cap, q)))


def freeze_torsions(pdbqt: str) -> str:
    """The same ligand with its torsion tree removed: every atom in one rigid ROOT, TORSDOF 0."""
    head = [l for l in pdbqt.splitlines() if l.startswith("REMARK")]
    atoms = [l for l in pdbqt.splitlines() if l.startswith(_ATOM)]
    return "\n".join(head + ["ROOT"] + atoms + ["ENDROOT", "TORSDOF 0"]) + "\n"


def torsion_count(pdbqt: str) -> int:
    m = re.search(r"^TORSDOF\s+(\d+)", pdbqt, re.M)
    return int(m.group(1)) if m else 0
