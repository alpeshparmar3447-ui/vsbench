"""Rebuild RDKit molecules from docked PDBQT poses.

Meeko records the atom correspondence between the PDBQT and a SMILES string in
``REMARK SMILES`` and ``REMARK SMILES IDX`` lines (pairs of SMILES atom index,
1-based, and PDBQT atom serial). Using that map, a pose becomes an RDKit
molecule whose heavy-atom conformer is the docked geometry, ready for
:func:`vsbench.rmsd.symmetry_rmsd`.
"""
from __future__ import annotations

import re

from rdkit import Chem
from rdkit.Geometry import Point3D

_PREFIXES = ("INPUT-LIGAND-PDBQT: ", "DOCKED: ")


def _strip(line: str) -> str:
    for p in _PREFIXES:
        if line.startswith(p):
            return line[len(p):]
    return line


def parse_mapping(text: str) -> tuple[str, dict[int, int]]:
    """Return (smiles, {pdbqt_serial: rdkit_atom_index}) from REMARK lines."""
    smiles = None
    pairs: list[int] = []
    for raw in text.splitlines():
        line = _strip(raw)
        if line.startswith("REMARK SMILES IDX"):
            pairs += [int(t) for t in line.split()[3:]]
        elif line.startswith("REMARK SMILES ") and smiles is None:
            smiles = line.split(None, 2)[2].strip()
    if smiles is None or not pairs or len(pairs) % 2:
        raise ValueError("no usable REMARK SMILES / SMILES IDX mapping found")
    return smiles, {pairs[i + 1]: pairs[i] - 1 for i in range(0, len(pairs), 2)}


def parse_coordinates(pose_lines) -> dict[int, tuple[float, float, float]]:
    """Atom serial -> (x, y, z) from the ATOM/HETATM records of one pose."""
    xyz = {}
    for raw in pose_lines:
        line = _strip(raw)
        if line.startswith(("ATOM", "HETATM")):
            xyz[int(line[6:11])] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    return xyz


def pose_to_mol(mapping_text: str, pose_lines) -> Chem.Mol:
    """Heavy-atom RDKit molecule carrying the docked coordinates of one pose."""
    smiles, serial_to_idx = parse_mapping(mapping_text)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("could not parse the SMILES recorded in the PDBQT")
    coords = parse_coordinates(pose_lines)
    conf = Chem.Conformer(mol.GetNumAtoms())
    seen = set()
    for serial, idx in serial_to_idx.items():
        if serial not in coords:
            raise ValueError(f"pose lacks atom serial {serial}")
        conf.SetAtomPosition(idx, Point3D(*coords[serial]))
        seen.add(idx)
    if len(seen) != mol.GetNumAtoms():
        raise ValueError(f"mapping covers {len(seen)} of {mol.GetNumAtoms()} heavy atoms")
    mol.AddConformer(conf, assignId=True)
    return mol


_MODEL_RE = re.compile(r"^(?:DOCKED: )?MODEL\s+(\d+)")


def split_models(text: str) -> list[list[str]]:
    """Split a multi-model PDBQT (or the DOCKED lines of a DLG) into per-model line lists."""
    models, cur = [], None
    for line in text.splitlines():
        if _MODEL_RE.match(line):
            cur = []
        elif line.strip().endswith("ENDMDL"):
            if cur is not None:
                models.append(cur)
            cur = None
        elif cur is not None:
            cur.append(line)
    return models
