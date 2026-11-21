"""In-place, symmetry-aware heavy-atom RMSD.

No superposition is performed: a pose is compared with the reference where it
sits, because a superposing metric measures conformational similarity and not
docking accuracy. Atom correspondence is found by graph isomorphism on the
element-labelled heavy-atom skeleton (bond orders ignored, so a SMILES-derived
molecule and a crystal-derived one can be compared), and the RMSD is minimised
over every isomorphism, which covers topological symmetry such as CF3 rotors,
para-substituted rings and symmetric piperazines.
"""
from __future__ import annotations

import warnings

import numpy as np
from rdkit import Chem


def _skeleton(mol: Chem.Mol, conf_id: int = -1):
    """Heavy-atom, single-bond, element-labelled copy plus its coordinates."""
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    index = {old: new for new, old in enumerate(heavy)}
    rw = Chem.RWMol()
    for old in heavy:
        rw.AddAtom(Chem.Atom(mol.GetAtomWithIdx(old).GetAtomicNum()))
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in index and j in index:
            rw.AddBond(index[i], index[j], Chem.BondType.SINGLE)
    skel = rw.GetMol()
    skel.UpdatePropertyCache(strict=False)
    pos = mol.GetConformer(conf_id).GetPositions()[heavy]
    return skel, pos


def symmetry_rmsd(ref: Chem.Mol, probe: Chem.Mol, max_matches: int = 100_000,
                  ref_conf: int = -1, probe_conf: int = -1) -> float:
    """Minimum in-place heavy-atom RMSD over all graph isomorphisms probe -> ref."""
    ref_skel, ref_xyz = _skeleton(ref, ref_conf)
    probe_skel, probe_xyz = _skeleton(probe, probe_conf)
    if ref_skel.GetNumAtoms() != probe_skel.GetNumAtoms():
        raise ValueError(
            f"heavy-atom counts differ: ref {ref_skel.GetNumAtoms()} vs probe {probe_skel.GetNumAtoms()}"
        )
    matches = ref_skel.GetSubstructMatches(probe_skel, uniquify=False, maxMatches=max_matches,
                                           useChirality=False)
    if not matches:
        raise ValueError("probe and reference are not the same molecular graph")
    if len(matches) >= max_matches:
        warnings.warn("isomorphism cap reached; RMSD is an upper bound", RuntimeWarning)
    m = np.asarray(matches)                      # m[k, i] = ref atom matched by probe atom i
    diff = ref_xyz[m] - probe_xyz[None, :, :]    # (n_matches, n_atoms, 3)
    msd = (diff ** 2).sum(axis=2).mean(axis=1)
    return float(np.sqrt(msd.min()))
