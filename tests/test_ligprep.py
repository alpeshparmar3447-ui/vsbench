import re

import numpy as np
import pytest
from rdkit import Chem

from vsbench import ligprep as L


def _types(pdbqt):
    return {l.split()[-1] for l in pdbqt.splitlines() if l.startswith(("ATOM", "HETATM"))}


def _charge_sum(pdbqt):
    return sum(float(l[70:76]) for l in pdbqt.splitlines() if l.startswith(("ATOM", "HETATM")))


def _coords(pdbqt):
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                     for l in pdbqt.splitlines() if l.startswith(("ATOM", "HETATM"))])


@pytest.mark.parametrize("smiles, expected_types", [
    ("CCO", {"C", "HD", "OA"}),
    ("c1ccccc1", {"A"}),
    ("CC(=O)Nc1ccc(O)cc1", {"A", "C", "N", "OA", "HD"}),
])
def test_autodock_atom_types(smiles, expected_types):
    lig = L.prepare_ligand(smiles, "x")
    assert lig.ok and _types(lig.pdbqt) == expected_types


@pytest.mark.parametrize("smiles, formal", [("CC(=O)[O-]", -1), ("C[N+](C)(C)C", 1), ("CCO", 0),
                                            ("C[NH3+]", 1), ("[O-][N+](=O)c1ccccc1", 0)])
def test_partial_charges_sum_to_the_formal_charge(smiles, formal):
    lig = L.prepare_ligand(smiles, "x")
    assert lig.ok and _charge_sum(lig.pdbqt) == pytest.approx(formal, abs=0.02)


def test_torsion_count_is_consistent_and_heavy_atoms_are_reported():
    lig = L.prepare_ligand("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "ibuprofen")
    assert lig.n_heavy == 15
    assert lig.n_torsions == sum(1 for l in lig.pdbqt.splitlines() if l.startswith("BRANCH"))
    assert f"TORSDOF {lig.n_torsions}" in lig.pdbqt
    assert lig.n_torsions >= 4
    assert Chem.MolFromMolBlock(lig.molblock, removeHs=False).GetNumHeavyAtoms() == 15


@pytest.mark.parametrize("smiles, fragment", [
    ("this is not a smiles", "parse"),
    ("[Na+].[Cl-]", "fragments"),
    ("C", "fewer than two heavy atoms"),
])
def test_failures_are_recorded_not_raised(smiles, fragment):
    lig = L.prepare_ligand(smiles, "bad")
    assert lig.ok is False and fragment in lig.error and lig.pdbqt == ""


def test_preparation_is_deterministic_and_seed_dependent():
    a = L.prepare_ligand("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "n1", 42)
    b = L.prepare_ligand("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "n1", 42)
    c = L.prepare_ligand("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "n1", 43)
    assert a.pdbqt == b.pdbqt
    assert a.pdbqt != c.pdbqt
    assert len(_coords(a.pdbqt)) == len(_coords(c.pdbqt))


def test_ligand_seed_depends_on_name_and_base_seed_only():
    assert L.ligand_seed("144203670", 42) == L.ligand_seed("144203670", 42)
    assert L.ligand_seed("144203670", 42) != L.ligand_seed("144203993", 42)
    assert L.ligand_seed("144203670", 42) != L.ligand_seed("144203670", 43)
    assert 0 <= L.ligand_seed("x", 42) <= 0x7FFFFFFF


def test_prepare_many_preserves_order_and_matches_serial_results():
    items = [(s, f"id{i}") for i, s in enumerate(
        ["CCO", "c1ccccc1O", "CC(=O)Nc1ccc(O)cc1", "bad smiles", "CCN(CC)CC", "C[N+](C)(C)C",
         "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "OC(=O)c1ccccc1", "CCCCCC", "c1ccncc1"] * 2)]
    items = [(s, f"{n}_{k}") for k, (s, n) in enumerate(items)]
    serial = L.prepare_many(items, 42, 1)
    parallel = L.prepare_many(items, 42, 4)
    assert [l.name for l in parallel] == [n for _, n in items]
    assert [l.pdbqt for l in serial] == [l.pdbqt for l in parallel]
    assert [l.ok for l in serial].count(False) == 2


def test_embedding_timeout_is_applied_and_does_not_break_normal_molecules():
    assert L.EMBED_TIMEOUT_S > 0
    lig = L.prepare_ligand("CC(=O)Nc1ccc(O)cc1", "t", 42, timeout=1)
    assert lig.ok
    mol = L.embed("CCO", 1, timeout=1)
    assert mol.GetNumConformers() == 1


def test_a_pathological_molecule_cannot_stall_preparation():
    # A polycyclic limonoid (epoxide, furan, eight stereocentres) from LIT-PCBA MAPK1 took
    # 390 s and then failed before the embedding was bounded; it must now return quickly,
    # either prepared or recorded as a failure.
    import time
    smiles = ("CC(=O)OC1CC(OC(=O)C)[C@]2(C)C(CC(OC(=O)C)[C@]3(C)C2CC[C@@]4(C)C(C[C@H]5O[C@]345)"
              "c6cocc6)C1(C)C")
    t0 = time.time()
    lig = L.prepare_ligand(smiles, "26750089", 42, timeout=1)
    assert time.time() - t0 < 30
    assert lig.ok or "embedding failed" in lig.error
