import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

from vsbench.box import MAX_GRID_POINTS, Box, box_from_coords
from vsbench.rmsd import symmetry_rmsd


def _embed(smiles, seed=7, hs=True):
    mol = Chem.MolFromSmiles(smiles)
    if hs:
        mol = Chem.AddHs(mol)
    assert AllChem.EmbedMolecule(mol, randomSeed=seed) == 0
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


def _move(mol, R=None, t=(0, 0, 0)):
    out = Chem.Mol(mol)
    conf = out.GetConformer()
    pos = conf.GetPositions()
    if R is not None:
        pos = pos @ R.T
    pos = pos + np.asarray(t)
    for i, p in enumerate(pos):
        conf.SetAtomPosition(i, p.tolist())
    return out


# ---------------------------------------------------------------- box
def test_box_center_and_cubic_edge():
    xyz = np.array([[0, 0, 0], [10, 0, 0], [0, 4, 0], [0, 0, 6.0]])
    box = box_from_coords(xyz, buffer=4.0, min_edge=10.0)
    assert box.center == pytest.approx(tuple(xyz.mean(axis=0)))
    assert box.size == (18.0, 18.0, 18.0)  # span 10 + 2*4


def test_box_min_edge_and_noncubic():
    xyz = np.array([[0, 0, 0], [3, 0, 0], [0, 2, 0], [0, 0, 1.0]])
    assert box_from_coords(xyz, buffer=1.0, min_edge=20.0).size == (20.0, 20.0, 20.0)
    box = box_from_coords(xyz, buffer=4.0, min_edge=5.0, cubic=False)
    assert box.size == (11.0, 10.0, 9.0)


def test_box_grid_points_are_even_and_limit_enforced():
    assert all(n % 2 == 0 for n in Box((0, 0, 0), (20.9, 20.9, 20.9)).npts)
    assert Box((0, 0, 0), (20.9, 20.9, 20.9)).npts == (56, 56, 56)
    huge = np.array([[0, 0, 0], [100, 0, 0], [0, 5, 0]])
    with pytest.raises(ValueError, match="grid points"):
        box_from_coords(huge)
    assert MAX_GRID_POINTS == 256


def test_box_needs_coordinates():
    with pytest.raises(ValueError):
        box_from_coords(np.zeros((2, 3)))


# --------------------------------------------------------------- RMSD
def test_identical_pose_is_zero_and_pure_translation_is_exact():
    ref = _embed("CC(=O)Nc1ccc(O)cc1")
    assert symmetry_rmsd(ref, ref) == pytest.approx(0.0, abs=1e-9)
    shifted = _move(ref, t=(3.0, 4.0, 0.0))
    assert symmetry_rmsd(ref, shifted) == pytest.approx(5.0, abs=1e-6)


def test_ring_flip_of_a_symmetric_ring_costs_nothing():
    # A 180-degree ring flip swaps the ortho pair and the meta pair together.
    # That permutation is a graph automorphism, so the symmetric RMSD is zero;
    # swapping only the ortho pair is not, and must register as a real change.
    ref = Chem.RemoveHs(_embed("c1ccc(cc1)C(=O)O"))
    ipso = next(a.GetIdx() for a in ref.GetAtoms() if a.GetIsAromatic()
                and any(not n.GetIsAromatic() for n in a.GetNeighbors()))
    o1, o2 = [n.GetIdx() for n in ref.GetAtomWithIdx(ipso).GetNeighbors() if n.GetIsAromatic()]
    m1 = next(n.GetIdx() for n in ref.GetAtomWithIdx(o1).GetNeighbors() if n.GetIdx() != ipso)
    m2 = next(n.GetIdx() for n in ref.GetAtomWithIdx(o2).GetNeighbors() if n.GetIdx() != ipso)

    def swapped(pairs):
        probe = Chem.Mol(ref)
        conf = probe.GetConformer()
        for i, j in pairs:
            pi, pj = np.array(conf.GetAtomPosition(i)), np.array(conf.GetAtomPosition(j))
            conf.SetAtomPosition(i, pj.tolist())
            conf.SetAtomPosition(j, pi.tolist())
        return probe

    assert symmetry_rmsd(ref, swapped([(o1, o2), (m1, m2)])) == pytest.approx(0.0, abs=1e-9)
    assert symmetry_rmsd(ref, swapped([(o1, o2)])) > 0.1


def test_cf3_fluorine_permutation_is_free_but_real_displacement_is_not():
    ref = _embed("FC(F)(F)c1ccccc1")
    probe = Chem.Mol(ref)
    fl = [a.GetIdx() for a in probe.GetAtoms() if a.GetAtomicNum() == 9]
    conf = probe.GetConformer()
    p = [np.array(conf.GetAtomPosition(i)) for i in fl]
    for i, pos in zip(fl, [p[1], p[2], p[0]]):
        conf.SetAtomPosition(i, pos.tolist())
    assert symmetry_rmsd(ref, probe) == pytest.approx(0.0, abs=1e-9)
    # Moving one fluorine 2 A should register: sqrt(4 / n_heavy)
    probe2 = Chem.Mol(ref)
    c2 = probe2.GetConformer()
    q = np.array(c2.GetAtomPosition(fl[0])) + np.array([2.0, 0, 0])
    c2.SetAtomPosition(fl[0], q.tolist())
    n_heavy = sum(a.GetAtomicNum() > 1 for a in ref.GetAtoms())
    assert 0 < symmetry_rmsd(ref, probe2) <= np.sqrt(4.0 / n_heavy) + 1e-9


def test_agrees_with_rdkit_calcrms_on_a_noisy_pose():
    ref = _embed("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
    rng = np.random.default_rng(0)
    probe = _move(ref, t=(0.5, -0.3, 0.2))
    conf = probe.GetConformer()
    for i in range(probe.GetNumAtoms()):
        pos = np.array(conf.GetAtomPosition(i)) + rng.normal(scale=0.3, size=3)
        conf.SetAtomPosition(i, pos.tolist())
    ours = symmetry_rmsd(Chem.RemoveHs(ref), Chem.RemoveHs(probe))
    theirs = rdMolAlign.CalcRMS(Chem.RemoveHs(probe), Chem.RemoveHs(ref))
    assert ours == pytest.approx(theirs, abs=1e-6)


def test_hydrogens_and_bond_orders_do_not_matter():
    withh = _embed("c1ccccc1O")
    noh = Chem.RemoveHs(withh)
    assert symmetry_rmsd(withh, noh) == pytest.approx(0.0, abs=1e-9)
    kek = Chem.Mol(noh)
    Chem.Kekulize(kek, clearAromaticFlags=True)
    assert symmetry_rmsd(noh, kek) == pytest.approx(0.0, abs=1e-9)


def test_different_molecules_are_rejected():
    a, b = _embed("CCO"), _embed("CCN")
    with pytest.raises(ValueError, match="not the same molecular graph"):
        symmetry_rmsd(a, b)
    with pytest.raises(ValueError, match="counts differ"):
        symmetry_rmsd(a, _embed("CCCO"))
