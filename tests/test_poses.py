import numpy as np
import pytest
from rdkit import Chem

from vsbench import poses as P
from vsbench.ligprep import embed, to_pdbqt
from vsbench.rmsd import symmetry_rmsd

SMILES = [
    "CC(=O)Nc1ccc(O)cc1",
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "FC(F)(F)c1ccc(cc1)N1CCN(CC1)C(=O)c1ccccc1",
    "C[C@H](N)C(=O)[O-]",
]


@pytest.mark.parametrize("smiles", SMILES)
def test_pdbqt_round_trip_reproduces_the_embedded_heavy_atoms(smiles):
    # Ground truth: the conformer that was written into the PDBQT.
    mol = embed(smiles, seed=3)
    pdbqt, _ = to_pdbqt(mol)
    rebuilt = P.pose_to_mol(pdbqt, pdbqt.splitlines())
    assert rebuilt.GetNumAtoms() == mol.GetNumHeavyAtoms()
    assert symmetry_rmsd(Chem.RemoveHs(mol), rebuilt) < 5e-3      # PDBQT keeps 3 decimals
    # The same heavy-atom positions must be recovered as a set, whatever the atom order.
    ref_xyz = np.sort(np.round(Chem.RemoveHs(mol).GetConformer().GetPositions(), 2), axis=0)
    new_xyz = np.sort(np.round(rebuilt.GetConformer().GetPositions(), 2), axis=0)
    assert np.abs(ref_xyz - new_xyz).max() <= 0.011


def test_translated_pose_gives_exact_translation_rmsd():
    mol = embed("CC(=O)Nc1ccc(O)cc1", seed=1)
    pdbqt, _ = to_pdbqt(mol)
    moved = []
    for line in pdbqt.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            line = line[:30] + f"{x + 3.0:8.3f}{y + 4.0:8.3f}{z:8.3f}" + line[54:]
        moved.append(line)
    rebuilt = P.pose_to_mol(pdbqt, moved)
    ref = P.pose_to_mol(pdbqt, pdbqt.splitlines())
    assert symmetry_rmsd(ref, rebuilt) == pytest.approx(5.0, abs=1e-3)


def test_dlg_style_prefixes_are_understood():
    mol = embed("CC(=O)Nc1ccc(O)cc1", seed=1)
    pdbqt, _ = to_pdbqt(mol)
    dlg_like = "\n".join("DOCKED: " + l for l in pdbqt.splitlines())
    header = "\n".join("INPUT-LIGAND-PDBQT: " + l for l in pdbqt.splitlines() if l.startswith("REMARK"))
    smiles, mapping = P.parse_mapping(header)
    assert smiles == "CC(=O)Nc1ccc(O)cc1"
    assert len(mapping) == 11
    assert P.pose_to_mol(dlg_like, dlg_like.splitlines()).GetNumAtoms() == 11


def test_missing_or_inconsistent_inputs_raise():
    with pytest.raises(ValueError, match="no usable"):
        P.parse_mapping("ATOM 1 C")
    mol = embed("CCO", seed=1)
    pdbqt, _ = to_pdbqt(mol)
    short = [l for l in pdbqt.splitlines() if not l.startswith("ATOM      3")]
    with pytest.raises(ValueError, match="lacks atom serial"):
        P.pose_to_mol(pdbqt, short)


def test_split_models_handles_pdbqt_and_dlg_forms():
    text = "MODEL 1\nATOM a\nENDMDL\nMODEL 2\nATOM b\nATOM c\nENDMDL\n"
    assert P.split_models(text) == [["ATOM a"], ["ATOM b", "ATOM c"]]
    dlg = "DOCKED: MODEL        1\nDOCKED: ATOM a\nDOCKED: ENDMDL\n"
    assert P.split_models(dlg) == [["DOCKED: ATOM a"]]
