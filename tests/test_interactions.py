from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from vsbench import interactions as I
from vsbench import sites as S

ROOT = Path(__file__).parents[1] / "data" / "raw" / "full_data" / "ESR1_ant"


def test_recovery_arithmetic():
    ref = frozenset({("hbond", "A", 1), ("hydrophobic", "A", 2), ("pi_stack", "A", 3), ("hbond", "A", 4)})
    pose = frozenset({("hbond", "A", 1), ("hydrophobic", "A", 2), ("hbond", "A", 9)})
    r = I.recovery(ref, pose)
    assert r["n_shared"] == 2 and r["recall"] == 0.5 and r["precision"] == pytest.approx(2 / 3)
    assert r["tanimoto"] == pytest.approx(2 / 5)
    assert np.isnan(I.recovery(frozenset(), pose)["recall"]) and np.isnan(I.recovery(ref, frozenset())["precision"])
    assert I.residue_level(ref) == frozenset({("A", 1), ("A", 2), ("A", 3), ("A", 4)})


def test_ligand_lines_use_the_fixed_residue_and_chain():
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    from rdkit.Chem import AllChem
    AllChem.EmbedMolecule(mol, randomSeed=1)
    lines = I.ligand_pdb_lines(mol)
    het = [l for l in lines if l.startswith("HETATM")]
    assert len(het) == mol.GetNumAtoms() and all(l[17:20] == "LIG" and l[21] == "Z" and l[22:26].strip() == "1" for l in het)
    assert any(l.startswith("CONECT") for l in lines)


@pytest.mark.skipif(not (ROOT / "1xqc_protein.mol2").exists(), reason="LIT-PCBA data not downloaded")
def test_plip_finds_the_crystal_contacts_and_a_displaced_ligand_loses_them(tmp_path):
    import subprocess
    from vsbench.tools import find_tool
    if find_tool("obabel") is None:
        pytest.skip("obabel missing")
    prot = tmp_path / "p.pdb"
    S.write_pdb_from_mol2(ROOT / "1xqc_protein.mol2", prot)
    sdf = tmp_path / "lig.sdf"
    subprocess.run([str(find_tool("obabel")), str(ROOT / "1xqc_ligand.mol2"), "-O", str(sdf)], capture_output=True)
    lig = Chem.MolFromMolFile(str(sdf), removeHs=False)
    ref = I.interaction_set(I.write_complex(prot, lig, tmp_path / "ref.pdb"))
    kinds = {k for k, _, _ in ref}
    assert len(ref) >= 6 and "hydrophobic" in kinds and "hbond" in kinds             # a real antagonist pocket
    again = I.interaction_set(I.write_complex(prot, lig, tmp_path / "ref2.pdb"))
    assert I.recovery(ref, again)["recall"] == 1.0                                  # deterministic
    moved = Chem.Mol(lig)
    conf = moved.GetConformer()
    for i in range(moved.GetNumAtoms()):
        p = np.array(conf.GetAtomPosition(i)) + np.array([25.0, 0.0, 0.0])
        conf.SetAtomPosition(i, p.tolist())
    far = I.interaction_set(I.write_complex(prot, moved, tmp_path / "far.pdb"))
    assert I.recovery(ref, far)["recall"] < 0.2                                     # moving it away removes the contacts


@pytest.mark.skipif(not (ROOT / "1xqc_protein.mol2").exists(), reason="LIT-PCBA data not downloaded")
def test_plip_leaves_no_temporary_file_in_the_system_temp_directory(tmp_path):
    """PLIP writes a repaired copy of the complex and never deletes it. Left at its default of
    the system temp directory that is one file per pose: a single screen left 10,000 behind,
    2.6 GB, on a tmpfs shared with everything else running."""
    import glob
    import subprocess
    import tempfile
    from vsbench.tools import find_tool
    if find_tool("obabel") is None:
        pytest.skip("obabel missing")
    prot = tmp_path / "p.pdb"
    S.write_pdb_from_mol2(ROOT / "1xqc_protein.mol2", prot)
    sdf = tmp_path / "lig.sdf"
    subprocess.run([str(find_tool("obabel")), str(ROOT / "1xqc_ligand.mol2"), "-O", str(sdf)],
                   capture_output=True)
    lig = Chem.MolFromMolFile(str(sdf), removeHs=False)

    pattern = str(Path(tempfile.gettempdir()) / "plipfixed.*")
    before = set(glob.glob(pattern))
    assert I.interaction_set(I.write_complex(prot, lig, tmp_path / "c.pdb"))
    assert set(glob.glob(pattern)) == before, "PLIP leaked a temporary file into the system temp dir"
    # It is written, just somewhere the caller can clean up.
    assert list(tmp_path.glob("plipfixed.*"))


def test_an_empty_reference_is_refused_rather_than_scored_as_nan(tmp_path, monkeypatch):
    """MTORC1/4jsx produced an empty reference set because the receptor had lost its residue names.
    Nothing raised, so every one of 980 poses was scored NaN and the target looked merely uninformative."""
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    from rdkit.Chem import AllChem
    AllChem.EmbedMolecule(mol, randomSeed=1)
    prot = tmp_path / "p.pdb"
    prot.write_text("ATOM      1  CA  UNK A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n")
    monkeypatch.setattr(I, "_interaction_job", lambda args: ("reference", frozenset(), ""))
    with pytest.raises(RuntimeError, match="no interaction"):
        I.reference_set(prot, mol)


MTOR = Path(__file__).parents[1] / "data" / "raw" / "full_data" / "MTORC1"


@pytest.mark.skipif(not (MTOR / "4jsx_protein.mol2").exists(), reason="LIT-PCBA data not downloaded")
def test_the_real_4jsx_receptor_keeps_mtor_numbering(tmp_path):
    """The structure the defect was found on: mTOR's kinase domain, residues 1385 to 2549 on chain B."""
    prot = tmp_path / "p.pdb"
    S.write_pdb_from_mol2(MTOR / "4jsx_protein.mol2", prot)
    text = prot.read_text()
    assert "UNK" not in text
    res = S.residues_from_pdb_text(text)
    chains = {c for c, _ in res}
    nums = [n for _, n in res]
    assert chains == {"B"} and min(nums) >= 1385 and max(nums) <= 2549
