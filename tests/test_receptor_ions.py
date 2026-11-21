"""Receptor-ion logic that needs no docking executables."""
import numpy as np
import pytest

from vsbench import engines as E
from vsbench import metals as K

PDB = """\
ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM      2  CA  ALA A   1      11.000  10.000  10.000  1.00  0.00           C
ATOM      3  C   ALA A   1      12.000  10.500  10.000  1.00  0.00           C
ATOM      4  CA BALA A   1      13.000  10.500  10.000  1.00  0.00           C
TER
HETATM  100 ZN    ZN A 301       9.000  14.000  10.000  1.00  0.00          ZN
HETATM  101 MG    MG A 302      50.000  50.000  50.000  1.00  0.00          MG
HETATM  102  O   HOH A 400       9.500  14.500  10.500  1.00  0.00           O
HETATM  103  C1  LIG A 500      12.000  16.000  10.000  1.00  0.00           C
HETATM  104  C2  LIG A 500      12.500  16.500  10.000  1.00  0.00           C
"""


def _ions():
    return K.read_pdb_ions(PDB)


def test_source_to_pdb_keeps_only_requested_ions_and_reports_everything_removed(tmp_path):
    src = tmp_path / "in.pdb"
    src.write_text(PDB)
    zn = _ions()[0]
    rep = E.source_to_pdb(src, tmp_path / "out.pdb", keep_ions=[zn])
    out = (tmp_path / "out.pdb").read_text().splitlines()
    assert rep.ions_kept == (zn.key(),) and rep.n_protein_atoms == 3
    assert rep.waters_removed == 1
    assert rep.hetero_removed == {"MG": 1, "LIG": 2}          # the kept Zn is not reported as removed
    het = [l for l in out if l.startswith("HETATM")]
    assert len(het) == 1 and het[0][17:20].strip() == "ZN"
    assert [l for l in out if l.startswith("ATOM")] and not any("HOH" in l or "LIG" in l for l in out)


def test_only_the_first_alternate_location_is_kept(tmp_path):
    src = tmp_path / "in.pdb"
    src.write_text(PDB)                                        # atom 4 is altloc B
    assert E.source_to_pdb(src, tmp_path / "b.pdb").n_protein_atoms == 3
    src.write_text(PDB.replace("CA BALA", "CA AALA"))          # altloc A is kept
    assert E.source_to_pdb(src, tmp_path / "a.pdb").n_protein_atoms == 4


def test_no_ions_requested_means_an_apo_receptor_and_the_report_says_so(tmp_path):
    src = tmp_path / "in.pdb"
    src.write_text(PDB)
    rep = E.source_to_pdb(src, tmp_path / "out.pdb")
    assert rep.ions_kept == () and rep.hetero_removed["ZN"] == 1 and rep.hetero_removed["MG"] == 1
    assert "HETATM" not in (tmp_path / "out.pdb").read_text()


def test_ion_survival_check_raises_when_an_ion_is_missing(tmp_path):
    pdbqt = tmp_path / "r.pdbqt"
    pdbqt.write_text("ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00    -0.350 N \n")
    with pytest.raises(RuntimeError, match="missing from the prepared receptor"):
        E._verify_ions_in_pdbqt(pdbqt, [_ions()[0]])
    pdbqt.write_text(pdbqt.read_text() + "ATOM      2  ZN  ZN  A 301       9.000  14.000  10.000  1.00  0.00     2.000 Zn\n")
    E._verify_ions_in_pdbqt(pdbqt, [_ions()[0]])                  # now present: no error
    with pytest.raises(RuntimeError):
        E._verify_ions_in_pdbqt(pdbqt, [_ions()[0], _ions()[1]])  # Mg still missing


def test_atom_type_casing_is_restored_after_the_zinc_script_upper_cases_it(tmp_path):
    original = tmp_path / "orig.pdbqt"
    corrupted = tmp_path / "bad.pdbqt"
    original.write_text(
        "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00    -0.350 N \n"
        "ATOM      2  ZN  ZN  A 301       9.000  14.000  10.000  1.00  0.00     2.000 Zn\n"
        "ATOM      3  MG  MG  A 302      50.000  50.000  50.000  1.00  0.00     2.000 Mg\n")
    corrupted.write_text(
        "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00    -0.350 N \n"
        "ATOM      2  ZN  ZN  A 301       9.000  14.000  10.000  1.00  0.00     0.000 ZN\n"
        "ATOM      3  MG  MG  A 302      50.000  50.000  50.000  1.00  0.00     2.000 MG\n"
        "ATOM      4   TZ ZN  A 301      10.900  14.300  10.300  1.00  0.00     0.000 TZ\n")
    E._restore_ion_type_case(original, corrupted)
    types = [l[77:79].strip() for l in corrupted.read_text().splitlines()]
    assert types == ["N", "Zn", "Mg", "TZ"]                       # casing fixed, pseudo-atom untouched


def test_zn_mode_without_a_zinc_ion_is_rejected_before_any_work(tmp_path):
    from vsbench.box import Box
    src = tmp_path / "in.pdb"
    src.write_text(PDB)
    with pytest.raises(ValueError, match="no zinc ion"):
        E.prepare_receptor(src, tmp_path / "w", Box((0, 0, 0), (24, 24, 24)), keep_ions=[_ions()[1]], zn_mode=True)
