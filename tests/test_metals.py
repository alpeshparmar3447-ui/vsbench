import numpy as np
import pytest
from rdkit import Chem

from vsbench import metals as K

PDB = """\
ATOM      1  N   HIS A  94      10.000  10.000  10.000  1.00  0.00           N
ATOM      2  CA  HIS A  94      11.000  10.000  10.000  1.00  0.00           C
ATOM      3  NE2 HIS A  94       9.000  12.000  10.000  1.00  0.00           N
HETATM  100 ZN    ZN A 301       9.000  14.000  10.000  1.00  0.00          ZN
HETATM  101 MG    MG A 302      50.000  50.000  50.000  1.00  0.00          MG
HETATM  102  K     K A 303      12.000  14.000  10.000  1.00  0.00           K
HETATM  103  O   HOH A 400       9.500  14.500  10.500  1.00  0.00           O
HETATM  104 CA    CA A 305      80.000  80.000  80.000  0.50  0.00          CA
HETATM  105 CA    CA B 305      81.000  80.000  80.000  0.50  0.00          CA
HETATM  106  C1  LIG A 500      12.000  16.000  10.000  1.00  0.00           C
"""

MOL2 = """\
@<TRIPOS>MOLECULE
x
@<TRIPOS>ATOM
      1 N          9.0000   12.0000   10.0000 N.am      1 HIS94  0.0
      2 CA        11.0000   10.0000   10.0000 C.3       1 HIS94  0.0
      3 ZN         9.0000   14.0000   10.0000 Zn        2 ZN301  0.0
      4 CA        80.0000   80.0000   80.0000 Ca        3 CA305  0.0
      5 CA        20.0000   20.0000   20.0000 C.3       4 ALA5   0.0
@<TRIPOS>BOND
"""


def test_pdb_ions_are_found_and_waters_and_ligands_are_ignored():
    ions = K.read_pdb_ions(PDB)
    assert [i.element for i in ions] == ["ZN", "MG", "K", "CA", "CA"]
    assert ions[0].xyz == (9.0, 14.0, 10.0) and ions[0].chain == "A" and ions[0].resnum == 301


def test_alternate_locations_keep_only_the_first_conformer():
    text = PDB + "HETATM  107 ZN   BZN A 301      30.000  30.000  30.000  0.50  0.00          ZN\n"
    assert len(K.read_pdb_ions(text)) == 5          # the 'B' location is not added


def test_mol2_ions_do_not_mistake_alpha_carbons_for_calcium():
    ions = K.read_mol2_ions(MOL2)
    assert [(i.element, i.xyz[0]) for i in ions] == [("ZN", 9.0), ("CA", 80.0)]   # atom named CA, type C.3, is not Ca


def test_cif_ions_and_ca_coordinates(tmp_path):
    cif = """data_x
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA . ALA A 1 1.0 2.0 3.0 1
HETATM 2 ZN ZN . ZN A 301 4.0 5.0 6.0 1
HETATM 3 O O . HOH A 400 7.0 8.0 9.0 1
#
"""
    p = tmp_path / "x.cif"
    p.write_text(cif)
    ions = K.read_structure_ions(p)
    assert len(ions) == 1 and ions[0].element == "ZN" and ions[0].xyz == (4.0, 5.0, 6.0)
    assert K.read_ca_coordinates(p).tolist() == [[1.0, 2.0, 3.0]]


def test_ca_coordinates_from_pdb_and_mol2(tmp_path):
    (tmp_path / "a.pdb").write_text(PDB)
    (tmp_path / "a.mol2").write_text(MOL2)
    assert K.read_ca_coordinates(tmp_path / "a.pdb").tolist() == [[11.0, 10.0, 10.0]]
    assert sorted(K.read_ca_coordinates(tmp_path / "a.mol2").tolist()) == [[11.0, 10.0, 10.0], [20.0, 20.0, 20.0]]


def test_frame_overlap_detects_a_shared_and_a_shifted_frame():
    a = np.random.default_rng(0).uniform(0, 50, (100, 3))
    assert K.frame_overlap(a, a + 0.1) == 1.0
    assert K.frame_overlap(a, a + 5.0) < 0.1
    assert K.frame_overlap(np.empty((0, 3)), a) == 0.0


def test_pocket_selection_keeps_functional_ions_near_the_ligand_only():
    ions = K.read_pdb_ions(PDB)
    lig = np.array([[9.0, 17.0, 10.0], [10.0, 18.0, 10.0]])
    kept = K.relevant_ions(ions, lig)
    assert [(i.element, round(d, 1)) for i, d in kept] == [("ZN", 3.0), ("K", 4.2)]   # Mg and Ca are remote
    assert K.relevant_ions(ions, lig, monovalent_cutoff=2.0)[0][0].element == "ZN"
    assert len(K.relevant_ions(ions, lig, monovalent_cutoff=2.0)) == 1                # K excluded when far
    far = K.ions_near(ions, lig, 200.0)
    assert [h[0].element for h in far][:2] == ["ZN", "K"] and far[0][1] <= far[1][1]


def test_coordination_counts_protein_donors_only():
    atoms = list(K.read_pdb_protein_atoms(PDB))
    zn = K.read_pdb_ions(PDB)[0]
    shell = K.coordination(zn.xyz, atoms, cutoff=3.0)
    assert [(a, r) for a, r, _, _ in shell] == [("NE2", "HIS")]          # protein donors; waters and ions excluded
    assert shell[0][3] == pytest.approx(2.0)
    assert K.coordination(zn.xyz, atoms, cutoff=1.0) == []
    with_water = K.coordination(zn.xyz, list(K.read_pdb_protein_atoms(PDB, include_hetero=True)), cutoff=3.0)
    assert [(a, r) for a, r, _, _ in with_water] == [("O", "HOH"), ("NE2", "HIS")]   # opt-in, nearest first


def test_ion_pdb_line_round_trips_through_the_reader():
    ion = K.read_pdb_ions(PDB)[0]
    line = K.ion_to_pdb_line(ion, 5)
    assert len(line) >= 78 and line.startswith("HETATM")
    back = K.read_pdb_ions(line)
    assert len(back) == 1 and back[0].element == "ZN" and back[0].xyz == ion.xyz and back[0].resnum == 301
    mg = K.ion_to_pdb_line(K.read_pdb_ions(PDB)[1], 6)
    assert K.read_pdb_ions(mg)[0].element == "MG"


@pytest.mark.parametrize("smiles, auto", [
    ("CC(=O)NO", ["hydroxamic_acid"]),
    ("Oc1ccccc1O", ["catechol"]),
    ("O=c1ccn(C)cc1O", ["hydroxypyridinone_or_tropolone"]),
    ("CCO", []),
])
def test_chelator_detection(smiles, auto):
    assert K.chelator_classes(Chem.MolFromSmiles(smiles))["auto"] == auto


@pytest.mark.parametrize("smiles, manual", [
    ("CCS", ["thiol"]), ("CP(=O)(O)O", ["phosphonate_or_phosphate_OH"]),
    ("Oc1cccc2cccnc12", ["8_hydroxyquinoline"]), ("NS(=O)(=O)c1ccccc1", ["primary_sulfonamide"]),
])
def test_groups_needing_manual_ionisation_are_flagged_not_changed(smiles, manual):
    mol = Chem.MolFromSmiles(smiles)
    assert K.chelator_classes(mol)["manual"] == manual
    same, classes = K.deprotonate_chelators(mol)
    assert classes == [] and Chem.MolToSmiles(same) == Chem.MolToSmiles(mol)


def test_deprotonation_gives_the_chelator_anion_and_keeps_the_molecule_valid():
    mol = Chem.MolFromSmiles("CC(=O)NO")
    out, classes = K.deprotonate_chelators(mol)
    assert classes == ["hydroxamic_acid"]
    assert Chem.GetFormalCharge(out) == -1
    assert Chem.MolToSmiles(out) == "CC(=O)N[O-]"
    cat, c2 = K.deprotonate_chelators(Chem.MolFromSmiles("Oc1ccccc1O"))
    assert c2 == ["catechol"] and Chem.GetFormalCharge(cat) == -1            # one oxygen, as in the tested module
    none, c3 = K.deprotonate_chelators(Chem.MolFromSmiles("CCO"))
    assert c3 == [] and Chem.GetFormalCharge(none) == 0


# ------------------------------------------------------------------ planning
from pathlib import Path  # noqa: E402

from vsbench.box import Box  # noqa: E402


def _pdb(tmp_path, name, ions=(), shift=(0, 0, 0)):
    """Tiny PDB: a few C-alpha atoms plus the given (element, x, y, z) ions."""
    lines = []
    for n in range(40):
        x, y, z = 10.0 + 3.0 * (n % 7) + shift[0], 5.0 + 3.1 * (n // 7) + shift[1], 2.0 + 1.7 * (n % 5) + shift[2]
        lines.append(f"ATOM  {n + 1:5d}  CA  ALA A{n + 1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C")
    for k, (el, x, y, z) in enumerate(ions):
        lines.append(f"HETATM{900 + k:5d} {el:>2s}   {el:>2s} A{300 + k:4d}    {x + shift[0]:8.3f}{y + shift[1]:8.3f}{z + shift[2]:8.3f}  1.00  0.00          {el:>2s}")
    p = tmp_path / name
    p.write_text("\n".join(lines) + "\n")
    return p


LIG = np.array([[20.0, 10.0, 5.0], [21.0, 10.5, 5.5], [22.0, 11.0, 6.0]])
BOX = Box(center=(21.0, 10.5, 5.5), size=(24.0, 24.0, 24.0))


def _fetch_returning(path):
    return lambda pdb_id, cache_dir, timeout=60.0: path


def test_plan_restores_a_pocket_ion_that_the_benchmark_file_lacks(tmp_path):
    bench = _pdb(tmp_path, "bench.pdb")
    dep = _pdb(tmp_path, "dep.pdb", ions=[("ZN", 20.0, 12.0, 5.0)])
    plan = K.plan_metals(bench, LIG, BOX, "1abc", tmp_path, fetch=_fetch_returning(dep))
    assert plan.source == "restored" and plan.audited and plan.zn
    assert [i.element for i in plan.ions] == ["ZN"] and plan.n_deposition_ions == 1
    assert "retained: ZN" in plan.note and plan.label.startswith("ZN")


def test_plan_does_not_duplicate_an_ion_present_in_both_files(tmp_path):
    bench = _pdb(tmp_path, "bench.pdb", ions=[("MG", 20.0, 12.0, 5.0)])
    dep = _pdb(tmp_path, "dep.pdb", ions=[("MG", 20.0, 12.0, 5.0)])
    plan = K.plan_metals(bench, LIG, BOX, "1abc", tmp_path, fetch=_fetch_returning(dep))
    assert len(plan.ions) == 1 and plan.source == "benchmark_file" and not plan.zn


def test_plan_reports_but_keeps_out_a_remote_structural_ion(tmp_path):
    bench = _pdb(tmp_path, "bench.pdb")
    dep = _pdb(tmp_path, "dep.pdb", ions=[("ZN", 50.0, 10.5, 5.5)])            # about 29 A away, outside the box
    plan = K.plan_metals(bench, LIG, BOX, "1abc", tmp_path, fetch=_fetch_returning(dep))
    assert plan.ions == () and plan.source == "none" and not plan.zn
    assert len(plan.remote) == 1 and plan.remote[0].startswith("ZN") and "kept out" in plan.note


def test_an_ion_inside_the_search_box_is_kept_even_beyond_the_pocket_cutoff(tmp_path):
    bench = _pdb(tmp_path, "bench.pdb")
    dep = _pdb(tmp_path, "dep.pdb", ions=[("MG", 30.5, 10.5, 5.5)])           # 8.5 A from the ligand, inside the 24 A box
    near = K.plan_metals(bench, LIG, BOX, "1abc", tmp_path, fetch=_fetch_returning(dep))
    assert [i.element for i in near.ions] == ["MG"]
    no_box = K.plan_metals(bench, LIG, Box((100.0, 100.0, 100.0), (1.0, 1.0, 1.0)), "1abc", tmp_path,
                           fetch=_fetch_returning(dep))
    assert no_box.ions == ()                                                  # beyond 8 A and outside a tiny box


def test_plan_refuses_a_deposition_in_a_different_frame(tmp_path):
    bench = _pdb(tmp_path, "bench.pdb")
    dep = _pdb(tmp_path, "dep.pdb", ions=[("ZN", 20.0, 12.0, 5.0)], shift=(40.0, 0, 0))
    plan = K.plan_metals(bench, LIG, BOX, "1abc", tmp_path, fetch=_fetch_returning(dep))
    assert plan.ions == () and not plan.audited and "frame overlap" in plan.note


def test_plan_survives_an_unavailable_deposition_and_says_so(tmp_path):
    bench = _pdb(tmp_path, "bench.pdb", ions=[("ZN", 20.0, 12.0, 5.0)])

    def boom(pdb_id, cache_dir, timeout=60.0):
        raise OSError("no network")

    plan = K.plan_metals(bench, LIG, BOX, "1abc", tmp_path, fetch=boom)
    assert [i.element for i in plan.ions] == ["ZN"] and not plan.audited
    assert "unavailable (OSError)" in plan.note


def test_plan_off_mode_and_no_pdb_code(tmp_path):
    bench = _pdb(tmp_path, "bench.pdb", ions=[("ZN", 20.0, 12.0, 5.0)])
    off = K.plan_metals(bench, LIG, BOX, "1abc", tmp_path, mode="off")
    assert off.ions == () and "disabled" in off.note
    plain = K.plan_metals(bench, LIG, BOX, None, None)
    assert [i.element for i in plain.ions] == ["ZN"] and plain.audited is False


def test_plan_with_nothing_near_the_ligand_says_the_deposition_was_checked(tmp_path):
    bench = _pdb(tmp_path, "bench.pdb")
    dep = _pdb(tmp_path, "dep.pdb")
    plan = K.plan_metals(bench, LIG, BOX, "1abc", tmp_path, fetch=_fetch_returning(dep))
    assert plan.ions == () and plan.audited and "deposition checked" in plan.note
