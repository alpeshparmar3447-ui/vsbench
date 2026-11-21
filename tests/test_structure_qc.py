import math

import numpy as np
import pytest

from vsbench import structure_qc as Q


def test_dihedral_known_geometries():
    p = [(1, 0, 0), (0, 0, 0), (0, 0, 1), (0, 1, 1)]
    assert Q.dihedral(*p) == pytest.approx(90.0)
    assert Q.dihedral((1, 0, 0), (0, 0, 0), (0, 0, 1), (1, 0, 1)) == pytest.approx(0.0)
    assert abs(Q.dihedral((1, 0, 0), (0, 0, 0), (0, 0, 1), (-1, 0, 1))) == pytest.approx(180.0)
    assert Q.dihedral((1, 0, 0), (0, 0, 0), (0, 0, 1), (0, -1, 1)) == pytest.approx(-90.0)


@pytest.mark.parametrize("phi, psi, expected", [
    (-60, -45, "core"), (-120, 130, "core"), (55, 35, "core"),
    (-95, 150, "core"), (-60, -10, "allowed"), (60, -150, "outside"), (0, 0, "outside"),
    (-180, 170, "core"),                        # periodic wrap: phi near +-180, psi near the beta region
])
def test_ramachandran_regions(phi, psi, expected):
    assert Q.classify_rama(phi, psi) in {expected, "allowed"} if expected == "core" and (phi, psi) == (-180, 170) \
        else Q.classify_rama(phi, psi) == expected


def _mol2(residues, breaks=()):
    """Straight poly-peptide-like mol2: each residue N, CA, C, O spaced so C-N is 1.33 A."""
    lines = ["@<TRIPOS>MOLECULE", "x", "@<TRIPOS>ATOM"]
    serial, x = 1, 0.0
    subs = []
    for i, (chain, num, name) in enumerate(residues, 1):
        if (i - 1) in breaks:
            x += 6.0                                       # a break: next residue 6 A away
        for aname, dx, typ in (("N", 0.0, "N.am"), ("CA", 1.4, "C.3"), ("C", 2.9, "C.2"), ("O", 2.9, "O.2")):
            lines.append(f"{serial:7d} {aname:<4s} {x + dx:9.4f} {0.3 * (aname == 'O'):9.4f} 0.0000 {typ:5s} {i:3d} {name}{num} 0.0")
            serial += 1
        x += 4.23
        subs.append(f"{i:6d} {name}{num}      1 RESIDUE           4 {chain}     {name}     1 ROOT {name} {chain} {num}")
    lines += ["@<TRIPOS>BOND", "@<TRIPOS>SUBSTRUCTURE"] + subs
    return "\n".join(lines) + "\n"


def test_residue_reader_takes_chain_and_original_numbers(tmp_path):
    p = tmp_path / "a.mol2"
    p.write_text(_mol2([("A", 307, "ALA"), ("A", 308, "LEU"), ("B", 10, "GLY")]))
    res = Q.read_mol2_residues(p)
    assert [(r.chain, r.number, r.name) for r in res] == [("A", 307, "ALA"), ("A", 308, "LEU"), ("B", 10, "GLY")]
    assert set(res[0].atoms) == {"N", "CA", "C", "O"}


def test_breaks_gaps_and_missing_atoms_are_counted_separately(tmp_path):
    seq = [("A", 1, "ALA"), ("A", 2, "ALA"), ("A", 5, "ALA"), ("A", 6, "ALA"), ("A", 7, "ALA")]
    p = tmp_path / "a.mol2"
    p.write_text(_mol2(seq, breaks={2}))                    # physical break between residues 2 and 5
    r = Q.backbone_report(Q.read_mol2_residues(p))
    assert r["n_chain_breaks"] == 1 and r["n_numbering_gaps"] == 1 and r["n_unexplained_breaks"] == 0
    p.write_text(_mol2(seq, breaks={3}))                    # break where numbering is continuous
    r = Q.backbone_report(Q.read_mol2_residues(p))
    assert r["n_unexplained_breaks"] == 1 and r["n_numbering_gaps"] == 1
    text = _mol2(seq).replace("O    ", "OXT  ", 1)         # rename one carbonyl O: backbone atom missing
    p.write_text(text)
    assert Q.backbone_report(Q.read_mol2_residues(p))["n_missing_backbone"] == 1


def test_two_chains_are_reported(tmp_path):
    p = tmp_path / "a.mol2"
    p.write_text(_mol2([("A", 1, "ALA"), ("A", 2, "ALA"), ("B", 1, "ALA"), ("B", 2, "ALA")]))
    assert Q.backbone_report(Q.read_mol2_residues(p))["n_chains"] == 2


PDB = """\
REMARK   2 RESOLUTION.    1.95 ANGSTROMS.
SEQADV 1ABC GLY A  300  UNP  P12345              EXPRESSION TAG
SEQADV 1ABC GLY A  301  UNP  P12345              EXPRESSION TAG
SEQADV 1ABC ALA A  120  UNP  P12345    VAL   120 ENGINEERED MUTATION
ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM      2  CA AALA A   1      11.000  10.000  10.000  0.50  0.00           C
ATOM      3  CA BALA A   1      11.100  10.000  10.000  0.50  0.00           C
HETATM  100  C1  LIG A 500      12.000  16.000  10.000  1.00  0.00           C
HETATM  101  C2  LIG A 500      12.500  16.500  10.000  1.00  0.00           C
HETATM  110  P   ATP A 600      12.000  19.000  10.000  1.00  0.00           P
HETATM  111  O1G ATP A 600      12.500  19.500  10.000  1.00  0.00           O
HETATM  120  C   GOL A 700      60.000  60.000  60.000  1.00  0.00           C
HETATM  130  O   HOH A 800      12.000  17.000  10.000  1.00  0.00           O
HETATM  140 ZN    ZN A 900      12.000  18.000  10.000  1.00  0.00          ZN
"""


def test_deposition_facts(tmp_path):
    p = tmp_path / "1abc.pdb"
    p.write_text(PDB)
    assert Q.parse_resolution(PDB, False) == pytest.approx(1.95)
    assert Q.count_engineered(PDB, False) == {"expression_tag": 2, "engineered_mutation": 1, "other_seqadv": 0}
    assert Q.count_altlocs(PDB, False) == 2
    lig = np.array([[12.0, 16.0, 10.0], [12.5, 16.5, 10.0]])
    assert Q.hetero_near_ligand(p, lig) == ["ATP"]            # the ligand itself, water, ions and the far GOL are excluded
    assert math.isnan(Q.parse_resolution("no remark", False))


def test_template_qc_flags(tmp_path):
    prot = tmp_path / "p.mol2"
    prot.write_text(_mol2([("A", 1, "ALA"), ("A", 2, "ALA"), ("A", 3, "ALA"), ("A", 4, "ALA")], breaks={2}))
    dep = tmp_path / "d.pdb"
    dep.write_text(PDB)
    qc = Q.qc_template("1abc", prot, np.array([[12.0, 16.0, 10.0], [12.5, 16.5, 10.0]]), dep)
    assert qc.resolution == pytest.approx(1.95) and qc.hetero_near_ligand == "ATP"
    flags = set(qc.flags.split(";"))
    assert {"unexplained_chain_break", "hetero_group_near_ligand", "engineered_mutation"} <= flags
    clean = Q.qc_template("x", tmp_path / "p.mol2", None, None)
    assert "hetero_group_near_ligand" not in clean.flags and math.isnan(clean.resolution)


def test_rama_outlier_flagging_is_relative_to_the_set():
    fr = {f"t{i}": 0.02 + 0.001 * i for i in range(10)}
    fr["bad"] = 0.20
    assert Q.flag_rama_outliers(fr) == ["bad"]
    assert Q.flag_rama_outliers({"a": 0.1, "b": 0.2}) == []      # too few templates to judge


def test_a_missing_terminal_oxygen_is_not_a_missing_backbone_atom(tmp_path):
    seq = [("A", 1, "ALA"), ("A", 2, "ALA"), ("A", 3, "ALA")]
    text = _mol2(seq)
    lines = text.splitlines()
    last_o = max(i for i, l in enumerate(lines) if l.split()[1:2] == ["O"])
    lines[last_o] = lines[last_o].replace(" O ", " OXT ", 1)           # terminal residue's carbonyl renamed
    p = tmp_path / "t.mol2"
    p.write_text("\n".join(lines) + "\n")
    assert Q.backbone_report(Q.read_mol2_residues(p))["n_missing_backbone"] == 0


def _res(chain, num, name, joined=True):
    import numpy as np
    from vsbench.structure_qc import Residue
    atoms = {"N": np.array([0.0, 0.0, float(num) * 3.6]), "CA": np.array([1.4, 0.0, float(num) * 3.6]),
             "C": np.array([2.4, 0.0, float(num) * 3.6]), "O": np.array([2.9, 0.0, float(num) * 3.6])}
    return Residue(chain, num, name, atoms)


def test_waters_are_not_counted_as_a_chain_or_a_numbering_gap():
    """A crystallographic water numbered in the thousands sat at the end of the protein chain. Counting
    it as polymer turned it into a numbering gap; losing its name turned it into a second chain."""
    import vsbench.structure_qc as Q
    prot = [_res("A", n, "ALA") for n in (301, 302, 303)]
    rep = Q.backbone_report(prot + [_res("A", 1003, "HOH")])
    assert rep["n_chains"] == 1 and rep["n_residues"] == 3 and rep["n_numbering_gaps"] == 0
    assert Q.backbone_report(prot + [_res("A", 1004, "ZN")])["n_residues"] == 3
    # A real internal gap is still reported.
    assert Q.backbone_report([_res("A", 301, "ALA"), _res("A", 305, "ALA")])["n_numbering_gaps"] == 1


def test_a_modified_residue_is_kept_in_the_backbone_report():
    import vsbench.structure_qc as Q
    rep = Q.backbone_report([_res("A", 301, "ALA"), _res("A", 302, "CSO"), _res("A", 303, "ALA")])
    assert rep["n_residues"] == 3                                 # CSO is polymer, not a heterogen
