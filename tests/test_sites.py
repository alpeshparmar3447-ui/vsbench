import json
from pathlib import Path

import numpy as np
import pytest

from vsbench import sites as S
from vsbench.structure_qc import read_mol2_residues

FIX = Path(__file__).parent / "fixtures"


def _mol2(n_res=6):
    lines = ["@<TRIPOS>MOLECULE", "x", "@<TRIPOS>ATOM"]
    serial, subs = 1, []
    for i in range(1, n_res + 1):
        x = 4.0 * i
        for name, dx, typ in (("N", 0, "N.am"), ("CA", 1.4, "C.3"), ("C", 2.9, "C.2"), ("O", 2.9, "O.2"), ("H", 0.2, "H")):
            lines.append(f"{serial:7d} {name:<4s} {x + dx:9.4f} 0.0000 0.0000 {typ:5s} {i:3d} ALA{300 + i} 0.0")
            serial += 1
        subs.append(f"{i:6d} ALA{300 + i}      1 RESIDUE           4 A     ALA     1 ROOT ALA A {300 + i}")
    return "\n".join(lines + ["@<TRIPOS>BOND", "@<TRIPOS>SUBSTRUCTURE"] + subs) + "\n"


def _mol2_four_digit(n_res=6, chain="B", first=1385):
    """The mol2 layout of a protein whose residue numbers reach four digits: the SUBSTRUCTURE line ends
    'GLU B1385', one field shorter than the usual 'ALA A 301', because the chain and number are joined."""
    lines = ["@<TRIPOS>MOLECULE", "x", "@<TRIPOS>ATOM"]
    serial, subs = 1, []
    for i in range(1, n_res + 1):
        x, num = 4.0 * i, first + i - 1
        for name, dx, typ in (("N", 0, "N.am"), ("CA", 1.4, "C.3"), ("C", 2.9, "C.2"), ("O", 2.9, "O.2")):
            lines.append(f"{serial:7d} {name:<4s} {x + dx:9.4f} 0.0000 0.0000 {typ:5s} {i:3d} GLU{num} 0.0")
            serial += 1
        subs.append(f"{i:6d} GLU{num}      1 RESIDUE           4 {chain}     GLU     1 ROOT GLU {chain}{num}")
    return "\n".join(lines + ["@<TRIPOS>BOND", "@<TRIPOS>SUBSTRUCTURE"] + subs) + "\n"


def test_a_four_digit_residue_number_keeps_its_chain_and_identity(tmp_path):
    """The 4jsx/4jt5 defect: a joined chain and residue number made the whole SUBSTRUCTURE section
    unparsable, so every residue was written UNK and renumbered from 1. Pocket detectors then returned
    about a hundred spurious cavities and PLIP found no interaction at all, neither raising an error."""
    p = tmp_path / "mtor.mol2"
    p.write_text(_mol2_four_digit())
    out = tmp_path / "mtor.pdb"
    n = S.write_pdb_from_mol2(p, out)
    text = out.read_text()
    assert n == 24 and "UNK" not in text
    assert S.residues_from_pdb_text(text) == frozenset(("B", 1385 + i) for i in range(6))
    assert all(l[17:20] == "GLU" for l in text.splitlines() if l.startswith("ATOM"))


def test_substructure_residue_reads_both_layouts():
    usual = "1 ALA301      1 RESIDUE           4 A     ALA     1 ROOT ALA A 301".split()
    joined = "1 GLU1385     1 RESIDUE           4 B     GLU     1 ROOT GLU B1385".split()
    assert S.substructure_residue(usual) == ("A", 301, "ALA")
    assert S.substructure_residue(joined) == ("B", 1385, "GLU")
    assert S.substructure_residue("1 X 1".split()) is None
    # An unusual residue is marked '****' and named in the following field, in either layout.
    star_usual = "1 CSO303      1 RESIDUE           4 A     ****     1 ROOT CSO A 303".split()
    star_joined = "1 CSO1303     1 RESIDUE           4 B     ****     1 ROOT CSO B1303".split()
    assert S.substructure_residue(star_usual) == ("A", 303, "CSO")
    assert S.substructure_residue(star_joined) == ("B", 1303, "CSO")


def test_a_mol2_with_no_readable_residue_is_refused_rather_than_written_as_unk(tmp_path):
    """Losing residue identity is silent in every downstream tool, so it has to fail here."""
    p = tmp_path / "broken.mol2"
    p.write_text(_mol2().split("@<TRIPOS>SUBSTRUCTURE")[0] + "@<TRIPOS>SUBSTRUCTURE\n")
    with pytest.raises(ValueError, match="residue"):
        S.write_pdb_from_mol2(p, tmp_path / "broken.pdb")


def test_native_site_selects_residues_within_the_cutoff(tmp_path):
    p = tmp_path / "a.mol2"
    p.write_text(_mol2())
    res = read_mol2_residues(p)
    lig = np.array([[8.0, 3.0, 0.0]])                       # 3.2 A from residue 1 (C) and 3.0 A from residue 2 (N)
    site = S.native_site(res, lig, cutoff=4.0)
    assert site == frozenset({("A", 301), ("A", 302)})
    assert S.native_site(res, lig, cutoff=0.5) == frozenset()


def test_pdb_written_from_mol2_keeps_original_numbering_and_drops_hydrogens(tmp_path):
    p = tmp_path / "a.mol2"
    p.write_text(_mol2())
    out = tmp_path / "a.pdb"
    n = S.write_pdb_from_mol2(p, out)
    text = out.read_text()
    assert n == 24 and text.count("ATOM") == 24 and " H " not in text.replace("HETATM", "")
    assert S.residues_from_pdb_text(text) == frozenset(("A", 300 + i) for i in range(1, 7))
    assert text.splitlines()[0][12:16] == " N  " and text.splitlines()[1][12:16] == " CA "
    coords = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in text.splitlines() if l.startswith("ATOM")])
    assert coords[0].tolist() == [4.0, 0.0, 0.0]


def test_recovery_and_ranked_site_recovery():
    native = frozenset({("A", 1), ("A", 2), ("A", 3), ("A", 4)})
    pockets = [S.Pocket(1, 0.9, 0.8, 500, frozenset({("A", 9)})),
               S.Pocket(2, 0.7, 0.4, 300, frozenset({("A", 1), ("A", 2), ("A", 3), ("A", 8)})),
               S.Pocket(3, 0.5, 0.2, 200, frozenset({("A", 4)}))]
    assert S.recovery(native, pockets[1].residues) == 0.75
    r = S.site_recovery(pockets, native)
    assert r["rank_of_site"] == 2 and not r["top1_hit"] and r["top3_hit"]
    assert r["best_recovery"] == 0.75 and r["site_pocket_druggability"] == 0.4 and r["site_pocket_volume"] == 300
    miss = S.site_recovery([pockets[0]], native)
    assert np.isnan(miss["rank_of_site"]) and not miss["top3_hit"] and np.isnan(miss["site_pocket_druggability"])
    assert np.isnan(S.site_recovery([], native)["best_recovery"])


def test_site_recovery_orders_by_rank_not_by_list_position():
    native = frozenset({("A", 1)})
    pockets = [S.Pocket(3, 0, 0, 0, frozenset({("A", 1)})), S.Pocket(1, 0, 0, 0, frozenset()), S.Pocket(2, 0, 0, 0, frozenset())]
    assert S.site_recovery(pockets, native)["rank_of_site"] == 3


# ----------------------------------------------------------------- real fpocket output
def test_fpocket_info_parses_real_output():
    pk = S.parse_fpocket_info((FIX / "fpocket_info.txt").read_text())
    assert [p["pocket"] for p in pk] == [1, 2]
    assert pk[0]["Score"] == pytest.approx(0.283) and pk[0]["Druggability Score"] == pytest.approx(0.883)
    assert pk[0]["Volume"] == pytest.approx(879.569) and pk[1]["Druggability Score"] == pytest.approx(0.012)
    assert pk[0]["Number of Alpha Spheres"] == 103


def test_pocket_atom_files_give_residue_sets():
    res = S.residues_from_pdb_text((FIX / "fpocket_pocket1_atm.pdb").read_text())
    assert res and all(isinstance(c, str) and isinstance(n, int) for c, n in res)


@pytest.mark.skipif(S.__dict__.get("run_tool") is None and __import__("vsbench.tools", fromlist=["x"]).find_tool("fpocket") is None,
                    reason="fpocket not installed")
def test_run_fpocket_on_a_real_template(tmp_path):
    root = Path(__file__).parents[1] / "data" / "raw" / "full_data" / "ESR1_ant"
    if not (root / "1xqc_protein.mol2").exists():
        pytest.skip("LIT-PCBA data not downloaded")
    pdb = tmp_path / "t.pdb"
    S.write_pdb_from_mol2(root / "1xqc_protein.mol2", pdb)
    pockets = S.run_fpocket(pdb, tmp_path / "fp")
    assert len(pockets) >= 3 and [p.rank for p in pockets] == sorted(p.rank for p in pockets)
    assert all(0.0 <= p.druggability <= 1.0 for p in pockets)
    # The co-crystal site of the ESR1 antagonist template must be among the top fpocket pockets.
    lig = np.array([[float(l.split()[2]), float(l.split()[3]), float(l.split()[4])]
                    for l in (root / "1xqc_ligand.mol2").read_text().split("@<TRIPOS>ATOM")[1].split("@<TRIPOS>BOND")[0].splitlines()
                    if l.strip() and not l.split()[5].startswith("H")])
    native = S.native_site(read_mol2_residues(root / "1xqc_protein.mol2"), lig)
    r = S.site_recovery(pockets, native)
    assert r["top3_hit"]


# ----------------------------------------------------------------- P2Rank
P2 = """  name  ,  rank  ,  score  ,  probability  ,  sas_points  ,  surf_atoms  ,  center_x  ,  center_y  ,  center_z  ,  residue_ids  ,  surf_atom_ids
  pocket1  ,   1  ,   18.77  ,   0.844  ,   128  ,   42  ,   12.3  ,   4.5  ,   -3.1  ,   A_307 A_308 A_310 B_12A  ,   1 2 3
  pocket2  ,   2  ,   6.10  ,   0.210  ,   40  ,   20  ,   30.0  ,   4.5  ,   -3.1  ,   A_400  ,   9 10
"""


def test_p2rank_predictions_parse():
    pk = S.parse_p2rank_predictions(P2)
    assert [p.rank for p in pk] == [1, 2] and pk[0].druggability == pytest.approx(0.844) and pk[0].score == pytest.approx(18.77)
    assert pk[0].residues == frozenset({("A", 307), ("A", 308), ("A", 310), ("B", 12)})   # insertion code ignored


# ----------------------------------------------------------------- DoGSite
def test_dogsite_table_keeps_main_pockets_only_and_reads_scores():
    rows = S.parse_dogsite_table((FIX / "dogsite_desc.txt").read_text())
    assert [r["name"] for r in rows] == ["P_0"]
    assert float(rows[0]["volume"]) == pytest.approx(1739.52) and float(rows[0]["drugScore"]) == pytest.approx(0.809707)


def _fake_http(calls, polls_before_ready=2, throttle_first=False):
    desc = ("name\tvolume\tsimpleScore\tdrugScore\tenclosure\tdepth\tsurface\n"
            "P_0\t900.0\t0.60\t0.30\t0.1\t10\t500\n"
            "P_1\t400.0\t0.70\t0.85\t0.2\t12\t300\n"
            "P_1_0\t100.0\t0.10\t0.10\t0.0\t1\t50\n")
    state = {"polls": 0, "posted": 0}

    def http(method, url, body=None, timeout=60.0):
        calls.append((method, url))
        if method == "POST":
            state["posted"] += 1
            if throttle_first and state["posted"] == 1:
                return 429, json.dumps({"status_code": 429})
            return 202, json.dumps({"location": "https://x/job1"})
        if url == "https://x/job1":
            state["polls"] += 1
            if state["polls"] <= polls_before_ready:
                return 202, json.dumps({"status_code": 202})
            return 200, json.dumps({"status_code": 200, "result_table": "https://x/desc.txt",
                                    "residues": ["https://x/1abc_P_0_res.pdb", "https://x/1abc_P_1_res.pdb", "https://x/1abc_P_1_0_res.pdb"]})
        if url == "https://x/desc.txt":
            return 200, desc
        if url.endswith("1abc_P_0_res.pdb"):
            return 200, "ATOM      1  CA  ALA A 301      0.0   0.0   0.0\n"
        if url.endswith("1abc_P_1_res.pdb"):
            return 200, "ATOM      1  CA  ALA A 302      0.0   0.0   0.0\nATOM      2  CA  ALA A 303      0.0   0.0   0.0\n"
        raise AssertionError(url)
    return http


def test_dogsite_client_polls_ranks_by_druggability_and_caches(tmp_path):
    calls = []
    pk = S.run_dogsite("1ABC", tmp_path, http=_fake_http(calls), sleep=lambda s: None)
    assert [p.extra["name"] for p in pk] == ["P_1", "P_0"] and [p.rank for p in pk] == [1, 2]    # best drugScore first
    assert pk[0].residues == frozenset({("A", 302), ("A", 303)}) and pk[0].druggability == pytest.approx(0.85)
    n = len(calls)
    again = S.run_dogsite("1ABC", tmp_path, http=lambda *a, **k: (_ for _ in ()).throw(AssertionError("server called")))
    assert [p.druggability for p in again] == [p.druggability for p in pk] and len(calls) == n    # served from cache


def test_dogsite_client_backs_off_when_throttled_and_gives_up_eventually(tmp_path):
    calls = []
    pk = S.run_dogsite("1ABC", tmp_path, http=_fake_http(calls, throttle_first=True), sleep=lambda s: None)
    assert len(pk) == 2
    bad = lambda *a, **k: (500, "server error")
    with pytest.raises(RuntimeError, match="submission failed"):
        S.run_dogsite("2XYZ", tmp_path, http=bad, sleep=lambda s: None)


def test_dogsite_timeout(tmp_path):
    def never_ready(method, url, body=None, timeout=60.0):
        return (202, json.dumps({"location": "https://x/j"})) if method == "POST" else (202, json.dumps({"status_code": 202}))
    with pytest.raises(TimeoutError):
        S.run_dogsite("3ABC", tmp_path, http=never_ready, sleep=lambda s: None, poll_s=100.0, max_wait_s=250.0)


# ----------------------------------------------------------------- burial
def test_ligand_buried_fraction_distinguishes_enclosed_from_exposed(tmp_path):
    shell = []
    rng = np.random.default_rng(0)
    pts = rng.normal(size=(400, 3))
    pts = 5.0 * pts / np.linalg.norm(pts, axis=1)[:, None]               # a sphere of protein atoms, radius 5 A
    for i, (x, y, z) in enumerate(pts):
        shell.append(f"ATOM  {i + 1:5d}  C   ALA A{i % 9999 + 1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C")
    p = tmp_path / "shell.pdb"
    p.write_text("\n".join(shell) + "\nEND\n")
    inside = S.ligand_buried_fraction(p, np.array([[0.0, 0.0, 0.0]]))
    outside = S.ligand_buried_fraction(p, np.array([[40.0, 0.0, 0.0]]))
    assert inside > 0.9 and outside < 0.05


def test_unusual_residues_marked_with_asterisks_do_not_break_the_pdb_columns(tmp_path):
    text = _mol2().replace("ALA RESIDUE", "ALA RESIDUE")
    text = text.replace("4 A     ALA     1 ROOT ALA A 303", "4 A     ****     1 ROOT CSO A 303")
    p = tmp_path / "u.mol2"
    p.write_text(text)
    out = tmp_path / "u.pdb"
    S.write_pdb_from_mol2(p, out)
    for l in out.read_text().splitlines():
        if l.startswith("ATOM"):
            assert l[17:20].strip().isalpha() and l[22:26].strip().isdigit(), l           # columns intact
    assert ("A", 303) in S.residues_from_pdb_text(out.read_text())
    assert "CSO" in out.read_text()
    assert [r.name for r in read_mol2_residues(p)][2] == "CSO"
