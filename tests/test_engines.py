from pathlib import Path

import numpy as np
import pytest

from vsbench import engines as E
from vsbench.poses import pose_to_mol

FIX = Path(__file__).parent / "fixtures"


@pytest.fixture(scope="module")
def dlg_text():
    return (FIX / "sample.dlg").read_text()


# ---------------------------------------------------------- DLG parsing
def test_parse_clusters_reads_the_real_histogram(dlg_text):
    clusters = E.parse_clusters(dlg_text)
    assert clusters == [
        {"rank": 1, "energy": -9.76, "run": 10, "mean_energy": -9.21, "size": 18},
        {"rank": 2, "energy": -8.77, "run": 3, "mean_energy": -8.74, "size": 2},
    ]
    assert sum(c["size"] for c in clusters) == 20  # every run is accounted for


def test_summarize_dlg(dlg_text):
    s = E.summarize_dlg(dlg_text)
    assert s == {"energy": -9.76, "run": 10, "n_clusters": 2, "largest_size": 18,
                 "largest_energy": -9.76, "largest_run": 10}


def test_summarize_picks_largest_cluster_and_breaks_ties_by_energy():
    def table(rows):
        return "\n".join(f"{r} | {e:9.2f} | {run:3d} | {e + 0.1:9.2f} | {n:3d} |####"
                         for r, e, run, n in rows)
    text = table([(1, -10.0, 5, 3), (2, -9.0, 8, 9), (3, -8.5, 2, 9), (4, -8.0, 1, 1)])
    s = E.summarize_dlg(text)
    assert s["energy"] == -10.0 and s["run"] == 5
    assert s["largest_size"] == 9 and s["largest_energy"] == -9.0 and s["largest_run"] == 8
    with pytest.raises(ValueError, match="no clustering histogram"):
        E.summarize_dlg("nothing here")


def test_run_blocks_are_keyed_by_run_and_rebuild_into_molecules(dlg_text):
    blocks = E.dlg_run_blocks(dlg_text)
    assert sorted(blocks) == [1, 2]
    mol = pose_to_mol(dlg_text, blocks[1])      # mapping comes from the INPUT-LIGAND-PDBQT remarks
    from vsbench.poses import parse_mapping
    from rdkit import Chem
    smiles, _ = parse_mapping(dlg_text)
    assert mol.GetNumAtoms() == Chem.MolFromSmiles(smiles).GetNumAtoms() == 25
    assert mol.GetConformer().GetPositions().std(axis=0).min() > 0.5


def test_dlg_complete(tmp_path, dlg_text):
    full = tmp_path / "a.dlg"
    full.write_text(dlg_text)
    cut = tmp_path / "b.dlg"
    cut.write_text(dlg_text.split("CLUSTERING HISTOGRAM")[0])
    assert E.dlg_complete(full)
    assert not E.dlg_complete(cut)
    assert not E.dlg_complete(tmp_path / "missing.dlg")


# --------------------------------------------------------------- gnina
def test_gnina_sdf_poses_and_top_row():
    poses = E.read_gnina_sdf(FIX / "gnina_sample.sdf")
    assert list(poses) == ["144203670", "144203993"]
    assert all(len(v) == 3 for v in poses.values())
    row = E.gnina_top_pose_row("144203670", poses["144203670"])
    assert row["gn_ok"] is True
    assert row["gn_vina"] == pytest.approx(-9.06585)
    assert row["gn_cnn_score"] == pytest.approx(0.7882131934)
    assert row["gn_cnn_affinity"] == pytest.approx(6.8496384621)
    assert row["gn_cnn_vs"] == pytest.approx(row["gn_cnn_score"] * row["gn_cnn_affinity"], rel=1e-6)
    clash = E.gnina_top_pose_row("144203993", poses["144203993"])
    assert clash["gn_vina"] > 0          # a real clashing pose is kept, not hidden


def test_gnina_missing_ligand_becomes_a_failure_row():
    row = E.gnina_top_pose_row("x", None)
    assert row["gn_ok"] is False and np.isnan(row["gn_cnn_vs"])


# ------------------------------------------------------------ utilities
def test_unsafe_names_are_rejected():
    E.check_names(["144203670", "a_b-1.2"])
    for bad in ("../x", "a b", "x;rm", ""):
        with pytest.raises(ValueError, match="unsafe"):
            E.check_names([bad])


def test_histidine_tie_message_is_recognised():
    msg = "... for residue_key='A:101', 2 have passed the filter and are tied for fewest atoms ..."
    m = E._HIE_TIE_RE.search(msg.replace("filter and are tied for fewest", "filter, tied for fewest"))
    assert m and m.group(1) == "A:101"


def test_clusters_with_positive_energies_are_parsed():
    """AutoDock-GPU prints positive energies with an explicit '+'. A parser that only accepts '-' drops those
    rows: an all-positive docking looked like a failed run, and a positive largest cluster vanished."""
    table = (
        "   1 |     -3.10 |   4 |     -3.05 |   2 |##\n"
        "   2 |     +0.50 |   7 |     +1.20 |   6 |######\n"          # the most populated cluster is positive
        "   3 |    +12.40 |   1 |    +12.40 |   1 |#\n")
    rows = E.parse_clusters(table)
    assert [r["energy"] for r in rows] == [-3.10, 0.50, 12.40]
    s = E.summarize_dlg(table)
    assert s["energy"] == -3.10 and s["n_clusters"] == 3
    assert s["largest_size"] == 6 and s["largest_energy"] == 0.50 and s["largest_run"] == 7
    only_positive = E.summarize_dlg("   1 |     +2.28 |  10 |     +2.28 |   1 |#\n")
    assert only_positive["energy"] == 2.28 and only_positive["run"] == 10


def _fake_prepare_env(monkeypatch, tmp_path, autogrid_results):
    """prepare_receptor with the preparation tools replaced; AutoGrid returns the given codes in order."""
    import subprocess

    from vsbench.box import Box
    calls = []

    def fake_run_tool(name, args, cwd=None, timeout=None, check=False):
        if name == "mk_prepare_receptor.py":
            (cwd / "receptor.pdbqt").write_text("ATOM      1  C   ALA A   1       1.000   2.000   3.000  1.00  0.00    +0.000 C\n")
            (cwd / "receptor.gpf").write_text("gpf")
            return subprocess.CompletedProcess([], 0, "", "")
        calls.append(name)
        code = autogrid_results[len(calls) - 1]
        if code == 0:
            (cwd / "receptor.maps.fld").write_text("fld")
        return subprocess.CompletedProcess([], code, "", "ERROR bad x,y,z" if code else "")

    src = tmp_path / "p.pdb"
    src.write_text("ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C\nEND\n")
    monkeypatch.setattr(E, "run_tool", fake_run_tool)
    return src, Box(center=(0.0, 0.0, 0.0), size=(24.0, 24.0, 24.0)), calls


def test_a_single_autogrid_failure_is_retried_and_recorded(monkeypatch, tmp_path):
    src, box, calls = _fake_prepare_env(monkeypatch, tmp_path, [1, 0])
    rec = E.prepare_receptor(src, tmp_path / "w", box)
    assert calls == ["autogrid4", "autogrid4"] and rec.report.autogrid_attempts == 2
    again = E.prepare_receptor(src, tmp_path / "w", box)                    # cached: no further AutoGrid call
    assert len(calls) == 2 and again.report.autogrid_attempts == 2


def test_a_repeated_autogrid_failure_raises_with_both_messages(monkeypatch, tmp_path):
    src, box, calls = _fake_prepare_env(monkeypatch, tmp_path, [1, 1])
    with pytest.raises(RuntimeError, match="failed twice.*bad x,y,z.*bad x,y,z"):
        E.prepare_receptor(src, tmp_path / "w", box)
    assert len(calls) == 2


def test_a_clean_run_records_one_attempt(monkeypatch, tmp_path):
    src, box, calls = _fake_prepare_env(monkeypatch, tmp_path, [0])
    assert E.prepare_receptor(src, tmp_path / "w", box).report.autogrid_attempts == 1


def test_reported_wall_time_ignores_a_wall_clock_jump(monkeypatch, tmp_path):
    """A laptop that sleeps mid-run moves the wall clock forward (observed: 9.7 h) while no compute happens.
    Engine timings feed the throughput table, so they must be measured on a monotonic clock."""
    import subprocess
    import time
    from types import SimpleNamespace

    from rdkit import Chem

    real, offset = time.time, [0.0]
    monkeypatch.setattr(time, "time", lambda: real() + offset[0])

    def fake_run_tool(name, args, cwd=None, timeout=None, check=False):
        offset[0] += 36000.0                                   # ten hours pass while the engine "runs"
        return subprocess.CompletedProcess([], 0, "", "")

    monkeypatch.setattr(E, "run_tool", fake_run_tool)
    block = Chem.MolToMolBlock(Chem.AddHs(Chem.MolFromSmiles("CCO")))
    lig = SimpleNamespace(name="a", ok=True, pdbqt="REMARK x\n", molblock=block)
    rec = SimpleNamespace(fld=tmp_path / "r.maps.fld", pdbqt=tmp_path / "r.pdbqt",
                          box=SimpleNamespace(center=(0.0, 0.0, 0.0), size=(20.0, 20.0, 20.0)))
    _, ad_info = E.run_autodock_gpu(rec, [lig], tmp_path / "ad")
    _, gn_info = E.run_gnina(rec, [lig], tmp_path / "gn")
    assert time.time() - real() > 70000                        # the simulated jump really happened (two calls)
    assert ad_info["wall_s"] < 60 and gn_info["wall_s"] < 60
