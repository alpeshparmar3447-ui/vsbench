import math
import shutil
from pathlib import Path

import numpy as np
import pytest

from vsbench import gridstrat as G
from vsbench import sites as S
from vsbench.box import MAX_GRID_POINTS
from vsbench.tools import find_tool


def test_wilson_interval_matches_reference_values():
    lo, hi = G.wilson(0, 10)
    assert lo == pytest.approx(0.0, abs=1e-9) and hi == pytest.approx(0.2775, abs=5e-4)
    lo, hi = G.wilson(5, 10)
    assert (lo, hi) == pytest.approx((0.2366, 0.7634), abs=5e-4)
    lo, hi = G.wilson(10, 10)
    assert hi == pytest.approx(1.0, abs=1e-9) and lo == pytest.approx(0.7225, abs=5e-4)
    assert math.isnan(G.wilson(0, 0)[0])


def test_wilson_interval_never_leaves_the_unit_range():
    """At k = 0 or k = n the closed form lands a hair outside [0, 1]; a reported rate interval of
    -0.00 or 1.0000000002 is wrong on its face, so the bounds are clamped."""
    for n in range(1, 40):
        for k in (0, n):
            lo, hi = G.wilson(k, n)
            assert 0.0 <= lo <= hi <= 1.0, (k, n, lo, hi)


def test_blind_box_covers_the_protein_and_respects_the_grid_limit():
    rng = np.random.default_rng(0)
    xyz = rng.uniform(-20, 20, size=(500, 3))
    b = G.blind_box(xyz, margin=2.0)
    for ax in range(3):
        assert b.center[ax] - b.size[ax] / 2 <= xyz[:, ax].min() + 1e-9
        assert b.center[ax] + b.size[ax] / 2 >= xyz[:, ax].max() - 1e-9
    assert b.spacing == 0.375 and max(b.npts) <= MAX_GRID_POINTS
    big = G.blind_box(np.array([[0, 0, 0], [150.0, 10, 10], [10, 10, 10]]))       # 154 A edge: needs coarser grid
    assert big.spacing > 0.375 and max(big.npts) <= MAX_GRID_POINTS


def test_pocket_box_is_cubic_at_the_given_centre():
    b = G.pocket_box((1.0, 2.0, 3.0), edge=20.0)
    assert b.center == (1.0, 2.0, 3.0) and b.size == (20.0, 20.0, 20.0)


def test_alpha_sphere_center_and_top_pocket_choice():
    pqr = ("ATOM      1  C   STP A   1       0.000   0.000   0.000  0.00  3.50\n"
           "ATOM      2  C   STP A   1       2.000   4.000   6.000  0.00  3.50\n")
    assert S.alpha_sphere_center(pqr) == pytest.approx((1.0, 2.0, 3.0))
    assert S.alpha_sphere_center("") is None
    pockets = [S.Pocket(2, 5.0, 0.5, 100.0, extra={"center": (9.0, 9.0, 9.0)}),
               S.Pocket(1, 9.0, 0.9, 200.0, extra={"center": (1.0, 2.0, 3.0)}),
               S.Pocket(3, 1.0, 0.1, 50.0)]
    assert G.top_pocket_center(pockets) == (1.0, 2.0, 3.0)                  # rank 1 wins
    assert G.top_pocket_center([S.Pocket(1, 1.0, 0.1, 5.0)]) is None


def test_summarize_counts_successes_per_strategy():
    rows = [{"template": f"t{i}", "strategy": s, "ad_success": ok, "gn_success": ok, "ad_lowest_rmsd": r,
             "gn_top_rmsd": r, "wall_s": 10.0}
            for s, vals in (("ligand", [(True, 0.5), (True, 1.0), (False, 3.0)]), ("blind", [(False, 9.0)] * 3))
            for i, (ok, r) in enumerate(vals)]
    t = G.summarize(rows).set_index(["strategy", "engine"])
    assert t.loc[("ligand", "gnina"), "successes"] == 2 and t.loc[("ligand", "gnina"), "n"] == 3
    assert t.loc[("blind", "autodock-gpu"), "success_rate"] == 0.0
    assert t.loc[("ligand", "gnina"), "median_rmsd"] == 1.0


ROOT = Path(__file__).parents[1] / "data" / "raw" / "full_data" / "ESR1_ant"
needs = pytest.mark.skipif(
    not (ROOT / "1xp1_protein.mol2").exists()
    or not all(find_tool(t) for t in ("autodock_gpu", "autogrid4", "gnina", "obabel", "mk_prepare_receptor.py"))
    or shutil.which("nvidia-smi") is None, reason="data, executables or GPU missing")


@pytest.mark.slow
@needs
def test_redocking_with_each_box_strategy_runs_and_a_wrong_site_is_penalised(tmp_path):
    from rdkit import Chem
    from vsbench import redock
    xtal, _ = redock.crystal_ligand(ROOT / "1xp1_ligand.mol2", tmp_path / "x")
    centre = Chem.RemoveHs(xtal).GetConformer().GetPositions().mean(axis=0)
    S.write_pdb_from_mol2(ROOT / "1xp1_protein.mol2", tmp_path / "p.pdb")
    prot_xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                         for l in (tmp_path / "p.pdb").read_text().splitlines() if l.startswith("ATOM")])
    kw = dict(seed=1, nrun=10, exhaustiveness=4, num_modes=2, cpu=8)
    rows = {}
    for name, box in (("pocket", G.pocket_box(centre)), ("blind", G.blind_box(prot_xyz)),
                      ("wrong_site", G.pocket_box(centre + np.array([45.0, 0, 0])))):
        rows[name] = G.run_strategy(name, ROOT / "1xp1_protein.mol2", ROOT / "1xp1_ligand.mol2", "1xp1",
                                    tmp_path / name, box, **kw)
        assert rows[name]["ok"], rows[name]["error"]
    assert rows["pocket"]["gn_top_rmsd"] < 4.0                              # box on the site finds it
    assert rows["wrong_site"]["gn_top_rmsd"] > 10.0 and not rows["wrong_site"]["gn_success"]
    assert np.isfinite(rows["blind"]["gn_top_rmsd"]) and rows["blind"]["box_edge"] > 40.0
