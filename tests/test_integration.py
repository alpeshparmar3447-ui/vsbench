"""End-to-end tests with the real executables, receptor files and GPU."""
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from vsbench import data, engines, ligprep, redock
from vsbench.tools import find_tool

ROOT = Path(__file__).parents[1] / "data" / "raw" / "full_data"
TARGET, TEMPLATE = "ESR1_ant", "1xp1"

needs_data = pytest.mark.skipif(not (ROOT / TARGET).is_dir(), reason="LIT-PCBA data not downloaded")
needs_tools = pytest.mark.skipif(
    not all(find_tool(t) for t in ("autodock_gpu", "autogrid4", "gnina", "obabel", "mk_prepare_receptor.py"))
    or shutil.which("nvidia-smi") is None,
    reason="docking executables or GPU not available",
)
pytestmark = [pytest.mark.slow, needs_data, needs_tools]


@pytest.fixture(scope="module")
def target():
    return data.load_target(ROOT, TARGET)


@pytest.fixture(scope="module")
def workdir(tmp_path_factory):
    return tmp_path_factory.mktemp("vsbench_it")


@pytest.fixture(scope="module")
def prepared(target, workdir):
    xtal, smiles = redock.crystal_ligand(target.ligand_mol2(TEMPLATE), workdir / "x")
    from rdkit import Chem
    heavy = Chem.RemoveHs(xtal)
    from vsbench.box import box_from_coords
    box = box_from_coords(heavy.GetConformer().GetPositions())
    rec = engines.prepare_receptor(target.protein_mol2(TEMPLATE), workdir / "rec", box)
    return heavy, smiles, box, rec


def test_crystal_ligand_is_read_with_the_expected_atoms(prepared):
    heavy, smiles, box, _ = prepared
    assert heavy.GetNumAtoms() == 34
    assert Chem_ok(smiles)
    assert box.size[0] >= 24.0                       # generous default box applies


def Chem_ok(smiles):
    from rdkit import Chem
    return Chem.MolFromSmiles(smiles) is not None


def test_receptor_preparation_produces_maps_and_is_cached(prepared, target, workdir):
    _, _, box, rec = prepared
    assert rec.pdbqt.exists() and rec.fld.exists()
    assert len(list(rec.workdir.glob("receptor.*.map"))) >= 10
    stamp = (rec.workdir / "stamp.json").stat().st_mtime_ns
    again = engines.prepare_receptor(target.protein_mol2(TEMPLATE), workdir / "rec", box)
    assert again.fld == rec.fld
    assert (rec.workdir / "stamp.json").stat().st_mtime_ns == stamp   # nothing was rebuilt


def test_both_engines_dock_a_small_batch_and_resume_without_redoing_work(prepared, target, workdir):
    _, _, _, rec = prepared
    pool = pd.concat([target.actives.head(2), target.inactives.head(2)])
    ligs = ligprep.prepare_many(list(zip(pool.smiles, pool.id)), n_jobs=4)
    assert all(l.ok for l in ligs)

    ad, info = engines.run_autodock_gpu(rec, ligs, workdir / "ad", nrun=10, seed=7)
    assert ad["ad_ok"].all() and np.isfinite(ad["ad_energy"]).all()
    assert (ad["ad_energy"] < 5).all() and (ad["ad_energy"] > -30).all()     # physical range
    assert (ad["ad_largest_size"] >= 1).all() and info["n_docked_now"] == 4

    ad2, info2 = engines.run_autodock_gpu(rec, ligs, workdir / "ad", nrun=10, seed=7)
    assert info2["n_docked_now"] == 0                                       # idempotent resume
    pd.testing.assert_frame_equal(ad, ad2)

    gn, ginfo = engines.run_gnina(rec, ligs, workdir / "gn", num_modes=2, seed=7, cpu=8)
    assert gn["gn_ok"].all() and ginfo["n_docked_now"] == 4
    assert gn["gn_cnn_vs"].between(0, 20).all()
    np.testing.assert_allclose(gn["gn_cnn_vs"], gn["gn_cnn_score"] * gn["gn_cnn_affinity"], rtol=1e-4)
    gn2, ginfo2 = engines.run_gnina(rec, ligs, workdir / "gn", num_modes=2, seed=7, cpu=8)
    assert ginfo2["n_docked_now"] == 0
    pd.testing.assert_frame_equal(gn, gn2)


def test_autodock_results_are_reproducible_for_a_fixed_seed(prepared, target, workdir):
    _, _, _, rec = prepared
    pool = target.actives.head(2)
    ligs = ligprep.prepare_many(list(zip(pool.smiles, pool.id)), n_jobs=2)
    a, _ = engines.run_autodock_gpu(rec, ligs, workdir / "rep_a", nrun=10, seed=11)
    b, _ = engines.run_autodock_gpu(rec, ligs, workdir / "rep_b", nrun=10, seed=11)
    assert np.abs(a["ad_energy"] - b["ad_energy"]).max() <= 0.05   # documented GPU non-associativity


def test_template_redocking_runs_end_to_end(target, workdir):
    res = redock.redock_template(target.protein_mol2(TEMPLATE), target.ligand_mol2(TEMPLATE),
                                 TEMPLATE, workdir / "redock", seed=42, nrun=20, cpu=8)
    assert res.ok, res.error
    assert res.n_heavy == 34
    for v in (res.ad_lowest_rmsd, res.ad_largest_rmsd, res.gn_top_rmsd, res.gn_best_rmsd):
        assert np.isfinite(v) and 0.0 <= v < 60.0
    assert res.gn_best_rmsd <= res.gn_top_rmsd + 1e-9              # best-of-N cannot exceed top-1
    assert res.ad_energy < 0


def test_a_broken_template_is_reported_not_raised(workdir, tmp_path):
    bad = tmp_path / "bad.mol2"
    bad.write_text("not a mol2")
    res = redock.redock_template(bad, bad, "bad", workdir / "bad", seed=1)
    assert res.ok is False and res.error


def test_gate_ablation_runs_all_four_conditions_and_the_gates_fire_where_they_should(prepared, target, workdir):
    from vsbench import ablation
    _, _, _, rec = prepared
    pool = pd.concat([target.actives.head(3), target.inactives.head(3)])
    ligs = ligprep.prepare_many(list(zip(pool.smiles, pool.id)), n_jobs=4)
    wide, summ = ablation.run_ablation(rec, ligs, workdir / "ablation", nrun=10, n_boot=50)
    assert list(summ.condition) == list(ablation.CONDITIONS) and len(wide) == 6
    s = summ.set_index("condition")
    assert s.loc["baseline", "n_charge_flagged"] == 0 and s.loc["charges_x4_uncapped", "n_charge_flagged"] > 0
    assert s.loc["charges_x4_capped", "n_charge_flagged"] == 0                      # the cap removes the flag
    assert s.loc["torsions_frozen", "n_torsion_frozen_flagged"] >= 1                 # the gate sees frozen flexible ligands
    assert s.loc["baseline", "mean_shift_vs_baseline"] == 0.0
    assert np.isfinite(wide.to_numpy()).all()
