import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from vsbench import perf as P
from vsbench.tools import find_tool


def test_parse_smi():
    assert P.parse_smi("1234 MiB, 56 %\n") == (1234.0, 56.0)
    assert P.parse_smi("7 MiB, 0 %") == (7.0, 0.0)


def test_summarize_speedup_is_relative_to_the_reference_at_its_best_thread_count():
    runs = pd.DataFrame([
        {"case": "vina (CPU)", "n": 10, "threads": 1, "wall_s": 100.0, "n_ok": 10, "peak_vram_mib": np.nan, "mean_gpu_util": np.nan, "repeat": 0},
        {"case": "vina (CPU)", "n": 10, "threads": 1, "wall_s": 120.0, "n_ok": 10, "peak_vram_mib": np.nan, "mean_gpu_util": np.nan, "repeat": 1},
        {"case": "vina (CPU)", "n": 10, "threads": 10, "wall_s": 20.0, "n_ok": 10, "peak_vram_mib": np.nan, "mean_gpu_util": np.nan, "repeat": 0},
        {"case": "autodock-gpu", "n": 10, "threads": 1, "wall_s": 5.0, "n_ok": 10, "peak_vram_mib": 300.0, "mean_gpu_util": 80.0, "repeat": 0},
    ])
    t = P.summarize(runs).set_index(["case", "threads"])
    assert t.loc[("vina (CPU)", 1), "wall_mean_s"] == 110.0
    assert t.loc[("vina (CPU)", 1), "wall_sd_s"] == pytest.approx(np.std([100, 120], ddof=1))
    assert t.loc[("autodock-gpu", 1), "s_per_compound"] == 0.5
    assert t.loc[("autodock-gpu", 1), "speedup_vs_ref"] == pytest.approx(2.0 / 0.5)          # 2 s/compound best CPU
    assert t.loc[("vina (CPU)", 10), "speedup_vs_ref"] == pytest.approx(1.0)
    assert t.loc[("autodock-gpu", 1), "compounds_per_hour"] == pytest.approx(7200.0)
    assert P.projected_hours(0.5, 7200) == pytest.approx(1.0)


def test_unknown_case_is_rejected(tmp_path):
    with pytest.raises(ValueError):
        P.run_case("not-a-case", None, [], tmp_path / "w", 1)


def test_gpu_sampler_without_nvidia_smi_reports_nan(monkeypatch):
    monkeypatch.setattr(P.shutil, "which", lambda _: None)
    with P.GpuSampler() as g:
        pass
    assert np.isnan(g.summary()["peak_vram_mib"])


ROOT = Path(__file__).parents[1] / "data" / "raw" / "full_data"
needs = pytest.mark.skipif(
    not (ROOT / "ESR1_ant").is_dir()
    or not all(find_tool(t) for t in ("autodock_gpu", "autogrid4", "gnina", "vina", "obabel", "mk_prepare_receptor.py"))
    or shutil.which("nvidia-smi") is None, reason="data, executables or GPU missing")


@pytest.mark.slow
@needs
def test_every_case_runs_on_real_inputs_and_reports_sane_numbers(tmp_path):
    from vsbench import data, engines, ligprep, redock
    from vsbench.box import box_from_coords
    from rdkit import Chem
    target = data.load_target(ROOT, "ESR1_ant")
    xtal, _ = redock.crystal_ligand(target.ligand_mol2("1xp1"), tmp_path / "x")
    box = box_from_coords(Chem.RemoveHs(xtal).GetConformer().GetPositions())
    rec = engines.prepare_receptor(target.protein_mol2("1xp1"), tmp_path / "rec", box)
    pool = pd.concat([target.actives.head(3), target.inactives.head(3)])
    ligs = ligprep.prepare_many(list(zip(pool.smiles, pool.id)), n_jobs=4)
    plan = [(c, 4, 4) for c in P.CASES]
    runs = P.run_suite(rec, ligs, tmp_path / "perf", plan, repeats=2, log=lambda *_: None, exhaustiveness=4, nrun=10)
    assert len(runs) == len(P.CASES) * 2 and (runs.n_ok == 4).all(), runs
    t = P.summarize(runs).set_index("case")
    assert (t.wall_mean_s > 0).all()
    assert t.loc["autodock-gpu", "peak_vram_mib"] >= 0
    # the GPU engine is faster than the CPU Vina reference on the same ligands
    assert t.loc["autodock-gpu", "speedup_vs_ref"] > 1.0
