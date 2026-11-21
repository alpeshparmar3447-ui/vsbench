"""The whole pipeline at tiny scale: real data, real engines, real GPU."""
import json
import re
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from vsbench import cli, pipeline
from vsbench.tools import find_tool

ROOT = Path(__file__).parents[1] / "data" / "raw" / "full_data"
pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(not (ROOT / "ESR1_ant").is_dir(), reason="LIT-PCBA data not downloaded"),
    pytest.mark.skipif(
        not all(find_tool(t) for t in ("autodock_gpu", "autogrid4", "gnina", "obabel", "mk_prepare_receptor.py"))
        or shutil.which("nvidia-smi") is None, reason="docking executables or GPU not available"),
]


@pytest.fixture(scope="module")
def config_path(tmp_path_factory):
    d = tmp_path_factory.mktemp("e2e")
    cfg = {
        "data_root": str(ROOT), "results_dir": str(d / "results"), "targets": ["ESR1_ant"],
        "n_actives": 4, "n_random": 8, "per_active_matched": 2, "seed": 5,
        "templates": {"ESR1_ant": ["1xp1", "1xqc"]},
        "autodock": {"nrun": 10, "autostop": True, "heuristics": True},
        "gnina": {"exhaustiveness": 4, "num_modes": 2, "cpu": 8},
        "n_boot": 50, "n_jobs": 4, "decoy_per_active": 2, "decoy_pool_per_target": 300,
        "seeds": [5, 6], "gnina_seed_study": {"n": 12, "seeds": [7]},
        "lead_selection": {"top_k": 10, "top_n": [3, 5], "admet": True, "interactions": True, "n_jobs": 4},
    }
    p = d / "cfg.json"
    p.write_text(json.dumps(cfg))
    return p


@pytest.fixture(scope="module")
def run_once(config_path):
    assert cli.main(["run", "--config", str(config_path)]) == 0
    return pipeline.Config.from_json(config_path)


def test_template_stage_validates_and_selects_without_labels(run_once):
    tdir = run_once.target_dir("ESR1_ant")
    redock = pd.read_csv(tdir / "redock.csv")
    assert set(redock.template) == {"1xp1", "1xqc"}
    assert redock.ok.all()
    sel = json.loads((tdir / "selected.json").read_text())
    worst = redock.assign(w=redock[["ad_lowest_rmsd", "gn_top_rmsd"]].max(axis=1)).sort_values(["w", "template"])
    assert sel["template"] == worst.iloc[0].template            # the documented rule, nothing else
    assert sel["n_templates"] == 2


def test_screen_table_has_the_designed_composition_and_finite_scores(run_once):
    df = pd.read_csv(run_once.target_dir("ESR1_ant") / "screen_scores.csv")
    assert (df.neg_set == "active").sum() == 4
    assert (df.neg_set == "random").sum() == 8
    assert (df.neg_set == "matched").sum() == 8                 # 4 actives x 2 per active
    assert (df.neg_set == "decoy").sum() == 8                   # DUD-E-protocol decoys, 2 per active
    assert {"ad_energy_s5", "ad_energy_s6", "ad_energy_sd", "gn_cnn_vs_s7", "nn_sim", "frozen_flexible"} <= set(df.columns)
    assert df["gn_cnn_vs_s7"].notna().sum() == 12               # the gnina seed study re-docked the subset
    assert df.name.is_unique
    assert df.prep_ok.all() and df.ad_ok.all() and df.gn_ok.all()
    assert df.ad_energy.between(-30, 10).all()
    assert df.gn_cnn_vs.between(0, 20).all()
    info = json.loads((run_once.target_dir("ESR1_ant") / "screen_info.json").read_text())
    assert info["autodock"]["wall_s"] > 0 and info["gnina"]["wall_s"] > 0


def test_analysis_outputs_exist_and_are_internally_consistent(run_once):
    adir = run_once.target_dir("ESR1_ant") / "analysis"
    for f in ("metrics", "baselines", "pose_selection", "engine_difference", "balance", "failures"):
        assert (adir / f"{f}.csv").exists(), f
    m = pd.read_csv(adir / "metrics.csv")
    assert set(m.negative_set) == {"random", "matched", "decoy"}
    for f in ("seed_stability", "gnina_seed_study", "negative_set_differences"):
        assert (adir / f"{f}.csv").exists(), f
    assert (adir / "similarity.json").exists()
    assert set(m.metric) == {"roc_auc", "pr_auc", "bedroc20", "ef1", "ef5", "nef1"}
    auc = m[m.metric == "roc_auc"]
    assert auc.value.between(0, 1).all() and (auc.ci_low <= auc.value + 1e-9).all()
    macro = pd.read_csv(Path(run_once.results_dir) / "metrics_macro.csv")
    assert (macro.n_targets == 1).all()


def test_lead_selection_outputs_are_complete_and_consistent(run_once):
    tdir = run_once.target_dir("ESR1_ant")
    adir = tdir / "analysis"
    adm = pd.read_csv(tdir / "admet.csv")
    assert len(adm) == 28 and adm["name"].is_unique
    assert set(adm.tier) <= {"PASS", "MODERATE", "HIGH_RISK", "EXCLUDE"}
    assert adm[["hERG", "AMES"]].apply(lambda c: c.between(0, 1)).all().all()
    ifp = pd.read_csv(tdir / "interactions.csv")
    assert set(ifp.name) >= set(pd.read_csv(tdir / "screen_scores.csv").query("neg_set == 'active'").name)
    ok = ifp[ifp.ok.astype(bool)]
    assert len(ok) > 0 and ok.recall.between(0, 1).all() and ok.precision.between(0, 1).all()
    assert (ifp.n_reference == ifp.n_reference.iloc[0]).all() and ifp.n_reference.iloc[0] >= 4
    pri = pd.read_csv(adir / "prioritization.csv")
    assert set(pri.negative_set) == {"random", "matched", "decoy"}
    assert {"docking", "docking+ADMET", "docking+IFP", "pareto"} <= set(pri.strategy)
    assert (pri.actives_in_top_n <= pri.top_n).all() and pri.hit_rate.between(0, 1).all()
    assert (adir / "admet_filters.csv").exists()


def test_rerunning_the_pipeline_does_no_new_docking(run_once, config_path):
    tdir = run_once.target_dir("ESR1_ant")
    before = json.loads((tdir / "screen_info.json").read_text())
    scores_before = pd.read_csv(tdir / "screen_scores.csv")
    assert cli.main(["run", "--config", str(config_path)]) == 0
    after = json.loads((tdir / "screen_info.json").read_text())
    assert after["autodock"]["n_docked_now"] == 0 and after["gnina"]["n_docked_now"] == 0
    n_compounds = len(scores_before)                            # 4 actives + 8 random + 8 matched + 8 decoys
    assert n_compounds == 28
    for eng in ("autodock", "gnina"):                           # first-run timings must survive a resume
        assert after[eng]["n_docked_cumulative"] == n_compounds == before[eng]["n_docked_cumulative"]
        assert after[eng]["wall_s_cumulative"] >= before[eng]["wall_s_cumulative"] > 1.0
    scores_after = pd.read_csv(tdir / "screen_scores.csv")
    pd.testing.assert_frame_equal(scores_before[["name", "ad_energy", "gn_cnn_vs"]],
                                  scores_after[["name", "ad_energy", "gn_cnn_vs"]])
    assert before["template"] == after["template"]


def test_cli_rejects_unknown_targets(config_path):
    with pytest.raises(SystemExit):
        cli.main(["run", "--config", str(config_path), "--targets", "NOPE"])


def _assert_gapless_numbering(text):
    """Main tables run 1..N and supplementary S1..SM with no gap and no repeat."""
    main = [int(n) for n in re.findall(r"^\*\*Table (\d+)\. ", text, re.M)]
    supp = [int(n) for n in re.findall(r"^\*\*Table S(\d+)\. ", text, re.M)]
    assert main == list(range(1, len(main) + 1)), main
    assert supp == list(range(1, len(supp) + 1)), supp
    return main, supp


def test_report_renders_every_table_from_real_output_files(run_once):
    from vsbench import report
    results = Path(run_once.results_dir)
    text = report.build_report(results, ["ESR1_ant"])
    # Asserted by analysis, not by number: the numbers follow position, so hard-coding them
    # here would break on any reordering without saying anything about the content.
    for analysis in ("Panel, template validation and compounds docked",
                     "ROC-AUC against each negative set",
                     "ROC-AUC of single properties used as scores",
                     "Pose selection", "Throughput", "Failures"):
        assert analysis in text, analysis
    main, supp = _assert_gapless_numbering(text)
    assert main and supp
    assert "pocket ions" in text and "deposition checked" in text
    assert (results / "tables.md").read_text() == text
    assert (results / "tables_main.md").read_text() in text
    assert (results / "tables_supplementary.md").read_text() in text
    assert "nan" not in text.lower().replace("nan%", "")
    panel = report.panel_table(Path(run_once.results_dir), ["ESR1_ant"])
    assert panel.loc[0, "actives"] == 4 and panel.loc[0, "matched inactives"] == 8
    runtime = report.runtime_table(Path(run_once.results_dir), ["ESR1_ant"])
    assert float(runtime.loc[0, "AutoDock-GPU (s/compound)"]) > 0


def test_figures_are_rendered_from_the_real_outputs(run_once, tmp_path):
    from vsbench import studies
    made = studies.make_figures(run_once)
    names = {p.name for p in made}
    assert {"auc_summary.png", "negative_sets_consensus.png", "redocking.png", "prioritization.png"} <= names
    assert all(p.exists() and p.stat().st_size > 5000 for p in made)
    from vsbench import report
    text = report.build_report(Path(run_once.results_dir), ["ESR1_ant"])
    for analysis in ("DUD-E-protocol decoys",
                     "Change in ROC-AUC between negative-set designs",
                     "Seed stability of AutoDock-GPU",
                     "Similarity of the negative sets to the actives",
                     "ADMET risk tiers", "Effect of ADMET filters on retrieval",
                     "Hit rate of the lead-selection strategies",
                     "Recovery of the co-crystal interactions by the top gnina pose"):
        assert analysis in text, analysis
    _assert_gapless_numbering(text)


def test_ablation_study_runs_on_the_screened_receptor(run_once):
    from vsbench import ablation, studies
    t = studies.ablation_study(run_once, ["ESR1_ant"], n_each=3, log=lambda *_: None)
    assert list(t.condition) == list(ablation.CONDITIONS) and (t.n_ligands == 6).all()
    s = t.set_index("condition")
    assert s.loc["baseline", "n_charge_flagged"] == 0 and s.loc["charges_x4_uncapped", "n_charge_flagged"] > 0
    assert (Path(run_once.results_dir) / "ablation_summary.csv").exists()


def test_grid_and_performance_studies_run_end_to_end(run_once):
    from dataclasses import replace
    from vsbench import studies
    small = replace(run_once, templates={"ESR1_ant": ["1xqc"]})
    from vsbench import gridstrat, receptors
    # p2rank only runs where P2Rank and a JRE are installed, so the expected set depends on
    # the machine; the ones needing nothing beyond the docking tools must always be there.
    expected = {"ligand", "pocket", "blind"}
    if receptors.options(small)["p2rank"]:
        expected.add("p2rank")
    g = studies.grid_study(small, ["ESR1_ant"], log=lambda *_: None)
    assert set(g.strategy) == expected and g.ok.all(), g[["strategy", "error"]]
    assert set(g.strategy) <= set(gridstrat.STRATEGIES)
    assert g.loc[g.strategy == "blind", "box_edge"].iloc[0] > g.loc[g.strategy == "ligand", "box_edge"].iloc[0]
    summ = pd.read_csv(Path(run_once.results_dir) / "grid_strategies_summary.csv")
    assert set(summ.strategy) == expected
    # the paired comparison is written beside the unpaired rates, not left to be done by hand
    paired = pd.read_csv(Path(run_once.results_dir) / "grid_strategies_paired.csv")
    assert set(paired.scope) == {"pairwise", "listwise"}
    assert paired.p_mcnemar_exact.between(0, 1).all()
    p = studies.perf_study(run_once, "ESR1_ant", n_values=(4,), cpu_threads=(2,), cpu_n=4, repeats=1, log=lambda *_: None)
    assert {"autodock-gpu", "vina (CPU)", "gnina (CNN, CPU)"} <= set(p.case) and (p.n_ok == 4).all()
