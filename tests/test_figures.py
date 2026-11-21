import numpy as np
import pandas as pd

from vsbench import figures as F


def _metrics():
    rng = np.random.default_rng(0)
    rows = []
    for t in ("T1", "T2", "T3"):
        for neg in ("random", "matched"):
            for eng in ("autodock-gpu", "gnina", "consensus"):
                v = rng.uniform(0.45, 0.8)
                rows.append(dict(target=t, negative_set=neg, score=eng, metric="roc_auc", value=v,
                                 ci_low=v - 0.07, ci_high=v + 0.07))
    return pd.DataFrame(rows)


def _baselines():
    return pd.DataFrame([dict(target=t, negative_set=n, property="n_heavy", roc_auc=0.6)
                         for t in ("T1", "T2", "T3") for n in ("random", "matched")])


def _png(path):
    assert path.exists() and path.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n" and path.stat().st_size > 5000


def test_auc_summary_with_and_without_baselines(tmp_path):
    _png(F.plot_auc_summary(_metrics(), _baselines(), tmp_path / "a.png"))
    _png(F.plot_auc_summary(_metrics(), None, tmp_path / "b.png"))


def test_roc_figure_from_a_scored_table(tmp_path):
    rng = np.random.default_rng(1)
    n_a, n_i = 30, 200
    df = pd.DataFrame({"label": [1] * n_a + [0] * n_i,
                       "neg_set": ["active"] * n_a + ["random"] * n_i})
    for c in ("ad_lowest", "gn_cnn_vs", "consensus"):
        df[c] = rng.normal(size=len(df)) + 0.6 * df["label"]
    _png(F.plot_roc(df, "random", tmp_path / "roc.png", "demo"))


def test_redock_figure_ignores_failed_templates(tmp_path):
    rows = [dict(target=t, template=f"x{i}", ok=(i != 2), ad_lowest_rmsd=0.5 + i, gn_top_rmsd=1.0 + 2 * i)
            for t in ("T1", "T2") for i in range(5)]
    _png(F.plot_redock(pd.DataFrame(rows), tmp_path / "r.png"))


def test_extended_figures_render_from_tables(tmp_path):
    from vsbench import figures as F
    rows = [dict(target=t, negative_set=n, score="consensus", metric="roc_auc", value=0.6 - 0.05 * i, ci_low=0.5, ci_high=0.7)
            for t in ("A", "B") for i, n in enumerate(("random", "matched", "decoy"))]
    sites = pd.DataFrame({"method": ["fpocket"] * 4 + ["dogsite"] * 3, "top1_hit": [1, 0, 1, 0, 1, 1, 0],
                          "top3_hit": [1, 1, 1, 0, 1, 1, 1]})
    prio = pd.DataFrame([dict(negative_set="matched", strategy=s, top_n=n, hit_rate=0.3, pool_active_fraction=0.25)
                         for s in ("docking", "docking+ADMET", "pareto") for n in (10, 20)] * 2)
    grid = pd.DataFrame([dict(strategy=s, engine=e, success_rate=0.5, ci_low=0.3, ci_high=0.7)
                         for s in ("ligand", "pocket", "blind") for e in ("autodock-gpu", "gnina")])
    perf = pd.DataFrame([dict(case="autodock-gpu", threads=1, n=10, compounds_per_hour=5000.0),
                         dict(case="vina (CPU)", threads=16, n=10, compounds_per_hour=300.0),
                         dict(case="vina (CPU)", threads=1, n=10, compounds_per_hour=40.0)])
    outs = [F.plot_negative_set_auc(pd.DataFrame(rows), "consensus", tmp_path / "a.png"),
            F.plot_site_recovery(sites, tmp_path / "b.png"),
            F.plot_prioritization(prio, tmp_path / "c.png"),
            F.plot_grid_strategies(grid, tmp_path / "d.png"),
            F.plot_performance(perf, tmp_path / "e.png")]
    assert all(o.exists() and o.stat().st_size > 5000 for o in outs)
