import numpy as np
import pandas as pd
import pytest

from vsbench import analysis as A
from vsbench.pipeline import merge_seed_runs, seed_study_subset


def _frame(energies, ok=None, largest=None):
    n = len(energies)
    return pd.DataFrame({"name": [f"c{i}" for i in range(n)], "ad_ok": ok if ok is not None else [True] * n,
                         "ad_energy": energies, "ad_largest_energy": largest if largest is not None else energies,
                         "ad_run": 1, "ad_n_clusters": 2, "ad_largest_size": 5, "ad_largest_run": 1})


def test_merge_seed_runs_takes_the_mean_and_keeps_per_seed_columns():
    m = merge_seed_runs({42: _frame([-8.0, -6.0, np.nan]), 43: _frame([-9.0, -6.5, -5.0], ok=[True, True, True])})
    assert m["ad_energy"].tolist()[:2] == [-8.5, -6.25]
    assert m.loc[2, "ad_energy"] == -5.0                                  # a seed that failed does not erase the compound
    assert m.loc[0, "ad_energy_s42"] == -8.0 and m.loc[0, "ad_energy_s43"] == -9.0
    assert m.loc[0, "ad_energy_sd"] == pytest.approx(np.std([-8.0, -9.0], ddof=1))
    assert bool(m.loc[2, "ad_ok"])


def test_a_single_seed_is_passed_through_without_a_spread():
    m = merge_seed_runs({42: _frame([-8.0, -6.0])})
    assert m["ad_energy"].tolist() == [-8.0, -6.0] and m["ad_energy_sd"].isna().all()


def test_a_compound_that_failed_in_every_seed_stays_failed():
    f = _frame([np.nan, -5.0], ok=[False, True])
    m = merge_seed_runs({1: f, 2: f})
    assert np.isnan(m.loc[0, "ad_energy"]) and not bool(m.loc[0, "ad_ok"]) and bool(m.loc[1, "ad_ok"])


def _seeded_table(noise, n=120, n_act=20, seed=0):
    rng = np.random.default_rng(seed)
    y = np.r_[np.ones(n_act), np.zeros(n - n_act)].astype(int)
    truth = -(7 + 1.0 * y + rng.normal(size=n))
    df = pd.DataFrame({"label": y})
    for s in (42, 43, 44):
        df[f"ad_energy_s{s}"] = truth + rng.normal(0, noise, n)
    return df


def test_seed_stability_is_high_for_a_stable_score_and_low_for_a_noisy_one():
    stable = A.seed_stability(_seeded_table(0.02), [42, 43, 44])
    noisy = A.seed_stability(_seeded_table(3.0), [42, 43, 44])
    assert stable["mean_spearman_between_seeds"] > 0.99 and stable["stable_top_n_fraction"] > 0.8
    assert noisy["mean_spearman_between_seeds"] < 0.4 and noisy["stable_top_n_fraction"] < 0.4
    assert stable["auc_range"] < noisy["auc_range"] + 0.2 and stable["mean_per_compound_sd"] < noisy["mean_per_compound_sd"]
    assert stable["n_seeds"] == 3


def test_gnina_seed_study_uses_only_the_redocked_subset():
    rng = np.random.default_rng(1)
    n = 80
    df = pd.DataFrame({"label": np.r_[np.ones(20), np.zeros(60)].astype(int)})
    df["gn_cnn_vs"] = rng.normal(size=n) + 0.8 * df.label
    for s in (43, 44):
        df[f"gn_cnn_vs_s{s}"] = np.where(np.arange(n) < 50, df.gn_cnn_vs + rng.normal(0, 0.1, n), np.nan)
    r = A.gnina_seed_study(df, [43, 44])
    assert r["n"] == 50 and r["mean_spearman_between_seeds"] > 0.95 and r["auc_range"] < 0.1
    assert A.gnina_seed_study(df.assign(label=0), [43, 44]) == {}


def _subset_table(n_act, n_neg):
    names = [f"a{i}" for i in range(n_act)] + [f"n{i}" for i in range(n_neg)]
    return pd.DataFrame({"name": names, "label": [1] * n_act + [0] * n_neg})


def _old_rule(table, n, seed):
    """The rule in use up to 2026-11-20, kept so the registered arm's draw stays pinned."""
    return set(table.loc[table.label == 1, "name"]) | set(
        table.loc[table.label == 0].sample(max(0, n - int((table.label == 1).sum())),
                                           random_state=seed)["name"])


@pytest.mark.parametrize("n_act", [27, 40])
def test_the_seed_study_subset_is_unchanged_when_the_actives_fit_in_the_budget(n_act):
    """The registered arm capped actives at 40 against a budget of 60, so every target had room for
    negatives. Those draws are published, and this fix must not move them."""
    table = _subset_table(n_act, 940)
    assert seed_study_subset(table, 60, 42) == _old_rule(table, 60, 42)


@pytest.mark.parametrize("n_act", [60, 79, 102, 308])
def test_an_actives_heavy_target_still_gets_negatives(n_act):
    """More actives than the budget used to leave zero negatives, so the subset was single-class and
    the per-seed ROC-AUC was undefined: four of five targets in the 1000-active arm produced no study
    at all. The subset now stays inside the budget and keeps both classes."""
    table = _subset_table(n_act, 1900)
    keep = seed_study_subset(table, 60, 42)
    labels = table.set_index("name").label
    assert len(keep) == 60
    assert labels[list(keep)].nunique() == 2
    assert sum(labels[x] == 0 for x in keep) >= 20
    assert _old_rule(table, 60, 42) != keep                      # the case this fix exists for


def test_the_negative_floor_scales_with_a_small_budget():
    """The end-to-end test runs a budget of 12, so a floor fixed at the registered arm's 20 would
    swamp it and leave one active."""
    keep = seed_study_subset(_subset_table(30, 30), 12, 42)
    labels = _subset_table(30, 30).set_index("name").label
    assert len(keep) == 12 and labels[list(keep)].nunique() == 2
    assert sum(labels[x] == 1 for x in keep) == 8 and sum(labels[x] == 0 for x in keep) == 4


def test_a_target_with_no_negatives_yields_actives_only_rather_than_raising():
    """Nothing can rescue a target that has one class, but it must degrade to a reportable subset the
    caller can warn about, not an exception or a silent empty set."""
    table = _subset_table(50, 0)
    keep = seed_study_subset(table, 60, 42)
    assert len(keep) == 40 and all(k.startswith("a") for k in keep)


def test_the_subset_is_reproducible_for_a_fixed_seed_and_moves_with_it():
    table = _subset_table(102, 1900)
    assert seed_study_subset(table, 60, 42) == seed_study_subset(table, 60, 42)
    assert seed_study_subset(table, 60, 42) != seed_study_subset(table, 60, 7)


def test_aggregate_covers_every_analysed_target_not_only_the_ones_this_run_touched(tmp_path):
    """`--targets X` sets cfg.targets to X alone, and the panel files were written from that list, so
    re-running one target after an interruption truncated `metrics_all_targets.csv` to it. report.py
    reads that file as its metrics source, so the loss was silent."""
    from vsbench.pipeline import Config, aggregate

    for t, auc in (("ALPHA", 0.70), ("BETA", 0.80), ("GAMMA", 0.90)):
        d = tmp_path / t / "analysis"
        d.mkdir(parents=True)
        pd.DataFrame([{"score": "gnina", "negative_set": "random", "metric": "roc_auc",
                       "value": auc}]).to_csv(d / "metrics.csv", index=False)

    cfg = Config(data_root=str(tmp_path), results_dir=str(tmp_path), targets=["BETA"])
    macro = aggregate(cfg)

    allm = pd.read_csv(tmp_path / "metrics_all_targets.csv")
    assert sorted(allm["target"]) == ["ALPHA", "BETA", "GAMMA"]       # not just BETA
    assert macro["n_targets"].iloc[0] == 3
    assert macro["value"].iloc[0] == pytest.approx(0.80)
