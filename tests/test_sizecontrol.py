import numpy as np
import pandas as pd
import pytest

from vsbench import sizecontrol as S


def test_one_stratum_reproduces_the_ordinary_auc():
    rng = np.random.default_rng(0)
    y = np.r_[np.ones(20, int), np.zeros(20, int)]
    s = rng.normal(size=40) + y * 0.8
    from vsbench.metrics import roc_auc
    pooled, per, _, _ = S.conditional_auc(y, s, np.zeros(40, int))
    assert pooled == pytest.approx(roc_auc(y, s))
    assert len(per) == 1 and per.loc[0, "pairs"] == 400


def test_a_score_that_is_only_size_loses_all_discrimination_when_size_is_controlled():
    """Actives are large, inactives small, and the score IS the size. Unstratified AUC is perfect;
    within a size class the score cannot tell them apart, so the conditional AUC must collapse."""
    size = np.r_[np.arange(30, 50), np.arange(10, 30)]           # actives 30-49, inactives 10-29
    y = np.r_[np.ones(20, int), np.zeros(20, int)]
    from vsbench.metrics import roc_auc
    assert roc_auc(y, size.astype(float)) == pytest.approx(1.0)
    strata = S.strata_from_quantiles(size, 4)
    pooled, per, _, _ = S.conditional_auc(y, size.astype(float), strata)
    assert np.isnan(pooled)                                      # no stratum holds both classes
    assert (per["pairs"] == 0).all()


def test_an_overlapping_size_confound_is_reduced_not_erased():
    rng = np.random.default_rng(1)
    size = np.r_[rng.integers(20, 40, 60), rng.integers(15, 35, 60)]
    y = np.r_[np.ones(60, int), np.zeros(60, int)]
    s = size + rng.normal(0, 1.0, 120)                           # mostly size, a little noise
    from vsbench.metrics import roc_auc
    plain = roc_auc(y, s.astype(float))
    pooled, per, lo, hi = S.conditional_auc(y, s.astype(float), S.strata_from_quantiles(size, 4),
                                            n_boot=200, seed=0)
    assert plain > pooled                                        # the confound inflated the plain value
    assert 0.3 < pooled < 0.75 and lo <= pooled <= hi
    assert (per["pairs"] > 0).sum() >= 3


def test_an_empty_stratum_is_reported_but_contributes_no_pairs():
    y = np.array([1, 1, 0, 0, 1, 1])
    s = np.array([5.0, 4.0, 1.0, 2.0, 9.0, 8.0])
    strata = np.array([0, 0, 0, 0, 1, 1])                        # stratum 1 has no inactive
    pooled, per, _, _ = S.conditional_auc(y, s, strata)
    assert np.isnan(per.set_index("stratum").loc[1, "roc_auc"])
    assert per.set_index("stratum").loc[1, "pairs"] == 0
    assert pooled == pytest.approx(1.0)                          # stratum 0 alone, actives on top


def test_ties_do_not_split_across_strata():
    x = np.r_[np.full(50, 10.0), np.full(50, 20.0)]
    st = S.strata_from_quantiles(x, 4)
    assert len(np.unique(st)) == 2
    for v in (10.0, 20.0):
        assert len(np.unique(st[x == v])) == 1


def test_non_finite_sizes_are_refused():
    with pytest.raises(ValueError):
        S.strata_from_quantiles([1.0, 2.0, np.nan])


def test_table_pairs_scores_on_one_set_of_strata_and_reports_the_shift():
    rng = np.random.default_rng(3)
    n = 80
    df = pd.DataFrame({"label": np.r_[np.ones(20, int), np.zeros(60, int)],
                       "n_heavy": rng.integers(15, 45, n).astype(float),
                       "good": rng.normal(size=n), "sizey": 0.0})
    df["good"] = df["good"] + df["label"] * 1.5
    df["sizey"] = df["n_heavy"] + rng.normal(0, 0.5, n)
    summary, detail = S.table(df, "random", {"good": "good", "sizey": "sizey"}, n_boot=100)
    assert list(summary.score) == ["good", "sizey"]
    assert (summary["shift"] == summary["conditional_auc"] - summary["roc_auc"]).all()
    assert set(detail.score) == {"good", "sizey"}
    g = detail[detail.score == "good"].set_index("stratum")["pairs"]
    z = detail[detail.score == "sizey"].set_index("stratum")["pairs"]
    pd.testing.assert_series_equal(g, z, check_names=False)      # same strata for every score


def test_size_profile_reports_the_gap_that_drives_the_shift():
    df = pd.DataFrame({"label": [1, 1, 0, 0], "n_heavy": [40.0, 42.0, 20.0, 22.0]})
    p = S.size_profile(df, "matched")
    assert p["median_active"] == 41 and p["median_inactive"] == 21 and p["median_gap"] == 20
