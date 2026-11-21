import itertools

import numpy as np
import pytest
from rdkit.ML.Scoring.Scoring import CalcBEDROC
from sklearn.metrics import average_precision_score, roc_auc_score

from vsbench import metrics as M


def _data(n=400, n_pos=40, shift=0.8, seed=1, decimals=None):
    rng = np.random.default_rng(seed)
    y = np.r_[np.ones(n_pos), np.zeros(n - n_pos)].astype(int)
    s = rng.normal(size=n) + shift * y
    if decimals is not None:
        s = np.round(s, decimals)
    return y, s


@pytest.mark.parametrize("decimals", [None, 1, 0])
def test_roc_auc_matches_sklearn_including_ties(decimals):
    y, s = _data(decimals=decimals)
    assert M.roc_auc(y, s) == pytest.approx(roc_auc_score(y, s), abs=1e-12)


@pytest.mark.parametrize("decimals", [None, 1, 0])
def test_average_precision_matches_sklearn_including_ties(decimals):
    y, s = _data(decimals=decimals, seed=5)
    assert M.average_precision(y, s) == pytest.approx(average_precision_score(y, s), abs=1e-12)


def test_perfect_and_inverted_ranking():
    y = np.array([1, 1, 0, 0, 0, 0])
    s = np.array([5, 4, 3, 2, 1, 0], float)
    assert M.roc_auc(y, s) == 1.0
    assert M.roc_auc(y, -s) == 0.0
    assert M.average_precision(y, s) == 1.0
    assert M.normalized_enrichment_factor(y, s, 0.34) == pytest.approx(1.0)


def test_constant_scores_give_chance_values():
    y = np.array([1, 0, 0, 0, 1, 0, 0, 0, 0, 0])
    s = np.zeros(10)
    assert M.roc_auc(y, s) == pytest.approx(0.5)
    assert M.enrichment_factor(y, s, 0.2) == pytest.approx(1.0)


def test_bedroc_matches_rdkit_when_there_are_no_ties():
    y, s = _data(seed=3)
    order = np.argsort(-s)
    ranked = [[int(y[i])] for i in order]
    for alpha in (20.0, 80.5):
        assert M.bedroc(y, s, alpha) == pytest.approx(CalcBEDROC(ranked, 0, alpha), abs=1e-9)


def test_bedroc_extremes():
    n, n_pos = 1000, 20
    y = np.r_[np.ones(n_pos), np.zeros(n - n_pos)].astype(int)
    best = -np.arange(n, dtype=float)
    assert M.bedroc(y, best, 20.0) > 0.99
    assert M.bedroc(y, -best, 20.0) < 0.01


def test_enrichment_matches_brute_force_tie_breaking():
    # Exact expectation over every permutation of tied compounds.
    y = np.array([1, 0, 1, 0, 0, 1, 0, 0])
    s = np.array([2, 2, 2, 1, 1, 1, 0, 0], float)
    k = M.n_top(len(y), 0.375)  # top 3
    hits = []
    for perm in itertools.permutations(range(len(y))):
        order = sorted(perm, key=lambda i: -s[i])  # stable: perm decides tie order
        hits.append(y[order[:k]].sum())
    expected = np.mean(hits)
    assert M._expected_hits(y, s, k) == pytest.approx(expected)
    ef = M.enrichment_factor(y, s, 0.375)
    assert ef == pytest.approx((expected / k) / (y.sum() / len(y)))


def test_n_top_rounds_up_and_floors_at_one():
    assert M.n_top(1050, 0.01) == 11
    assert M.n_top(100, 0.01) == 1
    assert M.n_top(10, 0.01) == 1
    assert M.n_top(1000, 0.05) == 50


def test_nef_is_bounded_by_one():
    for seed in range(5):
        y, s = _data(seed=seed)
        v = M.normalized_enrichment_factor(y, s, 0.01)
        assert 0.0 <= v <= 1.0 + 1e-12


def test_bootstrap_is_reproducible_and_brackets_the_point_estimate():
    y, s = _data()
    a = M.bootstrap_ci(M.roc_auc, y, s, n_boot=300, seed=7)
    b = M.bootstrap_ci(M.roc_auc, y, s, n_boot=300, seed=7)
    assert a == b
    point, lo, hi = a
    assert lo < point < hi


def test_bootstrap_interval_narrows_with_more_data():
    small = M.bootstrap_ci(M.roc_auc, *_data(n=200, n_pos=20), n_boot=300, seed=1)
    large = M.bootstrap_ci(M.roc_auc, *_data(n=2000, n_pos=200), n_boot=300, seed=1)
    assert (large[2] - large[1]) < (small[2] - small[1])


def test_bootstrap_coverage_for_a_known_auc():
    # Two unit-variance normals separated by d have AUC = Phi(d / sqrt(2)).
    from scipy.stats import norm
    true_auc = norm.cdf(0.8 / np.sqrt(2))
    covered = 0
    trials = 60
    for t in range(trials):
        y, s = _data(n=300, n_pos=60, shift=0.8, seed=100 + t)
        _, lo, hi = M.bootstrap_ci(M.roc_auc, y, s, n_boot=200, seed=t)
        covered += lo <= true_auc <= hi
    assert covered / trials >= 0.85


def test_paired_difference_identical_scores_and_clear_winner():
    y, s = _data()
    d, lo, hi, p = M.paired_bootstrap_difference(M.roc_auc, y, s, s, n_boot=200, seed=0)
    assert d == 0.0 and lo == 0.0 and hi == 0.0 and p == 1.0
    rng = np.random.default_rng(0)
    noise = rng.normal(size=len(y))
    d, lo, hi, p = M.paired_bootstrap_difference(M.roc_auc, y, s, noise, n_boot=300, seed=0)
    assert d > 0.15 and lo > 0 and p < 0.01


def test_consensus_zscore_fills_failures_with_worst_and_ignores_scale():
    a = np.array([1.0, 2.0, 3.0, np.nan])
    b = np.array([10.0, 20.0, 30.0, 40.0])
    c = M.consensus_zscore({"a": a, "b": b})
    assert c[2] > c[1] > c[0]
    assert c[3] == pytest.approx((np.nanmin((a - np.nanmean(a)) / np.nanstd(a)) + ((40 - 25) / b.std())) / 2)
    c2 = M.consensus_zscore({"a": a * 7 + 3, "b": b})
    np.testing.assert_allclose(c, c2)
    with pytest.raises(ValueError):
        M.consensus_zscore({"x": np.array([1.0, np.nan, np.nan])})


@pytest.mark.parametrize("bad", [
    (np.array([1, 1, 1]), np.array([1.0, 2.0, 3.0])),
    (np.array([0, 0, 0]), np.array([1.0, 2.0, 3.0])),
    (np.array([0, 1, 2]), np.array([1.0, 2.0, 3.0])),
    (np.array([0, 1]), np.array([1.0, np.nan])),
    (np.array([0, 1, 0]), np.array([1.0, 2.0])),
])
def test_invalid_inputs_raise(bad):
    with pytest.raises(ValueError):
        M.roc_auc(*bad)


def test_to_higher_better():
    np.testing.assert_array_equal(M.to_higher_better([-9.5, -7.0]), [9.5, 7.0])


# --------------------------------------------------- multiplicity-adjusted intervals
def test_the_adjusted_interval_contains_the_pointwise_one():
    """Bonferroni widens; it can never declare a difference the pointwise interval missed."""
    rng = np.random.default_rng(0)
    draws = rng.normal(0.05, 0.1, 4000)
    lo, hi = M.interval_from_draws(draws)
    alo, ahi = M.interval_from_draws(draws, n_comparisons=12)
    assert alo <= lo and ahi >= hi


def test_one_comparison_leaves_the_interval_unchanged():
    draws = np.linspace(-1, 1, 2001)
    assert M.interval_from_draws(draws) == M.interval_from_draws(draws, n_comparisons=1)


def test_a_non_positive_comparison_count_is_refused():
    with pytest.raises(ValueError):
        M.interval_from_draws(np.zeros(10), n_comparisons=0)


def test_the_draws_reproduce_the_summarised_difference():
    """paired_bootstrap_difference must stay a view on the same draws, not a second run."""
    y = np.array([1] * 20 + [0] * 40)
    rng = np.random.default_rng(1)
    s_a, s_b = rng.normal(size=60), rng.normal(size=60)
    d, lo, hi, _ = M.paired_bootstrap_difference(M.roc_auc, y, s_a, s_b, n_boot=300, seed=7)
    d2, draws = M.paired_bootstrap_draws(M.roc_auc, y, s_a, s_b, n_boot=300, seed=7)
    assert d == pytest.approx(d2)
    assert (lo, hi) == pytest.approx(M.interval_from_draws(draws))


def test_draws_resolve_interval_knows_when_a_tail_is_too_thin():
    """1000 draws put ~2 in each tail of a 12-way adjusted interval: not an estimate."""
    assert not M.draws_resolve_interval(1000, n_comparisons=12)
    assert M.draws_resolve_interval(5000, n_comparisons=12)
    assert M.draws_resolve_interval(1000)                      # unadjusted, 25 per tail


# ------------------------------------------------- docking against the best baseline


def _baseline_frame(n=120, sep=1.2, seed=3):
    """Actives separated from inactives on nn_sim but not on the docking score."""
    import numpy as np
    import pandas as pd
    rng = np.random.default_rng(seed)
    y = np.array([1] * (n // 3) + [0] * (n - n // 3))
    na, ni = int(y.sum()), int(len(y) - y.sum())
    return pd.DataFrame({
        "label": y,
        "neg_set": ["active"] * na + ["matched"] * ni,
        "nn_sim": np.r_[rng.normal(sep, 1.0, na), rng.normal(0.0, 1.0, ni)],
        "n_heavy": rng.normal(0.0, 1.0, len(y)),
        "mw": rng.normal(0.0, 1.0, len(y)),
        "clogp": rng.normal(0.0, 1.0, len(y)),
        "tpsa": rng.normal(0.0, 1.0, len(y)),
        "hba": rng.normal(0.0, 1.0, len(y)),
        "hbd": rng.normal(0.0, 1.0, len(y)),
        "rotb": rng.normal(0.0, 1.0, len(y)),
        "ad_lowest": rng.normal(0.0, 1.0, len(y)),
        "ad_largest": rng.normal(0.0, 1.0, len(y)),
        "gn_cnn_vs": rng.normal(0.0, 1.0, len(y)),
        "gn_vina_hb": rng.normal(0.0, 1.0, len(y)),
        "consensus": rng.normal(0.0, 1.0, len(y)),
    })


def test_baseline_comparison_picks_the_strongest_baseline_and_reports_docking_below_it():
    from vsbench import analysis
    out = analysis.baseline_comparison(_baseline_frame(), "matched", n_boot=400, seed=0)
    assert len(out) == len(analysis.SCORE_COLUMNS)
    # nn_sim carries the separation, so it must be the chosen baseline for every score
    assert set(out.best_baseline) == {"nn_sim"}
    assert out.n_baselines.eq(len(analysis.BASELINE_COLUMNS)).all()
    # the docking scores are pure noise here, so none may be reported above the baseline
    assert not out.beats_baseline.any()
    assert out.below_baseline.any()


def test_baseline_comparison_adjusted_interval_contains_the_pointwise_one():
    from vsbench import analysis
    out = analysis.baseline_comparison(_baseline_frame(), "matched", n_boot=400, seed=0)
    assert (out.adj_ci_low <= out.ci_low + 1e-12).all()
    assert (out.adj_ci_high >= out.ci_high - 1e-12).all()
    # a score can never be promoted by widening, only demoted
    assert (~out.beats_baseline | (out.ci_low > 0)).all()


def test_baseline_comparison_returns_empty_without_a_usable_split():
    import pandas as pd
    from vsbench import analysis
    d = _baseline_frame()
    assert analysis.baseline_comparison(d[d.label == 1], "matched", n_boot=50).empty
    assert analysis.baseline_comparison(pd.DataFrame({"label": []}), "matched", n_boot=50).empty
