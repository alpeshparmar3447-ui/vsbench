"""Retrieval metrics for virtual-screening benchmarks.

All functions take a binary label vector ``y`` (1 = active) and a score vector
``s`` where HIGHER means more likely active. Ties are handled explicitly:
rank-based metrics use mid-ranks and enrichment uses the exact expectation
under random tie-breaking, because docking energies are reported to two
decimals and ties are common.
"""
from __future__ import annotations

from typing import Callable, Mapping

import numpy as np
from scipy.stats import rankdata


def to_higher_better(energies) -> np.ndarray:
    """Convert docking energies (lower = better) to scores (higher = better)."""
    return -np.asarray(energies, dtype=float)


def _check(y, s):
    y = np.asarray(y)
    s = np.asarray(s, dtype=float)
    if y.ndim != 1 or y.shape != s.shape:
        raise ValueError("y and s must be 1-D arrays of equal length")
    if not np.isfinite(s).all():
        raise ValueError("scores must be finite; impute failed dockings before scoring")
    if not np.isin(y, (0, 1)).all():
        raise ValueError("labels must be 0 or 1")
    y = y.astype(int)
    n_pos = int(y.sum())
    n_neg = int(len(y) - n_pos)
    if n_pos == 0 or n_neg == 0:
        raise ValueError("both classes must be present")
    return y, s, n_pos, n_neg


def roc_auc(y, s) -> float:
    """Area under the ROC curve via the Mann-Whitney rank-sum identity."""
    y, s, n_pos, n_neg = _check(y, s)
    ranks = rankdata(s)  # average ranks, ascending
    return float((ranks[y == 1].sum() - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def average_precision(y, s) -> float:
    """Average precision (step-wise PR-AUC, the scikit-learn definition)."""
    y, s, n_pos, _ = _check(y, s)
    order = np.argsort(-s, kind="mergesort")
    s_sorted, y_sorted = s[order], y[order]
    tp = np.cumsum(y_sorted)
    fp = np.cumsum(1 - y_sorted)
    last_of_group = np.r_[np.flatnonzero(np.diff(s_sorted)), len(s) - 1]
    tp, fp = tp[last_of_group], fp[last_of_group]
    precision = tp / (tp + fp)
    recall = tp / n_pos
    recall_prev = np.r_[0.0, recall[:-1]]
    return float(np.sum((recall - recall_prev) * precision))


def n_top(n_total: int, fraction: float) -> int:
    """Number of compounds in the top ``fraction`` of a list (ceil, at least 1)."""
    if not 0 < fraction <= 1:
        raise ValueError("fraction must be in (0, 1]")
    return max(1, int(np.ceil(fraction * n_total - 1e-9)))


def _expected_hits(y: np.ndarray, s: np.ndarray, k: int) -> float:
    """Expected number of actives in the top k, with random tie-breaking."""
    if k >= len(s):
        return float(y.sum())
    order = np.argsort(-s, kind="mergesort")
    s_sorted, y_sorted = s[order], y[order]
    cutoff = s_sorted[k - 1]
    above = s_sorted > cutoff
    tied = s_sorted == cutoff
    need = k - int(above.sum())
    return float(y_sorted[above].sum() + y_sorted[tied].sum() * need / int(tied.sum()))


def enrichment_factor(y, s, fraction: float = 0.01) -> float:
    """EF at the given fraction: hit rate in the top slice over the base rate."""
    y, s, n_pos, _ = _check(y, s)
    k = n_top(len(y), fraction)
    return float((_expected_hits(y, s, k) / k) / (n_pos / len(y)))


def max_enrichment_factor(n_total: int, n_pos: int, fraction: float = 0.01) -> float:
    """Largest EF attainable for this list length, active count and fraction."""
    k = n_top(n_total, fraction)
    return float((min(n_pos, k) / k) / (n_pos / n_total))


def normalized_enrichment_factor(y, s, fraction: float = 0.01) -> float:
    """EF divided by its attainable maximum, so lists of different size compare."""
    y, s, n_pos, _ = _check(y, s)
    return enrichment_factor(y, s, fraction) / max_enrichment_factor(len(y), n_pos, fraction)


def bedroc(y, s, alpha: float = 20.0) -> float:
    """BEDROC (Truchon and Bayly 2007) using mid-ranks for tied scores."""
    y, s, n_pos, _ = _check(y, s)
    n = len(y)
    ranks = rankdata(-s, method="average")[y == 1]  # 1 = best
    ra = n_pos / n
    rie = (np.sum(np.exp(-alpha * ranks / n)) / n_pos) / (
        (1.0 / n) * (1.0 - np.exp(-alpha)) / (np.exp(alpha / n) - 1.0)
    )
    return float(
        rie * ra * np.sinh(alpha / 2.0) / (np.cosh(alpha / 2.0) - np.cosh(alpha / 2.0 - alpha * ra))
        + 1.0 / (1.0 - np.exp(alpha * (1.0 - ra)))
    )


METRICS: Mapping[str, Callable] = {
    "roc_auc": roc_auc,
    "pr_auc": average_precision,
    "bedroc20": lambda y, s: bedroc(y, s, 20.0),
    "ef1": lambda y, s: enrichment_factor(y, s, 0.01),
    "ef5": lambda y, s: enrichment_factor(y, s, 0.05),
    "nef1": lambda y, s: normalized_enrichment_factor(y, s, 0.01),
}


def _stratified_indices(y: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    pos = np.flatnonzero(y == 1)
    neg = np.flatnonzero(y == 0)
    return np.r_[rng.choice(pos, len(pos)), rng.choice(neg, len(neg))]


def bootstrap_ci(metric: Callable, y, s, n_boot: int = 1000, seed: int = 0, ci: float = 0.95):
    """Stratified percentile bootstrap. Returns (point, low, high).

    Actives and inactives are resampled separately so every replicate keeps the
    original class counts. The interval reflects compound sampling only.
    """
    y, s, _, _ = _check(y, s)
    rng = np.random.default_rng(seed)
    vals = np.empty(n_boot)
    for b in range(n_boot):
        idx = _stratified_indices(y, rng)
        vals[b] = metric(y[idx], s[idx])
    lo, hi = np.quantile(vals, [(1 - ci) / 2, 1 - (1 - ci) / 2])
    return float(metric(y, s)), float(lo), float(hi)


def paired_bootstrap_draws(metric: Callable, y, s_a, s_b, n_boot: int = 1000, seed: int = 0):
    """Paired difference metric(s_a) - metric(s_b), with the bootstrap draws behind it.

    The draws are returned rather than summarised so that several confidence levels can be read
    off one run: a multiplicity-adjusted interval then costs nothing beyond the pointwise one.
    """
    y, s_a, _, _ = _check(y, s_a)
    _, s_b, _, _ = _check(y, s_b)
    rng = np.random.default_rng(seed)
    diffs = np.empty(n_boot)
    for b in range(n_boot):
        idx = _stratified_indices(y, rng)
        diffs[b] = metric(y[idx], s_a[idx]) - metric(y[idx], s_b[idx])
    return float(metric(y, s_a) - metric(y, s_b)), diffs


def interval_from_draws(draws, ci: float = 0.95, n_comparisons: int = 1):
    """Percentile interval, Bonferroni-widened for ``n_comparisons`` simultaneous comparisons.

    One comparison leaves the interval unchanged. For m comparisons the level becomes
    1 - (1-ci)/m, which targets approximately ``ci`` simultaneous coverage across that family.
    Comparing many variants against one primary and reading only the pointwise intervals
    otherwise finds a winner at the nominal rate by chance alone.
    """
    if n_comparisons < 1:
        raise ValueError("n_comparisons must be at least 1")
    draws = np.asarray(draws, dtype=float)
    alpha = (1.0 - ci) / n_comparisons
    lo, hi = np.quantile(draws, [alpha / 2, 1 - alpha / 2])
    return float(lo), float(hi)


def draws_resolve_interval(n_boot: int, ci: float = 0.95, n_comparisons: int = 1,
                           min_draws: float = 5.0) -> bool:
    """Whether ``n_boot`` draws put enough mass in each tail for that interval to mean anything.

    An adjusted interval reads a far quantile: at 95% over 12 comparisons each tail holds
    0.2% of the draws, which is two of 1000. The endpoint is then a couple of extreme replicates
    rather than an estimate, so a caller that needs the adjusted interval needs more draws.
    """
    return n_boot * (1.0 - ci) / n_comparisons / 2.0 >= min_draws


def paired_bootstrap_difference(metric: Callable, y, s_a, s_b, n_boot: int = 1000,
                                seed: int = 0, ci: float = 0.95):
    """Bootstrap the paired difference metric(s_a) - metric(s_b).

    Both score vectors are evaluated on the same resampled compounds. Returns
    (difference, low, high, two_sided_p) where p is the bootstrap tail mass on
    the side of zero opposite to the point estimate, doubled and capped at 1.
    """
    point, diffs = paired_bootstrap_draws(metric, y, s_a, s_b, n_boot=n_boot, seed=seed)
    lo, hi = interval_from_draws(diffs, ci)
    tail = float(np.mean(diffs <= 0)) if point > 0 else float(np.mean(diffs >= 0))
    return point, lo, hi, float(min(1.0, 2 * tail))


def consensus_zscore(scores: Mapping[str, np.ndarray]) -> np.ndarray:
    """Average of per-engine z-scores (higher = better).

    NaN entries (failed dockings) receive that engine's worst finite z-score, so
    a failure is treated as a non-binder and cannot be silently dropped.
    """
    cols = []
    for name, v in scores.items():
        v = np.asarray(v, dtype=float)
        finite = np.isfinite(v)
        if finite.sum() < 2:
            raise ValueError(f"engine {name!r} has fewer than two finite scores")
        z = np.full(v.shape, np.nan)
        sd = v[finite].std()
        z[finite] = (v[finite] - v[finite].mean()) / (sd if sd > 0 else 1.0)
        z[~finite] = np.nanmin(z)
        cols.append(z)
    return np.mean(cols, axis=0)
