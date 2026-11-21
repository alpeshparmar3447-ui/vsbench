"""Conditional (size-controlled) AUC: does a score separate actives from inactives of their own size?

Docking energies grow with molecular size, so an AUC measured against a negative set that is smaller on
average is partly a size measurement. The statistic here is the probability that an active outranks an
inactive **in its own heavy-atom stratum**. Restricting the comparison to same-stratum pairs and pooling the
per-stratum AUCs by the number of pairs each contributes are the same operation, so this is one number, not
two, and it reduces to the ordinary AUC when every compound falls in one stratum.

A stratum with no active or no inactive contributes no pairs; it is reported with an empty AUC rather than
counted as chance, which would pull the pooled value towards 0.5 for free.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import metrics as M

N_STRATA = 4


def strata_from_quantiles(x, n_strata: int = N_STRATA) -> np.ndarray:
    """Label each value by which quantile bin of ``x`` it falls in (0 .. n-1).

    Ties are kept together: a heavy-atom count that spans a quantile edge lands wholly in one stratum, so
    the number of distinct strata can be smaller than requested.
    """
    x = np.asarray(x, dtype=float)
    if not np.isfinite(x).all():
        raise ValueError("strata require finite values")
    edges = np.unique(np.quantile(x, np.linspace(0, 1, n_strata + 1)[1:-1]))
    return np.searchsorted(edges, x, side="left")


def conditional_auc(y, s, strata, n_boot: int = 0, seed: int = 0, ci: float = 0.95):
    """Pair-weighted mean of the within-stratum AUCs, plus the per-stratum detail.

    Returns ``(pooled, per_stratum_df, low, high)``. ``low``/``high`` are NaN unless ``n_boot`` is given;
    the bootstrap resamples actives and inactives separately, as elsewhere in this project.
    """
    y = np.asarray(y, dtype=int)
    s = np.asarray(s, dtype=float)
    strata = np.asarray(strata)
    if not (len(y) == len(s) == len(strata)):
        raise ValueError("y, s and strata must have the same length")

    rows = []
    for k in np.unique(strata):
        m = strata == k
        n_pos, n_neg = int((y[m] == 1).sum()), int((y[m] == 0).sum())
        pairs = n_pos * n_neg
        auc = M.roc_auc(y[m], s[m]) if pairs else np.nan
        rows.append({"stratum": int(k), "n_active": n_pos, "n_inactive": n_neg,
                     "pairs": pairs, "roc_auc": auc})
    per = pd.DataFrame(rows)
    pooled = _pool(per)

    lo = hi = np.nan
    if n_boot:
        rng = np.random.default_rng(seed)
        pos, neg = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
        vals = np.empty(n_boot)
        for b in range(n_boot):
            idx = np.r_[rng.choice(pos, len(pos)), rng.choice(neg, len(neg))]
            vals[b] = _pool(_per_stratum(y[idx], s[idx], strata[idx]))
        good = vals[np.isfinite(vals)]
        if len(good):
            lo, hi = np.quantile(good, [(1 - ci) / 2, 1 - (1 - ci) / 2])
    return pooled, per, float(lo), float(hi)


def _per_stratum(y, s, strata) -> pd.DataFrame:
    rows = []
    for k in np.unique(strata):
        m = strata == k
        n_pos, n_neg = int((y[m] == 1).sum()), int((y[m] == 0).sum())
        pairs = n_pos * n_neg
        rows.append({"pairs": pairs, "roc_auc": M.roc_auc(y[m], s[m]) if pairs else np.nan})
    return pd.DataFrame(rows)


def _pool(per: pd.DataFrame) -> float:
    use = per[(per["pairs"] > 0) & per["roc_auc"].notna()]
    if not len(use):
        return float("nan")
    return float(np.average(use["roc_auc"], weights=use["pairs"]))


def table(df: pd.DataFrame, negative_set: str, score_columns: dict, size_column: str = "n_heavy",
          n_strata: int = N_STRATA, n_boot: int = 1000, seed: int = 0) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Unstratified and conditional AUC of every score, for one evaluation set.

    Returns ``(summary, per_stratum)``. The strata come from the pooled size distribution of this
    evaluation set, so they are the same for every score and the comparison between scores is paired.
    """
    y = df["label"].to_numpy(int)
    strata = strata_from_quantiles(df[size_column].to_numpy(float), n_strata)
    summary, detail = [], []
    for name, col in score_columns.items():
        s = df[col].to_numpy(float)
        plain, p_lo, p_hi = M.bootstrap_ci(M.roc_auc, y, s, n_boot=n_boot, seed=seed)
        cond, per, c_lo, c_hi = conditional_auc(y, s, strata, n_boot=n_boot, seed=seed)
        summary.append({"score": name, "negative_set": negative_set, "roc_auc": plain,
                        "ci_low": p_lo, "ci_high": p_hi, "conditional_auc": cond,
                        "cond_ci_low": c_lo, "cond_ci_high": c_hi, "shift": cond - plain,
                        "n_strata": int(per["pairs"].gt(0).sum())})
        detail.append(per.assign(score=name, negative_set=negative_set))
    return pd.DataFrame(summary), pd.concat(detail, ignore_index=True)


def size_profile(df: pd.DataFrame, negative_set: str, size_column: str = "n_heavy") -> dict:
    """Median size of actives and of inactives, and the gap. Explains the direction of any shift."""
    a = df.loc[df["label"] == 1, size_column].to_numpy(float)
    n = df.loc[df["label"] == 0, size_column].to_numpy(float)
    return {"negative_set": negative_set, "median_active": float(np.median(a)),
            "median_inactive": float(np.median(n)),
            "median_gap": float(np.median(a) - np.median(n)),
            "n_active": int(len(a)), "n_inactive": int(len(n))}
