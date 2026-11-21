"""Score variants built from values the screen already recorded: no re-docking, no new poses.

Three things live here.

*Scores that were computed and never evaluated.* gnina writes a CNN pose score and a CNN affinity estimate
alongside CNN_VS; the screen stored both and no table scored them.

*Fusion rules* other than the registered one (the mean of per-engine z-scores): the mean of per-engine
ranks, reciprocal-rank fusion, and a z-score mean weighted by how well each engine reproduced the crystal
pose on the template actually used. Rank fusion is scale-free, so an engine with a long score tail cannot
dominate by spread alone; the weighted variant asks whether trusting the engine with the better redocking
geometry pays off.

*A comparison against the registered consensus* by paired bootstrap on the same compounds.

Per PROTOCOL.md Amendment 3 the registered consensus stays the paper's primary score whatever this finds:
these variants are a secondary comparison, and the best of them is not promoted.
"""
from __future__ import annotations

from collections.abc import Mapping

import numpy as np
import pandas as pd

from . import metrics as M
from .analysis import impute_worst

RRF_K = 60                      # Cormack et al.'s constant; damps the influence of a single engine's top rank

# Scores recorded by the screen that the registered tables do not evaluate.
EXTRA_SCORES = {"gnina (CNN score)": "gn_cnn_score", "gnina (CNN affinity)": "gn_cnn_affinity"}

# Size normalisations of the energy-like scores (PROTOCOL.md Amendment 3, item C2). The energies are
# negative-is-better, so each is negated first and then divided by a positive function of ligand size;
# the result stays higher-is-better.
SIZE_NORMALISED = {
    "AutoDock-GPU per heavy atom": ("ad_lowest", 1.0),
    "AutoDock-GPU per heavy atom^(1/3)": ("ad_lowest", 1.0 / 3.0),
    "gnina (Vina-type) per heavy atom": ("gn_vina_hb", 1.0),
    "gnina (Vina-type) per heavy atom^(1/3)": ("gn_vina_hb", 1.0 / 3.0),
}


def size_normalise(score: np.ndarray, n_heavy: np.ndarray, power: float = 1.0) -> np.ndarray:
    """Divide a higher-is-better score by ``n_heavy ** power``.

    A compound with no recorded heavy-atom count, or none at all, cannot be normalised and keeps the worst
    finite value of the result rather than becoming infinite.
    """
    score = np.asarray(score, dtype=float)
    n = np.asarray(n_heavy, dtype=float)
    out = np.full(len(score), np.nan)
    ok = np.isfinite(n) & (n > 0) & np.isfinite(score)
    out[ok] = score[ok] / np.power(n[ok], power)
    return impute_worst(out)

# The two engine scores the registered consensus is built from, as higher-is-better columns.
ENGINE_COLUMNS = {"ad": "ad_lowest", "gn": "gn_cnn_vs"}

COMPARISON_METRICS = {"roc_auc": M.roc_auc, "ef1": lambda y, s: M.enrichment_factor(y, s, 0.01),
                      "bedroc20": lambda y, s: M.bedroc(y, s, 20.0)}


def _ranks(v: np.ndarray) -> np.ndarray:
    """1-based ranks, best (highest) first, ties averaged."""
    return pd.Series(-np.asarray(v, dtype=float)).rank(method="average").to_numpy()


def rank_mean(scores: Mapping[str, np.ndarray]) -> np.ndarray:
    """Mean of per-engine ranks, negated so that higher is better."""
    if not scores:
        raise ValueError("no scores to fuse")
    return -np.mean([_ranks(v) for v in scores.values()], axis=0)


def reciprocal_rank(scores: Mapping[str, np.ndarray], k: int = RRF_K) -> np.ndarray:
    """Reciprocal-rank fusion: the sum of 1/(k + rank) over engines. Higher is better."""
    if not scores:
        raise ValueError("no scores to fuse")
    if k <= 0:
        raise ValueError("k must be positive")
    return np.sum([1.0 / (k + _ranks(v)) for v in scores.values()], axis=0)


def weighted_zscore(scores: Mapping[str, np.ndarray], weights: Mapping[str, float]) -> np.ndarray:
    """Weighted mean of per-engine z-scores. Weights need not sum to one; they are normalised here."""
    missing = set(scores) - set(weights)
    if missing:
        raise ValueError(f"no weight for {sorted(missing)}")
    w = np.array([float(weights[k]) for k in scores])
    if (w < 0).any() or not w.sum():
        raise ValueError("weights must be non-negative and not all zero")
    z = []
    for v in scores.values():
        v = np.asarray(v, dtype=float)
        sd = v.std()
        z.append((v - v.mean()) / sd if sd else np.zeros_like(v))
    return np.average(np.array(z), axis=0, weights=w / w.sum())


def redocking_weights(ad_rmsd: float, gn_rmsd: float) -> dict[str, float]:
    """Weight each engine by its in-place redocking RMSD on the selected template: 1/(1 + RMSD).

    A failed or missing redocking gets the smallest weight the formula can give rather than being dropped,
    so a target whose template failed for both engines falls back to equal weighting.
    """
    def w(r):
        r = float(r)
        return 1.0 / (1.0 + r) if np.isfinite(r) else 0.0
    out = {"ad": w(ad_rmsd), "gn": w(gn_rmsd)}
    return {"ad": 1.0, "gn": 1.0} if not sum(out.values()) else out


def build_variants(df: pd.DataFrame, weights: Mapping[str, float] | None = None) -> pd.DataFrame:
    """Add every variant score column to a table that ``analysis.add_scores`` has already processed."""
    out = df.copy()
    for label, col in EXTRA_SCORES.items():
        if col in out.columns:
            out[_slug(label)] = impute_worst(out[col])
    engines = {k: out[c].to_numpy(float) for k, c in ENGINE_COLUMNS.items() if c in out.columns}
    if len(engines) < 2:
        raise ValueError(f"need both engine columns {sorted(ENGINE_COLUMNS.values())}")
    out["fusion_rank_mean"] = rank_mean(engines)
    out["fusion_rrf"] = reciprocal_rank(engines)
    out["fusion_weighted_z"] = weighted_zscore(engines, weights or {"ad": 1.0, "gn": 1.0})
    if "n_heavy" in out.columns:
        for label, (col, power) in SIZE_NORMALISED.items():
            if col in out.columns:
                out[_slug(label)] = size_normalise(out[col].to_numpy(float),
                                                   out["n_heavy"].to_numpy(float), power)
    return out


def _slug(label: str) -> str:
    keep = "".join(c if (c.isalnum() or c in " ^/") else "" for c in label.lower())
    return "score_" + keep.replace(" ", "_").replace("^", "p").replace("/", "_")


# Variant labels that name the same column as a registered score. The variant comparison needs
# them, because each engine alone is one of the alternatives the consensus is tested against;
# any table that also carries the registered scores must not list the column twice under two
# names, which would read as two findings with bit-identical numbers.
REGISTERED_ALIASES = {"AutoDock-GPU alone": "autodock-gpu", "gnina alone": "gnina"}


def variant_columns(df: pd.DataFrame) -> dict[str, str]:
    """Display name -> column, for every variant present, including the single engines alone."""
    out = {}
    for label in (*EXTRA_SCORES, *SIZE_NORMALISED):
        if _slug(label) in df.columns:
            out[label] = _slug(label)
    for label, col in (("AutoDock-GPU alone", "ad_lowest"), ("gnina alone", "gn_cnn_vs"),
                       ("fusion: rank mean", "fusion_rank_mean"),
                       ("fusion: reciprocal rank", "fusion_rrf"),
                       ("fusion: redocking-weighted z", "fusion_weighted_z"),
                       # Added by the caller, not by build_variants, and only where a full interaction
                       # pass exists for the set (PROTOCOL.md Amendment 3, item D).
                       ("interaction recall", "score_interaction_recall")):
        if col in df.columns:
            out[label] = col
    return out


def compare(df: pd.DataFrame, negative_set: str, primary: str = "consensus", n_boot: int = 1000,
            seed: int = 0) -> pd.DataFrame:
    """Every variant against the registered primary score, paired on the same compounds.

    Each row carries two intervals for the difference. The pointwise one describes that single
    comparison; the Bonferroni-adjusted one targets approximately 95% simultaneous coverage over
    the whole family of variants tried against the primary at that metric and negative set.
    Several variants are compared at once, so the pointwise interval alone would declare a winner
    at the nominal rate by chance; promotion reads the adjusted interval.
    """
    if primary not in df.columns:
        raise ValueError(f"primary score column {primary!r} is missing")
    y = df["label"].to_numpy(int)
    base = df[primary].to_numpy(float)
    variants = variant_columns(df)
    # The family is the variants compared within one metric and one negative set, matching how
    # the comparison is read. Adjusting across metrics and sets as well would be more severe
    # than the claim being made.
    m = max(1, len(variants))
    resolved = M.draws_resolve_interval(n_boot, n_comparisons=m)
    rows = []
    for label, col in variants.items():
        s = df[col].to_numpy(float)
        for mname, fn in COMPARISON_METRICS.items():
            point, lo, hi = M.bootstrap_ci(fn, y, s, n_boot=n_boot, seed=seed)
            d, draws = M.paired_bootstrap_draws(fn, y, s, base, n_boot=n_boot, seed=seed)
            d_lo, d_hi = M.interval_from_draws(draws)
            a_lo, a_hi = M.interval_from_draws(draws, n_comparisons=m)
            tail = float((draws <= 0).mean()) if d > 0 else float((draws >= 0).mean())
            rows.append({"score": label, "negative_set": negative_set, "metric": mname,
                         "value": point, "ci_low": lo, "ci_high": hi,
                         "vs_primary": d, "diff_ci_low": d_lo, "diff_ci_high": d_hi,
                         "p": float(min(1.0, 2 * tail)),
                         "n_comparisons": m, "diff_adj_low": a_lo, "diff_adj_high": a_hi,
                         "adj_resolved": resolved,
                         "better": bool(d_lo > 0), "worse": bool(d_hi < 0),
                         "better_adjusted": bool(a_lo > 0), "worse_adjusted": bool(a_hi < 0)})
    return pd.DataFrame(rows)


def promoted(table: pd.DataFrame, sets=("matched", "decoy"), metric: str = "roc_auc") -> pd.DataFrame:
    """Variants whose paired interval excludes zero on *every* listed negative set.

    Amendment 3's reporting rule: an improvement seen only against random inactives does not count, because
    that set is the one the paper shows to be separable without docking at all.

    The multiplicity-adjusted interval is used where ``compare`` recorded one, since several variants are
    tried against the same primary; an older table carrying only the pointwise flag still works.
    """
    t = table[(table["metric"] == metric) & table["negative_set"].isin(sets)]
    if not len(t):
        return t.iloc[0:0]
    column = "better_adjusted" if "better_adjusted" in t.columns else "better"
    g = t.groupby("score")[column].agg(["all", "count"])
    keep = g[(g["all"]) & (g["count"] == len(set(sets)))].index
    return t[t["score"].isin(keep)].sort_values(["score", "negative_set"])
