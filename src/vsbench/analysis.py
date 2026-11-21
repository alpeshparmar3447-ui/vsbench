"""Benchmark analysis: scores, bias audit, pose-selection effect, failure accounting.

Score columns follow the engine wrappers. All scores handed to the metrics are
"higher is better"; a failed docking (NaN) is imputed with that engine's worst
finite score, so a failure is counted against the compound and never dropped.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import metrics as M
from .data import DESCRIPTOR_COLUMNS, standardized_mean_difference

SCORE_COLUMNS = {
    "autodock-gpu": "ad_lowest",          # AD4 free energy of the lowest-energy pose
    "autodock-gpu (largest cluster)": "ad_largest",
    "gnina": "gn_cnn_vs",                 # CNN score x CNN affinity, gnina's virtual-screening score
    "gnina (Vina-type)": "gn_vina_hb",
    "consensus": "consensus",             # mean z-score of the two primary engine scores
}
PRIMARY = ("autodock-gpu", "gnina", "consensus")
BASELINE_COLUMNS = ["n_heavy", "mw", "clogp", "tpsa", "hba", "hbd", "rotb", "nn_sim"]


def impute_worst(v) -> np.ndarray:
    v = np.asarray(v, dtype=float).copy()
    finite = np.isfinite(v)
    if not finite.any():
        raise ValueError("no finite scores to impute from")
    v[~finite] = v[finite].min()
    return v


def add_scores(df: pd.DataFrame) -> pd.DataFrame:
    """Add higher-is-better score columns and the consensus to a screening table."""
    out = df.copy()
    out["ad_lowest"] = M.to_higher_better(out["ad_energy"])
    out["ad_largest"] = M.to_higher_better(out["ad_largest_energy"])
    out["gn_vina_hb"] = M.to_higher_better(out["gn_vina"])
    for c in ("ad_lowest", "ad_largest", "gn_cnn_vs", "gn_vina_hb"):
        out[c] = impute_worst(out[c])
    out["consensus"] = M.consensus_zscore({"ad": out["ad_lowest"].to_numpy(),
                                           "gn": out["gn_cnn_vs"].to_numpy()})
    return out


def failure_table(df: pd.DataFrame) -> pd.DataFrame:
    """Docking and preparation failures by class (counts and fractions)."""
    rows = []
    for label, name in ((1, "active"), (0, "inactive")):
        sub = df[df["label"] == label]
        rows.append({
            "class": name, "n": len(sub),
            "ligand_prep_failed": int((~sub["prep_ok"]).sum()),
            "autodock_failed": int((sub["prep_ok"] & ~sub["ad_ok"]).sum()),
            "gnina_failed": int((sub["prep_ok"] & ~sub["gn_ok"]).sum()),
            "gnina_clash_gt0": int((sub["gn_vina"] > 0).sum()),
            "autodock_positive": int((sub["ad_energy"] > 0).sum()),
            "frozen_flexible": int(sub["frozen_flexible"].fillna(False).astype(bool).sum()) if "frozen_flexible" in sub else 0,
            "charge_flag": int(sub["charge_flag"].fillna(False).astype(bool).sum()) if "charge_flag" in sub else 0,
            "chelator_needs_manual_ionisation": (int((sub["chelators_manual"].fillna("") != "").sum())
                                                 if "chelators_manual" in sub else 0),
            "chelator_deprotonated": int(sub["deprotonated"].fillna(False).astype(bool).sum())
            if "deprotonated" in sub else 0,
        })
    t = pd.DataFrame(rows)
    for c in ("ligand_prep_failed", "autodock_failed", "gnina_failed", "gnina_clash_gt0", "autodock_positive",
              "chelator_needs_manual_ionisation", "chelator_deprotonated", "frozen_flexible", "charge_flag"):
        t[c + "_frac"] = t[c] / t["n"]
    return t


def evaluation_sets(df: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Actives plus each negative set, keyed by negative-set name.

    The screening table carries a ``neg_set`` column: 'active', 'random', 'matched' and, when the
    decoy arm was run, 'decoy'.
    """
    act = df[df["neg_set"] == "active"]
    present = [n for n in ("random", "matched", "decoy") if (df["neg_set"] == n).any()]
    return {name: pd.concat([act, df[df["neg_set"] == name]]).reset_index(drop=True) for name in present}


def negative_set_difference(set_a: pd.DataFrame, set_b: pd.DataFrame, score: str, n_boot: int = 1000,
                            seed: int = 0) -> dict:
    """AUC(set_a) - AUC(set_b) for one score when the two sets share their actives but not their negatives.

    Actives are resampled together, so the shared actives are treated as paired; each set's negatives
    are resampled independently, because they are different compounds.
    """
    act_ids = sorted(set(set_a.loc[set_a["label"] == 1, "name"]) & set(set_b.loc[set_b["label"] == 1, "name"]))
    pa = set_a[set_a["label"] == 1].set_index("name").loc[act_ids, score].to_numpy(float)
    pb = set_b[set_b["label"] == 1].set_index("name").loc[act_ids, score].to_numpy(float)
    na = set_a.loc[set_a["label"] == 0, score].to_numpy(float)
    nb = set_b.loc[set_b["label"] == 0, score].to_numpy(float)

    def auc(pos, neg):
        return M.roc_auc(np.r_[np.ones(len(pos)), np.zeros(len(neg))], np.r_[pos, neg])

    rng = np.random.default_rng(seed)
    diffs = np.empty(n_boot)
    for b in range(n_boot):
        ai = rng.integers(0, len(act_ids), len(act_ids))
        diffs[b] = (auc(pa[ai], na[rng.integers(0, len(na), len(na))])
                    - auc(pb[ai], nb[rng.integers(0, len(nb), len(nb))]))
    point = auc(pa, na) - auc(pb, nb)
    lo, hi = np.quantile(diffs, [0.025, 0.975])
    tail = float(np.mean(diffs <= 0)) if point > 0 else float(np.mean(diffs >= 0))
    return {"delta_auc": float(point), "ci_low": float(lo), "ci_high": float(hi), "p_two_sided": min(1.0, 2 * tail)}


def balance_table(act: pd.DataFrame, neg: pd.DataFrame) -> pd.Series:
    return standardized_mean_difference(act[DESCRIPTOR_COLUMNS], neg[DESCRIPTOR_COLUMNS])


def evaluate(df: pd.DataFrame, negative_set: str, n_boot: int = 1000, seed: int = 0,
             columns: dict | None = None) -> pd.DataFrame:
    """Every metric with a bootstrap interval, for each score column."""
    columns = columns or SCORE_COLUMNS
    y = df["label"].to_numpy(int)
    rows = []
    for name, col in columns.items():
        s = df[col].to_numpy(float)
        for mname, fn in M.METRICS.items():
            point, lo, hi = M.bootstrap_ci(fn, y, s, n_boot=n_boot, seed=seed)
            rows.append({"score": name, "negative_set": negative_set, "metric": mname,
                         "value": point, "ci_low": lo, "ci_high": hi,
                         "n_active": int(y.sum()), "n_inactive": int(len(y) - y.sum())})
    return pd.DataFrame(rows)


def property_baselines(df: pd.DataFrame, negative_set: str, n_boot: int = 1000, seed: int = 0) -> pd.DataFrame:
    """ROC-AUC of single properties used as 'scores'. Values far from 0.5 mean the
    negative set is separable without any docking."""
    y = df["label"].to_numpy(int)
    rows = []
    for c in (c for c in BASELINE_COLUMNS if c in df.columns):
        s = df[c].to_numpy(float)
        if np.isnan(s).any():
            s = np.where(np.isnan(s), np.nanmin(s), s)
        point, lo, hi = M.bootstrap_ci(M.roc_auc, y, s, n_boot=n_boot, seed=seed)
        rows.append({"property": c, "negative_set": negative_set, "roc_auc": point,
                     "ci_low": lo, "ci_high": hi})
    return pd.DataFrame(rows)


def baseline_comparison(df: pd.DataFrame, negative_set: str, n_boot: int = 5000,
                        seed: int = 0) -> pd.DataFrame:
    """Paired bootstrap of each docking score against the best structure-free baseline.

    The baseline family is BASELINE_COLUMNS: seven physical descriptors plus the
    nearest-neighbour similarity to the actives. The best of the family is chosen by its
    own ROC-AUC on this negative set, then compared against each docking score on the
    same resampled compounds. Because the comparison is read against a family of m
    baselines, the difference carries a Bonferroni-widened interval over m alongside the
    pointwise one (PROTOCOL.md, Amendment 5); the adjusted interval governs both the
    claim that docking exceeds the baseline and the claim that it falls below.
    """
    y = df["label"].to_numpy(int)
    props = [c for c in BASELINE_COLUMNS if c in df.columns]
    if not props or y.sum() == 0 or (len(y) - y.sum()) == 0:
        return pd.DataFrame()
    m = len(props)

    def _clean(col: str) -> np.ndarray:
        s = df[col].to_numpy(float)
        return np.where(np.isnan(s), np.nanmin(s), s) if np.isnan(s).any() else s

    aucs = {c: M.roc_auc(y, _clean(c)) for c in props}
    best = max(aucs, key=aucs.get)
    bvals = _clean(best)

    rows = []
    for name, col in SCORE_COLUMNS.items():
        if col not in df.columns:
            continue
        s = _clean(col)
        d, draws = M.paired_bootstrap_draws(M.roc_auc, y, s, bvals, n_boot=n_boot, seed=seed)
        lo, hi = M.interval_from_draws(draws)
        alo, ahi = M.interval_from_draws(draws, n_comparisons=m)
        rows.append({"score": name, "negative_set": negative_set,
                     "best_baseline": best, "baseline_auc": aucs[best],
                     "score_auc": M.roc_auc(y, s), "delta_auc": d,
                     "ci_low": lo, "ci_high": hi,
                     "adj_ci_low": alo, "adj_ci_high": ahi,
                     "n_baselines": m,
                     "adj_resolved": M.draws_resolve_interval(n_boot, n_comparisons=m),
                     "beats_baseline": bool(alo > 0),
                     "below_baseline": bool(ahi < 0)})
    return pd.DataFrame(rows)


def pose_selection_effect(df: pd.DataFrame, n_boot: int = 1000, seed: int = 0) -> dict:
    """Does ranking by the lowest-energy pose differ from ranking by the most
    populated cluster's lowest energy? Paired bootstrap on ROC-AUC."""
    y = df["label"].to_numpy(int)
    d, lo, hi, p = M.paired_bootstrap_difference(
        M.roc_auc, y, df["ad_largest"].to_numpy(), df["ad_lowest"].to_numpy(), n_boot=n_boot, seed=seed)
    return {"delta_auc_largest_minus_lowest": d, "ci_low": lo, "ci_high": hi, "p_two_sided": p,
            "auc_lowest": M.roc_auc(y, df["ad_lowest"].to_numpy()),
            "auc_largest": M.roc_auc(y, df["ad_largest"].to_numpy())}


def seed_stability(df: pd.DataFrame, seeds, n_top: int = 20) -> dict:
    """How much the AutoDock-GPU result depends on the random seed (needs per-seed energy columns)."""
    from scipy.stats import spearmanr
    y = df["label"].to_numpy(int)
    cols = [f"ad_energy_s{s}" for s in seeds]
    E = np.column_stack([impute_worst(-df[c].to_numpy(float)) for c in cols])     # higher = better
    aucs = [M.roc_auc(y, E[:, j]) for j in range(E.shape[1])]
    rhos = [spearmanr(E[:, i], E[:, j])[0] for i in range(E.shape[1]) for j in range(i + 1, E.shape[1])]
    tops = [set(np.argsort(-E[:, j], kind="mergesort")[:n_top]) for j in range(E.shape[1])]
    shared = set.intersection(*tops) if tops else set()
    sd = np.nanstd(np.column_stack([df[c].to_numpy(float) for c in cols]), axis=1, ddof=1)
    return {"n_seeds": len(seeds), "auc_min": float(min(aucs)), "auc_max": float(max(aucs)),
            "auc_range": float(max(aucs) - min(aucs)), "mean_spearman_between_seeds": float(np.mean(rhos)) if rhos else float("nan"),
            "stable_top_n_fraction": len(shared) / n_top, "mean_per_compound_sd": float(np.nanmean(sd)),
            "max_per_compound_sd": float(np.nanmax(sd))}


def gnina_seed_study(df: pd.DataFrame, seeds) -> dict:
    """Seed variance of gnina on the subset that was re-docked (rows with every per-seed column)."""
    from scipy.stats import spearmanr
    base = "gn_cnn_vs"
    cols = [base] + [f"{base}_s{s}" for s in seeds]
    sub = df.dropna(subset=cols)
    if sub["label"].nunique() < 2:
        return {}
    y = sub["label"].to_numpy(int)
    aucs = [M.roc_auc(y, sub[c].to_numpy(float)) for c in cols]
    rhos = [spearmanr(sub[cols[i]], sub[cols[j]])[0] for i in range(len(cols)) for j in range(i + 1, len(cols))]
    return {"n": len(sub), "n_active": int(y.sum()), "auc_min": float(min(aucs)), "auc_max": float(max(aucs)),
            "auc_range": float(max(aucs) - min(aucs)), "mean_spearman_between_seeds": float(np.mean(rhos)),
            "mean_abs_cnn_vs_difference": float(np.mean(np.abs(sub[cols[0]] - sub[cols[1]])))}


def engine_difference(df: pd.DataFrame, n_boot: int = 1000, seed: int = 0) -> dict:
    """Paired comparison of the two primary engines (gnina minus AutoDock-GPU)."""
    y = df["label"].to_numpy(int)
    d, lo, hi, p = M.paired_bootstrap_difference(
        M.roc_auc, y, df["gn_cnn_vs"].to_numpy(), df["ad_lowest"].to_numpy(), n_boot=n_boot, seed=seed)
    return {"delta_auc_gnina_minus_autodock": d, "ci_low": lo, "ci_high": hi, "p_two_sided": p}
