"""Integrated lead prioritization (plan section 5.11).

Four strategies are compared on a pool fixed before any label is looked at (the top K compounds by
consensus docking score), by the fraction of known actives among the top N:
  docking        rank by consensus docking score
  docking+ADMET  drop EXCLUDE and HIGH_RISK compounds, then rank by docking
  docking+IFP    rank by the mean z-score of docking and co-crystal interaction recall
  pareto         non-dominated sorting on (docking up, interaction recall up, ADMET risk down)
A strategy that filters out actives is penalised automatically, because the metric counts actives.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import metrics as M


def zscore(v) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    sd = np.nanstd(v)
    return (v - np.nanmean(v)) / (sd if sd > 0 else 1.0)


def pareto_fronts(maximize: np.ndarray, minimize: np.ndarray | None = None) -> np.ndarray:
    """Front index (0 = non-dominated) of every row; columns of ``maximize`` are better when larger,
    columns of ``minimize`` better when smaller."""
    cols = [np.asarray(maximize, dtype=float)]
    if minimize is not None:
        cols.append(-np.asarray(minimize, dtype=float))
    X = np.column_stack(cols) if len(cols) > 1 else cols[0]
    X = X.reshape(len(X), -1)
    n = len(X)
    front = np.full(n, -1, dtype=int)
    remaining = set(range(n))
    k = 0
    while remaining:
        idx = sorted(remaining)
        nd = []
        for i in idx:
            dominated = any(np.all(X[j] >= X[i]) and np.any(X[j] > X[i]) for j in idx if j != i)
            if not dominated:
                nd.append(i)
        for i in nd:
            front[i] = k
        remaining -= set(nd)
        k += 1
    return front


def select_pool(df: pd.DataFrame, k: int, score: str = "consensus") -> pd.DataFrame:
    """The top-k compounds by docking score: chosen without reference to the labels."""
    return df.sort_values(score, ascending=False, kind="mergesort").head(k).reset_index(drop=True)


def strategies(pool: pd.DataFrame) -> dict[str, np.ndarray]:
    """Ordering (indices into ``pool``, best first) for each strategy."""
    dock = pool["consensus"].to_numpy(float)
    order = lambda s: np.argsort(-s, kind="mergesort")
    out = {"docking": order(dock)}
    keep = pool["tier"].isin(["PASS", "MODERATE"]).to_numpy() if "tier" in pool else np.ones(len(pool), bool)
    out["docking+ADMET"] = np.array([i for i in order(dock) if keep[i]], dtype=int)
    if "ifp_recall" in pool and pool["ifp_recall"].notna().all():
        out["docking+IFP"] = order(zscore(dock) + zscore(pool["ifp_recall"]))
        risk = pool["risk"].to_numpy(float) if "risk" in pool else np.zeros(len(pool))
        fronts = pareto_fronts(np.column_stack([dock, pool["ifp_recall"].to_numpy(float)]), risk[:, None])
        out["pareto"] = np.lexsort((-dock, fronts))
    return out


def hit_rate_table(pool: pd.DataFrame, n_list=(10, 20, 50), n_boot: int = 1000, seed: int = 0) -> pd.DataFrame:
    """Actives among the top N of each strategy, with a bootstrap interval over pool compounds."""
    y = pool["label"].to_numpy(int)
    base = float(y.mean())
    rng = np.random.default_rng(seed)
    rows = []
    for name, order in strategies(pool).items():
        for n in n_list:
            n_eff = min(n, len(order))
            if n_eff == 0:
                continue
            top = order[:n_eff]
            hits = int(y[top].sum())
            boots = []
            for _ in range(n_boot):
                # resample the pool; the strategy is re-applied to the compounds that were drawn
                present = np.zeros(len(pool), dtype=bool)
                present[rng.integers(0, len(pool), len(pool))] = True
                bo = order[present[order]][:n_eff]
                boots.append(float(y[bo].mean()) if len(bo) else np.nan)
            lo, hi = np.nanquantile(boots, [0.025, 0.975])
            rows.append({"strategy": name, "top_n": n, "n_evaluated": n_eff, "actives_in_top_n": hits,
                         "hit_rate": hits / n_eff, "ci_low": float(lo), "ci_high": float(hi),
                         "pool_active_fraction": base,
                         "enrichment_in_pool": (hits / n_eff) / base if base > 0 else np.nan})
    return pd.DataFrame(rows)


def filter_effect(df: pd.DataFrame, keep: pd.Series, name: str) -> dict:
    """What a filter does to the whole docked list: who it removes, and the retrieval metrics afterwards."""
    keep = keep.to_numpy(bool)
    y = df["label"].to_numpy(int)
    s = df["consensus"].to_numpy(float)
    out = {"filter": name, "active_retained": float(keep[y == 1].mean()), "inactive_retained": float(keep[y == 0].mean()),
           "actives_lost": int((~keep & (y == 1)).sum())}
    if keep.sum() and y[keep].sum() and (1 - y[keep]).sum():
        out.update({"roc_auc_after": M.roc_auc(y[keep], s[keep]), "ef5_after": M.enrichment_factor(y[keep], s[keep], 0.05)})
    else:
        out.update({"roc_auc_after": np.nan, "ef5_after": np.nan})
    out["roc_auc_before"] = M.roc_auc(y, s)
    out["ef5_before"] = M.enrichment_factor(y, s, 0.05)
    return out


def filter_table(df: pd.DataFrame) -> pd.DataFrame:
    rules = {
        "not EXCLUDE (no PAINS, Ames <= 0.5)": df["tier"] != "EXCLUDE",
        "not HIGH_RISK/EXCLUDE (also hERG <= 0.5)": ~df["tier"].isin(["EXCLUDE", "HIGH_RISK"]),
        "Veber pass": df["veber_pass"].astype(bool),
        "no PAINS alert": df["n_pains"] == 0,
    }
    return pd.DataFrame([filter_effect(df, keep, n) for n, keep in rules.items()])


def run(sets: dict, k: int, top_n=(10, 20, 50), n_boot: int = 1000, seed: int = 0) -> dict:
    """Prioritization tables for every negative set. ``sets`` maps a set name to its evaluation table
    (actives plus that negative set) carrying consensus, label and, optionally, the ADMET columns
    (tier, risk, n_pains, veber_pass) and ifp_recall."""
    hits, filters = [], []
    for name, sub in sets.items():
        pool = select_pool(sub, k)
        h = hit_rate_table(pool, top_n, n_boot, seed)
        h.insert(0, "negative_set", name)
        h["pool_size"] = len(pool)
        hits.append(h)
        if {"tier", "veber_pass", "n_pains"} <= set(sub.columns):
            f = filter_table(sub)
            f.insert(0, "negative_set", name)
            filters.append(f)
    out = {"prioritization": pd.concat(hits, ignore_index=True)}
    if filters:
        out["admet_filters"] = pd.concat(filters, ignore_index=True)
    return out
