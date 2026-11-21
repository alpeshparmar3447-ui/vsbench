"""Separate the retrospective ligand-similarity baseline from a deployable one.

The registered similarity baseline scores every compound by its maximum ECFP4 Tanimoto to the
*evaluation set's own actives*, with actives scored leave-one-out. That uses the labels the
benchmark is testing, so it is an oracle: an upper bound on what ligand similarity could extract
from this data, not a method a prospective screen could run.

This module computes the deployable counterpart. The reference set is the target's actives that
were **not** sampled into the evaluation set, so no evaluation label is used and the resulting
score is one anybody could compute before knowing the answer. Both baselines are then compared
against the docking scores on the same compounds, which is what decides whether the paper's
central claim rests on an oracle or survives without one.

A target whose actives were all sampled into the evaluation set has no held-out reference and is
reported as such rather than silently dropped.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from . import analysis, data, metrics as M, similarity
from .pipeline import Config

REGISTERED_SCORES = ("autodock-gpu", "gnina", "consensus")


def held_out_actives(cfg: Config, target: str, evaluated_ids) -> pd.DataFrame:
    """The target's actives that the evaluation set did not take."""
    tg = data.load_target(Path(cfg.data_root), target)
    used = {str(i) for i in evaluated_ids}
    acts = tg.actives.copy()
    acts["id"] = acts["id"].astype(str)
    return acts[~acts["id"].isin(used)]


def baselines_for_target(cfg: Config, target: str, n_boot: int = 5000) -> pd.DataFrame:
    """Oracle and deployable similarity baselines, and docking against each, for one target."""
    scored = pd.read_csv(Path(cfg.results_dir) / target / "analysis" / "scored.csv")
    act_ids = scored.loc[scored.label == 1, "name"]
    held = held_out_actives(cfg, target, act_ids)
    n_held = len(held)

    rows = []
    for neg_set, sub in analysis.evaluation_sets(scored).items():
        y = sub["label"].to_numpy(int)
        if y.sum() == 0 or (len(y) - y.sum()) == 0:
            continue
        oracle = sub["nn_sim"].to_numpy(float)
        deployable = (similarity.reference_similarity(list(sub["smiles"]), list(held["smiles"]))
                      if n_held else np.full(len(sub), np.nan))

        for label, score in (("similarity (oracle, evaluation-set actives)", oracle),
                             ("similarity (deployable, held-out actives)", deployable)):
            ok = np.isfinite(score)
            if ok.sum() < 2 or len(set(y[ok])) < 2:
                rows.append({"target": target, "negative_set": neg_set, "baseline": label,
                             "n_reference_actives": n_held, "roc_auc": np.nan,
                             "note": "no held-out actives"})
                continue
            point, lo, hi = M.bootstrap_ci(M.roc_auc, y[ok], score[ok], n_boot=1000, seed=cfg.seed)
            row = {"target": target, "negative_set": neg_set, "baseline": label,
                   "n_reference_actives": n_held if "deployable" in label else int(y.sum()) - 1,
                   "roc_auc": point, "ci_low": lo, "ci_high": hi,
                   "excludes_chance": bool(lo > 0.5 or hi < 0.5)}
            # docking against this baseline, same compounds, same adjustment rule as elsewhere
            for name in REGISTERED_SCORES:
                col = analysis.SCORE_COLUMNS[name]
                s = sub[col].to_numpy(float)
                both = ok & np.isfinite(s)
                d, draws = M.paired_bootstrap_draws(M.roc_auc, y[both], s[both], score[both],
                                                    n_boot=n_boot, seed=cfg.seed)
                alo, ahi = M.interval_from_draws(draws, n_comparisons=len(REGISTERED_SCORES))
                row[f"{name}_auc"] = M.roc_auc(y[both], s[both])
                row[f"{name}_minus_baseline"] = d
                row[f"{name}_beats"] = bool(alo > 0)
                row[f"{name}_below"] = bool(ahi < 0)
            rows.append(row)
    return pd.DataFrame(rows)


def study(cfg: Config, targets=None, log=print) -> pd.DataFrame:
    """Run both baselines for every target and write ``baseline_leakage.csv``."""
    targets = targets or cfg.targets
    out = pd.concat([baselines_for_target(cfg, t) for t in targets], ignore_index=True)
    dest = Path(cfg.results_dir) / "baseline_leakage.csv"
    out.to_csv(dest, index=False)

    for label, g in out.groupby("baseline"):
        hard = g[g.negative_set != "random"].dropna(subset=["roc_auc"])
        if hard.empty:
            continue
        beats = sum(int(hard[f"{s}_beats"].sum()) for s in REGISTERED_SCORES)
        below = sum(int(hard[f"{s}_below"].sum()) for s in REGISTERED_SCORES)
        n = len(hard) * len(REGISTERED_SCORES)
        log(f"{label}: mean AUC {hard.roc_auc.mean():.3f} over {len(hard)} target-set cells; "
            f"docking beats it in {beats}/{n}, falls below in {below}/{n}")
    log(f"-> {dest}")
    return out
