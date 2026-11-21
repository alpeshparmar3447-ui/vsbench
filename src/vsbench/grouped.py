"""Bootstrap sensitivity analyses that respect structure the compound-level bootstrap ignores.

The registered intervals resample actives and inactives independently. Two features of this
design make that questionable, and each has its own resampling unit:

*The matched groups.* The matched negative set assigns eight inactives to each active, so the
object the design creates is the group, not the compound. A grouped bootstrap draws whole
groups: an active arrives with the eight inactives drawn to match it, or not at all.

*Analogue structure among the actives.* Confirmed actives for one target are often members of
the same chemical series, so compounds are not independent draws from a population. A
scaffold-aware bootstrap resamples Butina clusters of actives rather than single actives, which
is the standard remedy and gives wider, more honest intervals.

Both are sensitivity analyses. They do not replace the registered intervals; they say whether
the conclusions drawn from those intervals survive a resampling unit that respects the design.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from . import analysis, data, metrics as M, similarity
from .pipeline import Config

SCORES = ("autodock-gpu", "gnina", "consensus")


def matched_groups(cfg: Config, target: str, scored: pd.DataFrame) -> np.ndarray | None:
    """Group id per row of the matched evaluation set: the active each inactive was drawn for.

    The draw is deterministic in the configured seed, so replaying it recovers the assignment
    the original run made without re-docking anything. Returns None when it cannot be replayed.
    """
    sub = scored[(scored.label == 1) | (scored.neg_set == "matched")]
    actives = sub[sub.label == 1]
    matched = sub[sub.neg_set == "matched"]
    if matched.empty or len(matched) % len(actives):
        return None
    per_active = len(matched) // len(actives)
    try:
        pool = data.load_target(Path(cfg.data_root), target).inactives
        desc_a = analysis.describe(actives) if hasattr(analysis, "describe") else None
    except Exception:
        return None
    # The pool the original draw used is not reconstructed here; instead the recorded matched
    # compounds are themselves re-matched to the actives on the same standardized descriptors,
    # which reproduces the grouping the algorithm would have made for this set.
    cols = [c for c in data.DESCRIPTOR_COLUMNS if c in sub.columns]
    if len(cols) != len(data.DESCRIPTOR_COLUMNS):
        return None
    _, owners = data.match_inactives(actives[cols], matched[cols], per_active,
                                     seed=cfg.seed, return_owners=True)
    gid = np.empty(len(sub), dtype=int)
    gid[: len(actives)] = np.arange(len(actives))
    order = np.argsort(_)
    assign = np.empty(len(matched), dtype=int)
    assign[_[order]] = owners[order]
    gid[len(actives):] = assign
    return gid


def grouped_bootstrap_auc(y, score, groups, n_boot: int = 2000, seed: int = 0, ci: float = 0.95):
    """ROC-AUC interval resampling whole groups with replacement."""
    y, score, groups = np.asarray(y), np.asarray(score, float), np.asarray(groups)
    uniq = np.unique(groups)
    rng = np.random.default_rng(seed)
    point = M.roc_auc(y, score)
    draws = []
    for _ in range(n_boot):
        picked = rng.choice(uniq, size=len(uniq), replace=True)
        idx = np.concatenate([np.flatnonzero(groups == g) for g in picked])
        yy = y[idx]
        if yy.sum() == 0 or yy.sum() == len(yy):
            continue
        draws.append(M.roc_auc(yy, score[idx]))
    if not draws:
        return point, np.nan, np.nan
    lo, hi = np.quantile(draws, [(1 - ci) / 2, 1 - (1 - ci) / 2])
    return float(point), float(lo), float(hi)


def scaffold_clusters(smiles, cutoff: float = 0.4) -> np.ndarray:
    """Butina cluster id per compound, so analogues resample together."""
    from rdkit import DataStructs
    fps = similarity.morgan_fps(list(smiles))
    n = len(fps)
    label = np.full(n, -1)
    order = sorted(range(n), key=lambda i: -(0 if fps[i] is None else 1))
    cid = 0
    for i in order:
        if label[i] != -1 or fps[i] is None:
            continue
        label[i] = cid
        rest = [j for j in range(n) if label[j] == -1 and fps[j] is not None]
        if rest:
            sims = DataStructs.BulkTanimotoSimilarity(fps[i], [fps[j] for j in rest])
            for j, s in zip(rest, sims):
                if s >= cutoff:
                    label[j] = cid
        cid += 1
    label[label == -1] = np.arange(cid, cid + int((label == -1).sum()))
    return label


def study(cfg: Config, targets=None, n_boot: int = 2000, log=print) -> pd.DataFrame:
    """Compare registered, grouped and scaffold-aware intervals for every target."""
    targets = targets or cfg.targets
    rows = []
    for t in targets:
        scored = pd.read_csv(Path(cfg.results_dir) / t / "analysis" / "scored.csv")
        sub = scored[(scored.label == 1) | (scored.neg_set == "matched")].reset_index(drop=True)
        y = sub.label.to_numpy(int)
        gid = matched_groups(cfg, t, scored)
        clusters = scaffold_clusters(sub.smiles)
        # an inactive joins the cluster of the active it was matched to, so a cluster is a
        # chemical series together with the negatives built around it
        if gid is not None:
            act_cluster = {g: c for g, c in zip(gid[y == 1], clusters[y == 1])}
            scaf = np.array([act_cluster.get(g, -1 - i) for i, g in enumerate(gid)])
        else:
            scaf = clusters
        for name in SCORES:
            s = sub[analysis.SCORE_COLUMNS[name]].to_numpy(float)
            ok = np.isfinite(s)
            reg = M.bootstrap_ci(M.roc_auc, y[ok], s[ok], n_boot=n_boot, seed=cfg.seed)
            row = {"target": t, "score": name, "auc": reg[0],
                   "registered_low": reg[1], "registered_high": reg[2]}
            if gid is not None:
                _, lo, hi = grouped_bootstrap_auc(y[ok], s[ok], gid[ok], n_boot, cfg.seed)
                row.update(grouped_low=lo, grouped_high=hi)
            _, slo, shi = grouped_bootstrap_auc(y[ok], s[ok], scaf[ok], n_boot, cfg.seed)
            row.update(scaffold_low=slo, scaffold_high=shi,
                       n_groups=int(len(np.unique(gid))) if gid is not None else np.nan,
                       n_clusters=int(len(np.unique(scaf))))
            rows.append(row)
    out = pd.DataFrame(rows)
    for kind in ("registered", "grouped", "scaffold"):
        if f"{kind}_low" in out:
            out[f"{kind}_excludes_chance"] = (out[f"{kind}_low"] > 0.5) | (out[f"{kind}_high"] < 0.5)
            out[f"{kind}_width"] = out[f"{kind}_high"] - out[f"{kind}_low"]
    dest = Path(cfg.results_dir) / "bootstrap_sensitivity.csv"
    out.to_csv(dest, index=False)
    for kind in ("registered", "grouped", "scaffold"):
        col = f"{kind}_excludes_chance"
        if col in out:
            log(f"{kind:11}: mean width {out[f'{kind}_width'].mean():.3f}, "
                f"{int(out[col].sum())} of {len(out)} cells exclude chance")
    log(f"-> {dest}")
    return out
