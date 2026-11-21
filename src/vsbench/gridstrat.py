"""Docking-box strategies (plan section 4.7).

  ligand  box centred on the co-crystal ligand: the oracle choice used by the screening campaign
  pocket  box (fixed edge) centred on the top-ranked detected pocket, using no ligand information
  blind   one box over the whole protein

Each strategy redocks the deposited ligand; the RMSD is measured in place against the deposited pose, so a
strategy that finds the wrong site is penalised rather than aligned away.
"""
from __future__ import annotations

import math
import time

import numpy as np

from .box import MAX_GRID_POINTS, Box
from . import redock as R

#: Box-placement strategies. ``pocket`` is seeded from fpocket's top-ranked cavity and
#: ``p2rank`` from P2Rank's, so the pair separates how much of a pocket-seeded box's
#: performance comes from the detector rather than from the idea of seeding a box at all.
STRATEGIES = ("ligand", "pocket", "p2rank", "blind")


def pocket_box(center, edge: float = 24.0, spacing: float = 0.375) -> Box:
    return Box(center=tuple(float(c) for c in center), size=(edge, edge, edge), spacing=spacing)


def blind_box(protein_xyz, margin: float = 2.0, spacing: float = 0.375) -> Box:
    """Box over the protein's bounding box plus a margin. If that needs more than the 256 grid points per
    axis AutoGrid allows, the spacing is coarsened just enough, and the returned box records it."""
    xyz = np.asarray(protein_xyz, dtype=float)
    lo, hi = xyz.min(axis=0) - margin, xyz.max(axis=0) + margin
    size = tuple(float(v) for v in (hi - lo))
    spacing = max(spacing, math.ceil(max(size) / (MAX_GRID_POINTS - 2) * 1000) / 1000)
    return Box(center=tuple(float(v) for v in (lo + hi) / 2), size=size, spacing=spacing)


def top_pocket_center(pockets) -> tuple | None:
    """Centre of the highest-ranked pocket that has one (fpocket order: rank 1 first)."""
    for p in sorted(pockets, key=lambda p: p.rank):
        c = p.extra.get("center")
        if c is not None:
            return c
    return None


def strategy_row(strategy: str, res: R.RedockResult, wall_s: float, box: Box | None) -> dict:
    return {"template": res.template, "strategy": strategy, "ok": res.ok, "error": res.error,
            "box_edge": res.box_edge, "spacing": box.spacing if box else np.nan,
            "ad_lowest_rmsd": res.ad_lowest_rmsd, "gn_top_rmsd": res.gn_top_rmsd,
            "ad_success": bool(res.ok and res.ad_lowest_rmsd <= R.SUCCESS_RMSD),
            "gn_success": bool(res.ok and res.gn_top_rmsd <= R.SUCCESS_RMSD),
            "wall_s": wall_s}


def run_strategy(strategy: str, protein_mol2, ligand_mol2, template: str, workdir, box: Box | None, **kw) -> dict:
    t0 = time.time()
    res = R.redock_template(protein_mol2, ligand_mol2, template, workdir, box_override=box, **kw)
    return strategy_row(strategy, res, time.time() - t0, box)


def summarize(rows) -> "pd.DataFrame":
    """Success rate per strategy and engine with exact (Wilson) intervals."""
    import pandas as pd
    df = pd.DataFrame(rows)
    out = []
    # Report in the declared order, so a strategy added later still appears beside the one it
    # is meant to be compared with rather than tacked on at the end.
    order = [s for s in STRATEGIES if s in set(df["strategy"])]
    df = df.assign(_o=df["strategy"].map({s: i for i, s in enumerate(order)}))
    for strat, g in sorted(df.groupby("strategy", sort=False), key=lambda kv: kv[1]["_o"].iloc[0]):
        n = len(g)
        for eng in ("ad", "gn"):
            k = int(g[f"{eng}_success"].sum())
            lo, hi = wilson(k, n)
            out.append({"strategy": strat, "engine": {"ad": "autodock-gpu", "gn": "gnina"}[eng], "n": n,
                        "successes": k, "success_rate": k / n if n else np.nan, "ci_low": lo, "ci_high": hi,
                        "median_rmsd": float(np.nanmedian(g[f"{eng}_lowest_rmsd" if eng == "ad" else "gn_top_rmsd"])),
                        "median_wall_s": float(np.nanmedian(g["wall_s"]))})
    return pd.DataFrame(out)


def paired_strategies(rows, engines=("ad", "gn"), scope: str = "pairwise") -> "pd.DataFrame":
    """Every pair of box strategies compared on the templates where both ran, exact McNemar.

    Success rates in :func:`summarize` are unpaired, so two of them can differ because the
    strategies were applied to different template sets. This pairs them per template and tests
    the discordant counts, which is what a claim that one box beats another has to rest on.

    Preparation fails unevenly across strategies, so only templates that actually ran are
    compared. ``scope`` decides how that is read. ``"pairwise"`` compares each pair on the
    templates those two share, which keeps every comparison as powered as the data allow and,
    importantly, leaves a comparison unchanged when a further strategy is added to the study.
    ``"listwise"`` restricts every comparison to the templates on which *all* strategies ran,
    which puts them on one common set at the cost of making each result depend on which other
    strategies happen to be included.
    """
    import pandas as pd
    from scipy.stats import binomtest

    df = pd.DataFrame(rows)
    df = df[df.get("ok", True).astype(bool)] if "ok" in df else df
    strategies = [s for s in STRATEGIES if s in set(df["strategy"])]

    if scope == "listwise":
        counts = df.groupby("template").strategy.nunique()
        df = df[df.template.isin(set(counts[counts == len(strategies)].index))]
    elif scope != "pairwise":
        raise ValueError(f"scope must be 'pairwise' or 'listwise', not {scope!r}")

    out = []
    for engine in engines:
        wide = df.pivot_table(index="template", columns="strategy",
                              values=f"{engine}_success", aggfunc="first")
        for i, a in enumerate(strategies):
            for b in strategies[i + 1:]:
                if a not in wide or b not in wide:
                    continue
                both = wide[[a, b]].dropna().astype(bool)
                n_b = int((both[a] & ~both[b]).sum())     # a succeeds, b fails
                n_c = int((~both[a] & both[b]).sum())     # b succeeds, a fails
                n_disc = n_b + n_c
                p = float(binomtest(min(n_b, n_c), n_disc, 0.5).pvalue) if n_disc else 1.0
                out.append({
                    "engine": {"ad": "autodock-gpu", "gn": "gnina"}[engine],
                    "strategy_a": a, "strategy_b": b, "scope": scope,
                    "n_paired": len(both), "a_only": n_b, "b_only": n_c,
                    "discordant": n_disc, "p_mcnemar_exact": p,
                    "favours": ((a if n_b > n_c else b) if p < 0.05 and n_b != n_c
                                else "not resolved"),
                })
    return pd.DataFrame(out)


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    # Clamp: at k = 0 or k = n the arithmetic lands a hair outside [0, 1], and a reported rate
    # interval of -0.00 is wrong on its face.
    return (max(0.0, c - h), min(1.0, c + h))
