"""Publication figures from the benchmark CSVs (matplotlib, Agg backend)."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from . import metrics as M  # noqa: E402

ENGINE_COLORS = {"autodock-gpu": "#0072B2", "gnina": "#D55E00", "consensus": "#009E73"}
ENGINES = tuple(ENGINE_COLORS)


VECTOR_FORMATS = ("pdf", "svg")


def _finish(fig, outfile: Path, dpi: int = 300) -> Path:
    """Write the raster figure and a vector copy of it beside it.

    Journals want vector artwork; the raster path is what the report links. Both are
    written from the same figure object, so they cannot diverge.
    """
    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outfile, dpi=dpi, bbox_inches="tight")
    for ext in VECTOR_FORMATS:
        fig.savefig(outfile.with_suffix(f".{ext}"), bbox_inches="tight")
    plt.close(fig)
    return outfile


def plot_auc_summary(metrics_all: pd.DataFrame, baselines_all: pd.DataFrame | None, outfile: Path) -> Path:
    """ROC-AUC per target: engines with bootstrap intervals, on random and matched negatives.

    ``baselines_all`` (columns target, negative_set, property, roc_auc) adds the
    heavy-atom-count baseline as a black marker, so any docking AUC that sits no
    higher than that marker is visibly no better than counting atoms.
    """
    auc = metrics_all[metrics_all["metric"] == "roc_auc"]
    targets = list(dict.fromkeys(auc["target"]))
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    for ax, neg in zip(axes, ("random", "matched")):
        for j, eng in enumerate(ENGINES):
            sub = auc[(auc["negative_set"] == neg) & (auc["score"] == eng)].set_index("target").reindex(targets)
            x = np.arange(len(targets)) + (j - 1) * 0.22
            ax.errorbar(x, sub["value"], yerr=[sub["value"] - sub["ci_low"], sub["ci_high"] - sub["value"]],
                        fmt="o", ms=5, capsize=2, color=ENGINE_COLORS[eng], label=eng)
        if baselines_all is not None:
            b = baselines_all[(baselines_all["negative_set"] == neg) & (baselines_all["property"] == "n_heavy")]
            b = b.set_index("target").reindex(targets)
            ax.scatter(np.arange(len(targets)), b["roc_auc"], marker="_", s=300, color="black",
                       label="heavy-atom count alone", zorder=3)
        ax.axhline(0.5, ls="--", lw=0.8, color="grey")
        ax.set_xticks(range(len(targets)), targets, rotation=30, ha="right")
        ax.set_title(f"{neg} inactives")
        ax.set_ylim(0.2, 1.0)
    axes[0].set_ylabel("ROC-AUC")
    axes[1].legend(loc="upper right", fontsize=8, frameon=False)
    fig.tight_layout()
    return _finish(fig, outfile)


def plot_roc(scored: pd.DataFrame, negative_set: str, outfile: Path, title: str = "") -> Path:
    """ROC curves of the three primary scores on one target and negative set."""
    act = scored[scored["neg_set"] == "active"]
    sub = pd.concat([act, scored[scored["neg_set"] == negative_set]])
    y = sub["label"].to_numpy(int)
    cols = {"autodock-gpu": "ad_lowest", "gnina": "gn_cnn_vs", "consensus": "consensus"}
    fig, ax = plt.subplots(figsize=(3.8, 3.8))
    for name, col in cols.items():
        s = sub[col].to_numpy(float)
        order = np.argsort(-s, kind="mergesort")
        tp = np.cumsum(y[order]) / y.sum()
        fp = np.cumsum(1 - y[order]) / (len(y) - y.sum())
        ax.plot(np.r_[0, fp], np.r_[0, tp], color=ENGINE_COLORS[name],
                label=f"{name} ({M.roc_auc(y, s):.2f})")
    ax.plot([0, 1], [0, 1], ls="--", lw=0.8, color="grey")
    ax.set(xlabel="False-positive rate", ylabel="True-positive rate", title=title or negative_set)
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    fig.tight_layout()
    return _finish(fig, outfile)


def plot_redock(redock_all: pd.DataFrame, outfile: Path, threshold: float = 2.0) -> Path:
    """Redocking RMSD of every template: AutoDock-GPU against gnina, one panel per target."""
    targets = list(dict.fromkeys(redock_all["target"]))
    fig, axes = plt.subplots(1, len(targets), figsize=(2.6 * len(targets), 2.8), sharex=True, sharey=True,
                             squeeze=False)
    for ax, t in zip(axes[0], targets):
        d = redock_all[(redock_all["target"] == t) & redock_all["ok"]]
        ax.scatter(d["ad_lowest_rmsd"], d["gn_top_rmsd"], s=18, color="#444444")
        ax.axvline(threshold, ls="--", lw=0.8, color="grey")
        ax.axhline(threshold, ls="--", lw=0.8, color="grey")
        ax.set_title(t, fontsize=9)
        ax.set_xlim(0, 14)
        ax.set_ylim(0, 14)
    axes[0][0].set_ylabel("gnina top pose RMSD (A)")
    for ax in axes[0]:
        ax.set_xlabel("AutoDock-GPU RMSD (A)")
    fig.tight_layout()
    return _finish(fig, outfile)


# ---------------------------------------------------------------- extended-pipeline figures
SET_COLORS = {"random": "#999999", "matched": "#0072B2", "decoy": "#D55E00"}
STRATEGY_COLORS = {"docking": "#999999", "docking+ADMET": "#0072B2", "docking+IFP": "#E69F00", "pareto": "#009E73"}


def plot_negative_set_auc(metrics_all: pd.DataFrame, score: str, outfile: Path) -> Path:
    """ROC-AUC of one score per target for each negative-set design (bars with bootstrap intervals)."""
    d = metrics_all[(metrics_all["metric"] == "roc_auc") & (metrics_all["score"] == score)]
    targets = sorted(d["target"].unique())
    sets = [s for s in SET_COLORS if s in set(d["negative_set"])]
    fig, ax = plt.subplots(figsize=(max(6.0, 0.55 * len(targets) * len(sets) / 2 + 3), 4.0))
    w = 0.8 / max(len(sets), 1)
    for j, s in enumerate(sets):
        g = d[d["negative_set"] == s].set_index("target").reindex(targets)
        x = np.arange(len(targets)) + (j - (len(sets) - 1) / 2) * w
        ax.bar(x, g["value"], w, color=SET_COLORS[s], label=s,
               yerr=[g["value"] - g["ci_low"], g["ci_high"] - g["value"]], capsize=1.5, error_kw={"lw": 0.7})
    ax.axhline(0.5, color="k", lw=0.7, ls="--")
    ax.set_xticks(np.arange(len(targets)), targets, rotation=60, ha="right")
    ax.set_ylabel(f"ROC-AUC ({score})")
    ax.set_ylim(0, 1)
    ax.legend(title="negatives", frameon=False, ncol=len(sets))
    return _finish(fig, outfile)


def plot_site_recovery(sites: pd.DataFrame, outfile: Path) -> Path:
    """Fraction of templates whose top-1 / top-3 pocket recovers the ligand site, per detector."""
    from .gridstrat import wilson
    methods = list(dict.fromkeys(sites["method"]))
    fig, ax = plt.subplots(figsize=(4.2, 3.6))
    for j, col in enumerate(("top1_hit", "top3_hit")):
        for i, m in enumerate(methods):
            g = sites[sites["method"] == m][col].astype(bool)
            k, n = int(g.sum()), len(g)
            lo, hi = wilson(k, n)
            x = i + (j - 0.5) * 0.36
            ax.bar(x, k / n, 0.34, color=["#56B4E9", "#0072B2"][j], label=("top-1" if j == 0 else "top-3") if i == 0 else None,
                   yerr=[[k / n - lo], [hi - k / n]], capsize=2, error_kw={"lw": 0.8})
    ax.set_xticks(range(len(methods)), methods)
    ax.set_ylabel("templates with the ligand site recovered")
    ax.set_ylim(0, 1.05)
    ax.legend(frameon=False)
    return _finish(fig, outfile)


def plot_prioritization(prio: pd.DataFrame, outfile: Path, negative_set: str = "matched") -> Path:
    """Mean hit rate (actives among the top N) of each lead-selection strategy, across targets."""
    d = prio[prio["negative_set"] == negative_set]
    top_ns = sorted(d["top_n"].unique())
    strategies = [s for s in STRATEGY_COLORS if s in set(d["strategy"])]
    fig, ax = plt.subplots(figsize=(5.0, 3.6))
    w = 0.8 / max(len(strategies), 1)
    for j, s in enumerate(strategies):
        g = d[d["strategy"] == s].groupby("top_n")["hit_rate"].agg(["mean", "sem"]).reindex(top_ns)
        x = np.arange(len(top_ns)) + (j - (len(strategies) - 1) / 2) * w
        ax.bar(x, g["mean"], w, color=STRATEGY_COLORS[s], label=s, yerr=g["sem"].fillna(0), capsize=2, error_kw={"lw": 0.8})
    base = d.groupby("top_n")["pool_active_fraction"].mean().reindex(top_ns)
    ax.plot(np.arange(len(top_ns)), base, "k_", ms=22, mew=1.2, label="pool active fraction")
    ax.set_xticks(np.arange(len(top_ns)), [f"top {n}" for n in top_ns])
    ax.set_ylabel("fraction of selected compounds that are active")
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, fontsize=7)
    return _finish(fig, outfile)


def plot_grid_strategies(summary: pd.DataFrame, outfile: Path) -> Path:
    """Redocking success rate (RMSD <= 2 A, measured in place) for each docking-box strategy and engine."""
    # "pocket" is fpocket's top cavity; name the detector on the axis so the two pocket-seeded
    # boxes cannot be read as one generic "predicted pocket" strategy.
    labels = {"ligand": "ligand-\ncentred", "pocket": "fpocket\ntop pocket",
              "p2rank": "P2Rank\ntop pocket", "blind": "whole\nprotein"}
    engines = list(dict.fromkeys(summary["engine"]))
    strategies = list(dict.fromkeys(summary["strategy"]))
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    w = 0.8 / max(len(engines), 1)
    for j, e in enumerate(engines):
        g = summary[summary["engine"] == e].set_index("strategy").reindex(strategies)
        x = np.arange(len(strategies)) + (j - (len(engines) - 1) / 2) * w
        ax.bar(x, g["success_rate"], w, color=ENGINE_COLORS.get(e, "#777777"), label=e,
               yerr=[g["success_rate"] - g["ci_low"], g["ci_high"] - g["success_rate"]], capsize=2, error_kw={"lw": 0.8})
    ax.set_xticks(np.arange(len(strategies)), [labels.get(s, s) for s in strategies])
    ax.set_ylabel("redocking success (RMSD <= 2 A)")
    ax.set_ylim(0, 1.05)
    ax.legend(frameon=False)
    return _finish(fig, outfile)


def plot_performance(perf: pd.DataFrame, outfile: Path) -> Path:
    """Throughput (compounds per hour, log scale) of each case at its largest ligand count."""
    d = perf.sort_values("n").groupby(["case", "threads"], as_index=False).last()
    d = d.sort_values("compounds_per_hour")
    labels = [f"{c}\n{int(t)} thr" if "GPU" not in c and c != "autodock-gpu" else c for c, t in zip(d["case"], d["threads"])]
    fig, ax = plt.subplots(figsize=(6.0, 0.45 * len(d) + 1.2))
    colors = ["#D55E00" if ("GPU" in c or c == "autodock-gpu") else "#999999" for c in d["case"]]
    ax.barh(range(len(d)), d["compounds_per_hour"], color=colors)
    ax.set_xscale("log")
    ax.set_yticks(range(len(d)), labels, fontsize=7)
    ax.set_xlabel("compounds per hour")
    return _finish(fig, outfile)


def plot_size_control(macro: pd.DataFrame, gaps: pd.DataFrame | None, outfile: Path,
                      scores=("autodock-gpu", "gnina", "consensus")) -> Path:
    """Plain against size-controlled AUC per negative set, with the heavy-atom gap that explains it.

    The left panel pairs each score's ordinary AUC with its AUC among compounds of the same size class; the
    arrow between them is the part of the apparent performance that was molecular size. The right panel
    shows the median heavy-atom gap between actives and each negative set, which is what drives the arrow.
    """
    sets = [s for s in ("random", "matched", "decoy") if (macro["negative_set"] == s).any()]
    use = macro[macro["score"].isin(scores)]
    fig, axes = plt.subplots(1, 2 if gaps is not None else 1, figsize=(11 if gaps is not None else 7, 4.2),
                             gridspec_kw={"width_ratios": [2, 1]} if gaps is not None else None)
    ax = axes[0] if gaps is not None else axes
    width, x0 = 0.26, np.arange(len(sets), dtype=float)
    for i, sc in enumerate([s for s in scores if s in set(use["score"])]):
        d = use[use["score"] == sc].set_index("negative_set")
        xs = x0 + (i - 1) * width
        plain = [d.loc[s, "roc_auc"] if s in d.index else np.nan for s in sets]
        cond = [d.loc[s, "conditional_auc"] if s in d.index else np.nan for s in sets]
        c = ENGINE_COLORS.get(sc, "#555555")
        ax.bar(xs, plain, width * 0.9, color=c, alpha=0.35,
               label="as measured (bars)" if i == 0 else None)
        ax.plot(xs, cond, "o", color=c, ms=7, mec="black", mew=0.6, label=sc)
        for x, p, q in zip(xs, plain, cond):
            if np.isfinite(p) and np.isfinite(q) and abs(p - q) > 0.004:
                ax.annotate("", xy=(x, q), xytext=(x, p),
                            arrowprops=dict(arrowstyle="->", color=c, lw=1.1))
    ax.axhline(0.5, color="black", lw=0.8, ls=":")
    ax.set_xticks(x0)
    ax.set_xticklabels([{"random": "random\ninactives", "matched": "property-matched\ninactives",
                         "decoy": "DUD-E-protocol\ndecoys"}[s] for s in sets])
    ax.set_ylabel("ROC-AUC (mean over targets)")
    ax.set_ylim(0.4, max(0.75, float(np.nanmax(use[["roc_auc", "conditional_auc"]].to_numpy())) + 0.04))
    ax.set_title("Bars as measured, markers controlled for molecular size")
    ax.legend(fontsize=7, loc="upper right", framealpha=0.9)

    if gaps is not None:
        ax2 = axes[1]
        g = gaps.pivot_table(index="target", columns="negative_set", values="median_gap")
        for i, s in enumerate(sets):
            if s in g.columns:
                ax2.scatter(np.full(len(g), i) + np.linspace(-0.12, 0.12, len(g)), g[s],
                            s=26, color="#444444", alpha=0.8)
        ax2.axhline(0.0, color="black", lw=0.8, ls=":")
        ax2.set_xticks(range(len(sets)))
        ax2.set_xticklabels([s[:8] for s in sets], fontsize=8)
        ax2.set_ylabel("median heavy-atom gap\n(actives - inactives)")
        ax2.set_title("The size confound itself", fontsize=9)
    return _finish(fig, outfile)


def plot_template_difficulty(templates: pd.DataFrame, outfile: Path, threshold: float = 2.0) -> Path:
    """Worst-engine redocking RMSD against the co-crystal ligand's largest ring size.

    Macrocyclic templates are the ones rigid-receptor docking cannot reproduce, which is why a target built
    on them cannot be validated whatever its enrichment.
    """
    d = templates[templates["sanitized"] & np.isfinite(templates["worst_rmsd"])]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for flag, color, label in ((False, "#0072B2", "non-macrocyclic"), (True, "#D55E00", "macrocyclic")):
        sub = d[d["macrocycle"] == flag]
        ax.scatter(sub["max_ring_size"], sub["worst_rmsd"], s=34, color=color, alpha=0.8,
                   edgecolor="black", linewidth=0.4, label=f"{label} (n = {len(sub)})")
    ax.axhline(threshold, color="black", lw=0.9, ls="--")
    ax.annotate(f"validated at {threshold:g} A", xy=(ax.get_xlim()[1], threshold), xytext=(-4, 4),
                textcoords="offset points", ha="right", fontsize=8)
    ax.set_xlabel("largest ring in the co-crystal ligand (atoms)")
    ax.set_ylabel("worst-engine redocking RMSD (A)")
    ax.set_yscale("log")
    ax.legend(fontsize=8)
    return _finish(fig, outfile)
