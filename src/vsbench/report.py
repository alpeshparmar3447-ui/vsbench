"""Render the manuscript tables from the result files, so no number is typed by hand."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from .rescore import REGISTERED_ALIASES

ENGINE_ORDER = ["autodock-gpu", "gnina", "consensus"]


_NEG_LABEL = {"random": "random inactives", "matched": "property-matched inactives",
              "decoy": "DUD-E-protocol decoys"}

# Where a table is published. The main text carries the registered comparison and the controls
# it depends on; everything else is the evidence package, which belongs in the additional file
# rather than being dropped. Numbering follows from this and from position in the spec list, so
# an optional table that is absent cannot leave a gap in the sequence.
MAIN, SUPP = "main", "supp"


def to_markdown(df: pd.DataFrame) -> str:
    cols = [str(c) for c in df.columns]
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join("" if (isinstance(v, float) and np.isnan(v)) else str(v)
                                         for v in row) + " |")
    return "\n".join(lines)


def _ci(row) -> str:
    return f"{row['value']:.2f} [{row['ci_low']:.2f}, {row['ci_high']:.2f}]"


def panel_table(results_dir: Path, targets) -> pd.DataFrame:
    rows = []
    for t in targets:
        tdir = Path(results_dir) / t
        sel = json.loads((tdir / "selected.json").read_text())
        red = pd.read_csv(tdir / "redock.csv")
        chosen = red[red["template"] == sel["template"]].iloc[0]
        scr = pd.read_csv(tdir / "screen_scores.csv")
        rows.append({
            "target": t,
            "actives": int((scr.neg_set == "active").sum()),
            "random inactives": int((scr.neg_set == "random").sum()),
            "matched inactives": int((scr.neg_set == "matched").sum()),
            "decoys": int((scr.neg_set == "decoy").sum()),
            "templates validated": f"{int(red['success'].sum())}/{len(red)}",
            "pocket ions": (str(chosen.get("ions", "")) if str(chosen.get("ions", "")) not in ("", "nan") else "none"),
            "deposition checked": "yes" if bool(chosen.get("metal_audited", False)) else "no",
            "selected": sel["template"],
            "RMSD AutoDock-GPU (A)": f"{chosen['ad_lowest_rmsd']:.2f}",
            "RMSD gnina (A)": f"{chosen['gn_top_rmsd']:.2f}",
            "validated": "yes" if sel["validated"] else "NO",
        })
    return pd.DataFrame(rows)


def auc_table(metrics_all: pd.DataFrame, negative_set: str) -> pd.DataFrame:
    auc = metrics_all[(metrics_all["metric"] == "roc_auc") & (metrics_all["negative_set"] == negative_set)]
    wide = auc.assign(cell=auc.apply(_ci, axis=1)).pivot(index="target", columns="score", values="cell")
    cols = [c for c in ENGINE_ORDER if c in wide.columns]
    return wide[cols].reset_index()


def baseline_table(baselines_all: pd.DataFrame) -> pd.DataFrame:
    sub = baselines_all[baselines_all["property"].isin(["n_heavy", "mw", "tpsa"])]
    sub = sub.assign(cell=sub.apply(lambda r: f"{r['roc_auc']:.2f}", axis=1))
    wide = sub.pivot_table(index="target", columns=["negative_set", "property"], values="roc_auc")
    wide.columns = [f"{n}: {p}" for n, p in wide.columns]
    return wide.round(2).reset_index()


def baseline_comparison_table(results_dir: Path, targets) -> pd.DataFrame | None:
    """Each docking score against the best structure-free baseline, per target and set."""
    rows = []
    for t in targets:
        f = results_dir / t / "analysis" / "baseline_comparison.csv"
        if f.exists():
            rows.append(pd.read_csv(f).assign(target=t))
    if not rows:
        return None
    d = pd.concat(rows, ignore_index=True)
    d = d[d.score.isin(("autodock-gpu", "gnina", "consensus"))]
    out = pd.DataFrame({
        "target": d["target"],
        "negatives": d["negative_set"].map(lambda s: _NEG_LABEL.get(s, s)),
        "score": d["score"],
        "best baseline": d["best_baseline"],
        "baseline AUC": d["baseline_auc"].round(3),
        "score AUC": d["score_auc"].round(3),
        "difference": d["delta_auc"].round(3),
        "95% CI": [f"[{lo:+.3f}, {hi:+.3f}]" for lo, hi in zip(d["ci_low"], d["ci_high"])],
        "adjusted CI": [f"[{lo:+.3f}, {hi:+.3f}]" for lo, hi in zip(d["adj_ci_low"], d["adj_ci_high"])],
        "verdict": ["above" if b else "below" if w else "not resolved"
                    for b, w in zip(d["beats_baseline"], d["below_baseline"])],
    })
    return out.sort_values(["target", "negatives", "score"]).reset_index(drop=True)


def pose_table(results_dir: Path, targets) -> pd.DataFrame:
    rows = []
    for t in targets:
        p = pd.read_csv(Path(results_dir) / t / "analysis" / "pose_selection.csv")
        for _, r in p.iterrows():
            rows.append({"target": t, "negative set": r["negative_set"],
                         "AUC lowest energy": f"{r['auc_lowest']:.2f}",
                         "AUC largest cluster": f"{r['auc_largest']:.2f}",
                         "difference [95% CI]": f"{r['delta_auc_largest_minus_lowest']:+.3f} "
                                                f"[{r['ci_low']:+.3f}, {r['ci_high']:+.3f}]"})
    return pd.DataFrame(rows)


def runtime_table(results_dir: Path, targets) -> pd.DataFrame:
    rows = []
    for t in targets:
        info = json.loads((Path(results_dir) / t / "screen_info.json").read_text())
        n = info["n"]
        ad, gn = info["autodock"], info["gnina"]
        rows.append({
            "target": t, "compounds": n,
            "ligand prep (s/compound)": f"{info['prep_s'] / n:.3f}",
            "AutoDock-GPU (s/compound)": f"{ad['wall_s_cumulative'] / max(ad['n_docked_cumulative'], 1):.2f}",
            "gnina (s/compound)": f"{gn['wall_s_cumulative'] / max(gn['n_docked_cumulative'], 1):.2f}",
        })
    return pd.DataFrame(rows)


def failure_table(results_dir: Path, targets) -> pd.DataFrame:
    rows = []
    for t in targets:
        f = pd.read_csv(Path(results_dir) / t / "analysis" / "failures.csv")
        for _, r in f.iterrows():
            rows.append({"target": t, "class": r["class"], "n": int(r["n"]),
                         "prep failed": int(r["ligand_prep_failed"]),
                         "AutoDock-GPU failed": int(r["autodock_failed"]),
                         "gnina failed": int(r["gnina_failed"]),
                         "gnina positive energy": f"{r['gnina_clash_gt0_frac']:.0%}",
                         "chelator needs manual ionisation": int(r.get("chelator_needs_manual_ionisation", 0)),
                         "chelator deprotonated": int(r.get("chelator_deprotonated", 0))})
    return pd.DataFrame(rows)


def metal_audit_table(audit: pd.DataFrame) -> pd.DataFrame:
    """Per-target summary of the template metal audit (one row per target)."""
    rows = []
    for target, g in audit.drop_duplicates(["target", "template"]).groupby("target"):
        pocket = audit[(audit["target"] == target) & audit["element"].notna()] if "element" in audit else audit.iloc[0:0]
        remote = []
        for v in g.get("nearest_functional_ion_within_30A", pd.Series(dtype=object)).fillna(""):
            if v and "@" in v:
                el, d = v.split("@")
                remote.append((float(d), f"{el} {float(d):.1f} A"))
        rows.append({
            "target": target, "templates": len(g),
            "min frame overlap": f"{g['frame_overlap'].min():.3f}",
            "ions in benchmark files": int(g["ions_in_benchmark_file"].sum()),
            "templates with a pocket ion": int(pocket["template"].nunique()),
            "pocket ion (nearest)": "; ".join(sorted({f"{r.element} {r.dist_to_ligand_A:.1f} A" for r in pocket.itertuples()})),
            "nearest functional ion within 30 A": min(remote)[1] if remote else "",
        })
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ additional tables
def _fmt(series: pd.Series, fmt: str) -> pd.Series:
    """Format numbers; a missing value renders blank, never the string 'nan'."""
    return series.map(lambda v: "" if pd.isna(v) else fmt.format(v))


def _stack(results_dir: Path, targets, fname: str, folder: str = "analysis") -> pd.DataFrame | None:
    frames = []
    for t in targets:
        p = Path(results_dir) / t / folder / fname if folder else Path(results_dir) / t / fname
        if p.exists():
            frames.append(pd.read_csv(p).assign(target=t))
    return pd.concat(frames, ignore_index=True) if frames else None


def negset_table(results_dir: Path, targets) -> pd.DataFrame | None:
    """Change in ROC-AUC between negative-set designs, averaged over targets."""
    d = _stack(results_dir, targets, "negative_set_differences.csv")
    if d is None:
        return None
    d["excludes_zero"] = (d["ci_low"] > 0) | (d["ci_high"] < 0)
    g = d.groupby(["set_a", "set_b", "score"], sort=False)
    out = g.agg(mean_delta_auc=("delta_auc", "mean"), targets=("target", "nunique"),
                targets_with_ci_excluding_zero=("excludes_zero", "sum")).reset_index()
    out["mean_delta_auc"] = _fmt(out["mean_delta_auc"], "{:+.3f}")
    return out.rename(columns={"set_a": "negatives A", "set_b": "negatives B", "score": "score"})


def seed_table(results_dir: Path, targets) -> pd.DataFrame | None:
    d = _stack(results_dir, targets, "seed_stability.csv")
    if d is None:
        return None
    d = d[d["negative_set"] == d["negative_set"].iloc[0]] if d["negative_set"].nunique() else d
    out = d[["target", "n_seeds", "auc_min", "auc_max", "auc_range", "mean_spearman_between_seeds",
             "stable_top_n_fraction", "mean_per_compound_sd"]].copy()
    for c in out.columns[2:]:
        out[c] = _fmt(out[c], "{:.3f}")
    return out


def similarity_table(results_dir: Path, targets) -> pd.DataFrame | None:
    rows = []
    for t in targets:
        p = Path(results_dir) / t / "analysis" / "similarity.json"
        if not p.exists():
            continue
        sim = json.loads(p.read_text())
        a = sim["actives"]
        row = {"target": t, "actives": a["n"], "median NN similarity among actives": f"{a['median_nn_similarity']:.2f}",
               "clusters at 0.4": a["n_clusters_at_0.4"]}
        for name in ("random", "matched", "decoy"):
            if name in sim:
                row[f"analogue AUC vs {name}"] = f"{sim[name]['analogue_baseline_auc']:.2f}"
        rows.append(row)
    return pd.DataFrame(rows) if rows else None


def admet_class_table(results_dir: Path, targets) -> pd.DataFrame | None:
    """Share of each ADMET tier among actives and among the inactive/decoy compounds, per target."""
    rows = []
    for t in targets:
        pa, ps = Path(results_dir) / t / "admet.csv", Path(results_dir) / t / "screen_scores.csv"
        if not (pa.exists() and ps.exists()):
            continue
        a = pd.read_csv(pa, dtype={"name": str})
        sc = pd.read_csv(ps, dtype={"name": str})[["name", "label"]]
        m = a.merge(sc, on="name")
        row = {"target": t}
        for lab, nm in ((1, "actives"), (0, "negatives")):
            g = m[m["label"] == lab]
            for tier in ("PASS", "MODERATE", "HIGH_RISK", "EXCLUDE"):
                row[f"{nm} {tier}"] = f"{(g['tier'] == tier).mean():.0%}"
        rows.append(row)
    return pd.DataFrame(rows) if rows else None


def filter_macro_table(results_dir: Path, targets) -> pd.DataFrame | None:
    d = _stack(results_dir, targets, "admet_filters.csv")
    if d is None:
        return None
    g = d.groupby(["negative_set", "filter"], sort=False)
    out = g.agg(actives_retained=("active_retained", "mean"), inactives_retained=("inactive_retained", "mean"),
                actives_lost=("actives_lost", "sum"), auc_before=("roc_auc_before", "mean"),
                auc_after=("roc_auc_after", "mean"), targets=("target", "nunique")).reset_index()
    for c in ("actives_retained", "inactives_retained"):
        out[c] = _fmt(out[c], "{:.0%}")
    for c in ("auc_before", "auc_after"):
        out[c] = _fmt(out[c], "{:.3f}")
    return out


def prioritization_table(results_dir: Path, targets) -> pd.DataFrame | None:
    d = _stack(results_dir, targets, "prioritization.csv")
    if d is None:
        return None
    g = d.groupby(["negative_set", "strategy", "top_n"], sort=False)
    out = g.agg(hit_rate=("hit_rate", "mean"), sd_across_targets=("hit_rate", "std"),
                pool_active_fraction=("pool_active_fraction", "mean"), targets=("target", "nunique")).reset_index()
    for c in ("hit_rate", "sd_across_targets", "pool_active_fraction"):
        out[c] = _fmt(out[c], "{:.3f}")
    return out


def interaction_table(results_dir: Path, targets) -> pd.DataFrame | None:
    rows = []
    for t in targets:
        pi, ps = Path(results_dir) / t / "interactions.csv", Path(results_dir) / t / "screen_scores.csv"
        if not (pi.exists() and ps.exists()):
            continue
        i = pd.read_csv(pi, dtype={"name": str})
        sc = pd.read_csv(ps, dtype={"name": str})[["name", "label"]]
        m = i.merge(sc, on="name")
        m = m[m["ok"].astype(bool)]
        row = {"target": t, "crystal contacts": int(i["n_reference"].iloc[0])}
        for lab, nm in ((1, "actives"), (0, "top-ranked negatives")):
            g = m[m["label"] == lab]
            row[f"recall, {nm}"] = f"{g['recall'].mean():.2f}" if len(g) else ""
            row[f"n, {nm}"] = len(g)
        rows.append(row)
    return pd.DataFrame(rows) if rows else None


def structure_qc_table(qc: pd.DataFrame) -> pd.DataFrame:
    g = qc.groupby("target")
    return g.agg(templates=("template", "count"), median_resolution=("resolution", "median"),
                 with_chain_break=("n_unexplained_breaks", lambda v: int((v > 0).sum())),
                 multiple_chains=("flags", lambda v: int(v.fillna("").str.contains("multiple_chains").sum())),
                 engineered_mutation=("engineered_mutation", lambda v: int((v > 0).sum())),
                 median_rama_core=("rama_core", "median")).round(3).reset_index()


def site_summary_table(sites: pd.DataFrame) -> pd.DataFrame:
    """Top-1 and top-3 site recovery per detector, with Wilson intervals."""
    from .gridstrat import wilson
    rows = []
    for method, g in sites.groupby("method", sort=False):
        n = len(g)
        row = {"detector": method, "templates": n}
        for col, lab in (("top1_hit", "top-1"), ("top3_hit", "top-3")):
            k = int(g[col].astype(bool).sum())
            lo, hi = wilson(k, n)
            row[f"{lab} recovery"] = f"{k}/{n} ({k / n:.0%}) [{lo:.0%}, {hi:.0%}]"
        row["mean best-pocket recovery"] = f"{g['best_recovery'].mean():.2f}"
        row["median pockets"] = int(g["n_pockets"].median())
        rows.append(row)
    return pd.DataFrame(rows)


def site_control_table(controls: pd.DataFrame) -> pd.DataFrame:
    out = controls[["template", "control_role", "description", "method", "rank_of_site", "site_pocket_druggability",
                    "ligand_buried_fraction"]].copy()
    out["site_pocket_druggability"] = _fmt(out["site_pocket_druggability"], "{:.2f}")
    out["ligand_buried_fraction"] = _fmt(out["ligand_buried_fraction"], "{:.2f}")
    return out


def table_from_csv(path: Path, columns: list[str] | None = None, fmt: dict | None = None) -> pd.DataFrame | None:
    path = Path(path)
    if not path.exists():
        return None
    d = pd.read_csv(path)
    if columns:
        d = d[[c for c in columns if c in d.columns]]
    for c, f in (fmt or {}).items():
        if c in d:
            d[c] = d[c].map(f.format)
    # A CSV carries full float precision; a table does not. Printing -7.402666666666666 in a
    # six-column-wide cell wraps it over three lines and reads as noise, so columns not given
    # an explicit format are rounded to a readable number of significant figures.
    explicit = set(fmt or {})
    for c in d.columns:
        if c in explicit or not pd.api.types.is_float_dtype(d[c]):
            continue
        finite = d[c].abs().replace(0, np.nan).dropna()
        if finite.empty:
            continue
        # keep small quantities (p-values, fractions) usefully precise, large ones tidy
        decimals = 4 if finite.min() < 0.01 else 3
        d[c] = d[c].round(decimals)
    return d


# ------------------------------------------------- Amendment 3 tables
def size_control_table(results_dir: Path, negative_set: str) -> pd.DataFrame | None:
    """Plain against size-controlled AUC, averaged over targets, for one negative set."""
    p = Path(results_dir) / "size_control_macro.csv"
    if not p.exists():
        return None
    d = pd.read_csv(p)
    # A result set written before the alias fix carries one score twice, once under its
    # registered name and once under a variant label for the same column. Keep the registered
    # row, so the table lists each score once whichever version wrote the file.
    present = set(d["score"])
    d = d[~d["score"].isin([a for a, reg in REGISTERED_ALIASES.items() if reg in present])]
    d = d[d["negative_set"] == negative_set].sort_values("conditional_auc", ascending=False)
    if not len(d):
        return None
    return pd.DataFrame({
        "score": d["score"],
        "ROC-AUC": _fmt(d["roc_auc"], "{:.3f}"),
        "size-controlled AUC": _fmt(d["conditional_auc"], "{:.3f}"),
        "change": _fmt(d["mean_shift"], "{:+.3f}"),
        "SD of change": _fmt(d["sd_shift"], "{:.3f}"),
        "targets": d["n_targets"],
    }).reset_index(drop=True)


def size_gap_table(results_dir: Path) -> pd.DataFrame | None:
    """Median heavy-atom gap between actives and each negative set: the size confound itself."""
    p = Path(results_dir) / "size_profile.csv"
    if not p.exists():
        return None
    d = pd.read_csv(p)
    w = d.pivot_table(index="target", columns="negative_set", values="median_gap").reset_index()
    for c in w.columns:
        if c != "target":
            w[c] = _fmt(w[c], "{:+.1f}")
    return w


def score_variant_table(results_dir: Path, negative_set: str, metric: str = "roc_auc") -> pd.DataFrame | None:
    """Every score variant against the registered consensus, averaged over targets."""
    p = Path(results_dir) / "score_variants_macro.csv"
    if not p.exists():
        return None
    d = pd.read_csv(p)
    d = d[(d["negative_set"] == negative_set) & (d["metric"] == metric)].sort_values("value", ascending=False)
    if not len(d):
        return None
    return pd.DataFrame({
        "score": d["score"],
        metric: _fmt(d["value"], "{:.3f}" if metric == "roc_auc" else "{:.2f}"),
        "vs registered consensus": _fmt(d["mean_vs_primary"], "{:+.3f}" if metric == "roc_auc" else "{:+.2f}"),
        # The adjusted counts govern every claim of improvement; the pointwise ones are shown
        # beside them so the reader can see how many apparent wins were multiple comparisons.
        "targets better": d.get("n_better_adjusted", d["n_better"]).astype(int),
        "targets worse": d.get("n_worse_adjusted", d["n_worse"]).astype(int),
        "better unadjusted": d["n_better"].astype(int),
        "targets": d["n_targets"].astype(int),
    }).reset_index(drop=True)


def template_difficulty_table(results_dir: Path) -> pd.DataFrame | None:
    """Per target: how many templates validated, and how large their co-crystal ligands are."""
    p = Path(results_dir) / "template_difficulty_per_target.csv"
    if not p.exists():
        return None
    d = pd.read_csv(p)
    return pd.DataFrame({
        "target": d["target"], "templates": d["templates"], "redocked": d["ran"],
        "validated (<=2 A)": d["validated"], "macrocyclic ligands": d["macrocyclic"],
        "median heavy atoms": _fmt(d["median_heavy_atoms"], "{:.0f}"),
        "median MW": _fmt(d["median_mol_weight"], "{:.0f}"),
        "best RMSD (A)": _fmt(d["best_rmsd"], "{:.2f}"),
    })


def template_property_table(results_dir: Path) -> pd.DataFrame | None:
    """Which ligand property predicts redocking failure."""
    p = Path(results_dir) / "template_difficulty_correlations.csv"
    if not p.exists():
        return None
    d = pd.read_csv(p).sort_values("spearman_rho", ascending=False)
    names = {"heavy_atoms": "heavy atoms", "mol_weight": "molecular weight",
             "rotatable_bonds": "rotatable bonds", "n_rings": "number of rings",
             "max_ring_size": "largest ring size"}
    return pd.DataFrame({
        "ligand property": d["property"].map(lambda v: names.get(v, v)),
        "Spearman rho vs RMSD": _fmt(d["spearman_rho"], "{:+.3f}"),
        "p (permutation)": _fmt(d["p_permutation"], "{:.3f}"),
        "templates": d["n_templates"].astype(int),
    }).reset_index(drop=True)


def macrocycle_table(results_dir: Path) -> pd.DataFrame | None:
    p = Path(results_dir) / "template_difficulty_by_macrocycle.csv"
    if not p.exists():
        return None
    d = pd.read_csv(p)
    return pd.DataFrame({
        "template group": d["group"], "templates": d["n_templates"].astype(int),
        "validated": d["validated"].astype(int),
        "validation rate": _fmt(d["rate"], "{:.1%}"),
        "95% CI": [("" if pd.isna(lo) else f"{lo:.1%} to {hi:.1%}") for lo, hi in zip(d["ci_low"], d["ci_high"])],
    })


def box_diagnostic_table(results_dir: Path) -> pd.DataFrame | None:
    """Redocking RMSD against box size for any template put through the box diagnostic."""
    files = sorted(Path(results_dir).glob("box_diagnostic_*.csv"))
    if not files:
        return None
    d = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    return pd.DataFrame({
        "target": d["target"], "template": d["template"],
        "box edge (A)": _fmt(d["box_edge"], "{:.0f}"),
        "AutoDock-GPU RMSD (A)": _fmt(d["ad_rmsd"], "{:.2f}"),
        "gnina RMSD (A)": _fmt(d["gn_rmsd"], "{:.2f}"),
        "best gnina pose (A)": _fmt(d["gn_best_rmsd"], "{:.2f}"),
        "recovered": d["recovered"].map({True: "yes", False: "no"}),
    })

# ------------------------------------------------- merged across negative set
# One analysis reported per negative set is one finding, not three. These wrappers stack the
# per-set tables into a single block with a negative-set column, which is what a reader
# compares across; the per-set functions stay as they are because the figures use them.
def _stack_by_negative_set(make, sets=("random", "matched", "decoy")) -> pd.DataFrame | None:
    frames = []
    for ns in sets:
        t = make(ns)
        if t is not None and len(t):
            frames.append(t.assign(**{"negative set": _NEG_LABEL[ns]}))
    if not frames:
        return None
    d = pd.concat(frames, ignore_index=True)
    return d[["negative set"] + [c for c in d.columns if c != "negative set"]]


def auc_table_all(metrics_all: pd.DataFrame) -> pd.DataFrame | None:
    """ROC-AUC against every negative set in one table."""
    return _stack_by_negative_set(
        lambda ns: auc_table(metrics_all, ns) if (metrics_all["negative_set"] == ns).any() else None)


def size_control_table_all(results_dir: Path) -> pd.DataFrame | None:
    """Size-controlled ROC-AUC against every negative set in one table."""
    return _stack_by_negative_set(lambda ns: size_control_table(results_dir, ns))


def score_variant_table_all(results_dir: Path, metric: str = "roc_auc") -> pd.DataFrame | None:
    """Score variants against every matched design in one table, for one metric.

    Kept per metric: ROC-AUC and EF1 are different scales, so a shared value column would
    invite comparing numbers that are not comparable.
    """
    return _stack_by_negative_set(lambda ns: score_variant_table(results_dir, ns, metric),
                                  sets=("matched", "decoy"))


def number_tables(entries):
    """Label tables by position: main text ``Table N``, additional file ``Table SN``."""
    out, n, s = [], 0, 0
    for placement, title, df in entries:
        if placement == MAIN:
            n += 1
            label = f"Table {n}"
        else:
            s += 1
            label = f"Table S{s}"
        out.append((placement, f"{label}. {title}", df))
    return out


def build_report(results_dir: Path, targets, audit_csv: Path | None = None) -> str:
    """All tables as one markdown document (written to ``tables.md``).

    Main-text and supplementary tables are also written separately, to
    ``tables_main.md`` and ``tables_supplementary.md``.
    """
    results_dir = Path(results_dir)
    m = pd.read_csv(results_dir / "metrics_all_targets.csv")
    b = pd.concat([pd.read_csv(results_dir / t / "analysis" / "baselines.csv").assign(target=t)
                   for t in targets], ignore_index=True)
    top = results_dir.parent

    def own(fname, **kw):
        """A file this run produced. Scoped to the run directory: a same-named file in the
        parent would be another run's, and silently reporting it is the failure this avoids."""
        f = results_dir / fname
        return table_from_csv(f, **kw) if f.exists() else None

    def shared(fname, fn=None):
        """A study output, which may sit beside the run or one level up (studies covering
        every target write there)."""
        for base in (results_dir, top):
            if (base / fname).exists():
                return fn(pd.read_csv(base / fname)) if fn else table_from_csv(base / fname)
        return None

    audit = Path(audit_csv) if audit_csv else top / "metal_audit.csv"
    spec = [
        # ---- main text: the registered comparison and the controls it rests on
        (MAIN, "Panel, template validation and compounds docked",
         panel_table(results_dir, targets)),
        (MAIN, "ROC-AUC against each negative set", auc_table_all(m)),
        (MAIN, "ROC-AUC of single properties used as scores", baseline_table(b)),
        (MAIN, "Size-controlled ROC-AUC against each negative set",
         size_control_table_all(results_dir)),
        (MAIN, "Docking against the best structure-free baseline",
         baseline_comparison_table(results_dir, targets)),
        (MAIN, "Resampling sensitivity: compound, matched-group and scaffold-aware bootstraps",
         own("bootstrap_sensitivity.csv")),
        (MAIN, "Retrospective and deployable ligand-similarity baselines",
         own("baseline_leakage.csv")),
        (MAIN, "Docking-box strategies", shared("grid_strategies_summary.csv")),
        (MAIN, "Docking-box strategies, paired per template (exact McNemar)",
         shared("grid_strategies_paired.csv")),
        # ---- additional file: the rest of the evidence package, disclosed not deleted
        (SUPP, "Pose selection", pose_table(results_dir, targets)),
        (SUPP, "Throughput", runtime_table(results_dir, targets)),
        (SUPP, "Throughput, CPU against GPU", shared("perf_summary.csv")),
        (SUPP, "Failures", failure_table(results_dir, targets)),
        (SUPP, "Metal-ion audit of the templates against their original depositions",
         metal_audit_table(pd.read_csv(audit)) if audit.exists() else None),
        (SUPP, "Change in ROC-AUC between negative-set designs",
         negset_table(results_dir, targets)),
        (SUPP, "Seed stability of AutoDock-GPU", seed_table(results_dir, targets)),
        (SUPP, "Similarity of the negative sets to the actives",
         similarity_table(results_dir, targets)),
        (SUPP, "Median heavy-atom gap, actives against each negative set",
         size_gap_table(results_dir)),
        (SUPP, "ADMET risk tiers", admet_class_table(results_dir, targets)),
        (SUPP, "Effect of ADMET filters on retrieval", filter_macro_table(results_dir, targets)),
        (SUPP, "Hit rate of the lead-selection strategies",
         prioritization_table(results_dir, targets)),
        (SUPP, "Recovery of the co-crystal interactions by the top gnina pose",
         interaction_table(results_dir, targets)),
        (SUPP, "Structure quality of the templates", shared("structure_qc.csv", structure_qc_table)),
        (SUPP, "Pocket-detector recovery of the ligand site", shared("sites.csv", site_summary_table)),
        (SUPP, "Druggability controls", shared("site_controls.csv", site_control_table)),
        (SUPP, "Quality-gate ablation", shared("ablation_summary.csv")),
        (SUPP, "Template difficulty by target", template_difficulty_table(results_dir)),
        (SUPP, "Which co-crystal ligand property predicts redocking failure",
         template_property_table(results_dir)),
        (SUPP, "Redocking validation rate by macrocyclicity", macrocycle_table(results_dir)),
        (SUPP, "Redocking RMSD against box size", box_diagnostic_table(results_dir)),
        (SUPP, "Score variants by ROC-AUC against each matched design",
         score_variant_table_all(results_dir, "roc_auc")),
        (SUPP, "Score variants by EF1 against each matched design",
         score_variant_table_all(results_dir, "ef1")),
        # Every baseline separately, so the selection of the strongest can be checked rather
        # than trusted; the primary claim is read from this table, not from the selected winner.
        (SUPP, "Docking against each structure-free baseline separately",
         own("baseline_comparison_all.csv")),
    ]
    numbered = number_tables([(p, t, d) for p, t, d in spec if d is not None and len(d)])

    def render(placement):
        return "\n\n".join(f"**{t}**\n\n{to_markdown(d)}"
                           for p, t, d in numbered if p == placement) + "\n"

    main_text, supp_text = render(MAIN), render(SUPP)
    (results_dir / "tables_main.md").write_text(main_text)
    (results_dir / "tables_supplementary.md").write_text(supp_text)
    text = ("## Main-text tables\n\n" + main_text
            + "\n## Additional file: supplementary tables\n\n" + supp_text)
    (results_dir / "tables.md").write_text(text)
    return text
