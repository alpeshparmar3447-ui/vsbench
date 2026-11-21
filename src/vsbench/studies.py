"""Studies that run on top of the screening results: gate ablation, docking-box strategies, performance,
and the figure set. Each writes CSV files under ``results/`` and is resumable where it is expensive."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

from . import ablation, data, engines, gridstrat, ligprep, perf, receptors, redock, sites
from .box import box_from_coords
from .pipeline import Config


def _selected(cfg: Config, target_name: str) -> str:
    return json.loads((cfg.target_dir(target_name) / "selected.json").read_text())["template"]


def _receptor(cfg: Config, target_name: str):
    """Target, selected template, prepared receptor and metal plan: exactly as the screening stage built them."""
    template = _selected(cfg, target_name)
    target = data.load_target(Path(cfg.data_root), target_name)
    tdir = cfg.target_dir(target_name)
    xtal, _ = redock.crystal_ligand(target.ligand_mol2(template), tdir / "redock" / template)
    heavy = Chem.RemoveHs(xtal).GetConformer().GetPositions()
    box = box_from_coords(heavy, **cfg.box)
    plan = redock.plan_for(target.protein_mol2(template), heavy, box, template, cfg.metals)
    rec = engines.prepare_receptor(target.protein_mol2(template), tdir / "redock" / template / "receptor", box,
                                   keep_ions=plan.ions)
    return target, template, rec, plan


def _sample_ligands(cfg: Config, target_name: str, n_each: int, plan, neg_set: str = "matched"):
    scr = pd.read_csv(cfg.target_dir(target_name) / "screen_scores.csv", dtype={"name": str})
    act = scr[scr["neg_set"] == "active"].sample(min(n_each, int((scr["neg_set"] == "active").sum())), random_state=cfg.seed)
    neg = scr[scr["neg_set"] == neg_set].sample(min(n_each, int((scr["neg_set"] == neg_set).sum())), random_state=cfg.seed)
    pool = pd.concat([act, neg])
    return ligprep.prepare_many(list(zip(pool["smiles"], pool["name"])), cfg.seed, cfg.n_jobs, chelate_zn=plan.zn)


# -------------------------------------------------------------------- ablation
def ablation_study(cfg: Config, targets, n_each: int = 30, log=print) -> pd.DataFrame:
    out_dir = Path(cfg.results_dir) / "ablation"
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for t in targets:
        _, template, rec, plan = _receptor(cfg, t)
        ligs = _sample_ligands(cfg, t, n_each, plan)
        wide, summ = ablation.run_ablation(rec, ligs, out_dir / t, nrun=cfg.autodock["nrun"], seed=cfg.seed,
                                           n_boot=cfg.n_boot)
        wide.to_csv(out_dir / f"{t}_energies.csv")
        rows.append(summ.assign(target=t, template=template, n_ligands=len(ligs)))
        log(f"[{t}] ablation on {len(ligs)} ligands")
    table = pd.concat(rows, ignore_index=True)
    table.to_csv(Path(cfg.results_dir) / "ablation_summary.csv", index=False)
    return table


# -------------------------------------------------------------- grid strategies
def _protein_xyz(pdb: Path) -> np.ndarray:
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                     for l in Path(pdb).read_text().splitlines() if l.startswith("ATOM")])


def grid_study(cfg: Config, targets, log=print) -> pd.DataFrame:
    """Redock every template's co-crystal ligand with the ligand, pocket and blind boxes."""
    out_dir = Path(cfg.results_dir) / "grid"
    rows_path = Path(cfg.results_dir) / "grid_strategies.csv"
    rows = pd.read_csv(rows_path, dtype={"template": str}).to_dict("records") if rows_path.exists() else []
    done = {(r["template"], r["strategy"]) for r in rows}
    kw = dict(seed=cfg.seed, nrun=cfg.autodock["nrun"], exhaustiveness=cfg.gnina["exhaustiveness"],
              num_modes=cfg.gnina["num_modes"], cpu=cfg.gnina["cpu"], box_kwargs=cfg.box, metals_cfg=cfg.metals)
    for t in targets:
        target = data.load_target(Path(cfg.data_root), t)
        wanted = (cfg.templates or {}).get(t) or list(target.templates)
        for tpl in wanted:
            work = out_dir / t / tpl
            pdb = work / "protein.pdb"
            work.mkdir(parents=True, exist_ok=True)
            sites.write_pdb_from_mol2(target.protein_mol2(tpl), pdb)
            boxes = {"ligand": None}
            try:
                centre = gridstrat.top_pocket_center(sites.run_fpocket(pdb, work / "fpocket"))
                boxes["pocket"] = gridstrat.pocket_box(centre) if centre else None
            except Exception as exc:
                log(f"[{t}/{tpl}] fpocket failed: {exc}")
                boxes["pocket"] = None
            # The same box seeded from the learned detector instead of the geometric one, so the
            # comparison measures box seeding rather than one detector's ranking accuracy.
            opt = receptors.options(cfg)
            if opt["p2rank"] and opt.get("p2rank_bin") and opt.get("java_bin"):
                try:
                    pk = sites.run_p2rank(pdb, work / "p2rank", Path(opt["p2rank_bin"]),
                                          Path(opt["java_bin"]))
                    centre = gridstrat.top_pocket_center(pk)
                    boxes["p2rank"] = gridstrat.pocket_box(centre) if centre else None
                except Exception as exc:
                    log(f"[{t}/{tpl}] p2rank failed: {exc}")
                    boxes["p2rank"] = None
            boxes["blind"] = gridstrat.blind_box(_protein_xyz(pdb))
            for strat, box in boxes.items():
                if (tpl, strat) in done:
                    continue
                if strat != "ligand" and box is None:
                    rows.append({"target": t, "template": tpl, "strategy": strat, "ok": False,
                                 "error": "no pocket centre", "ad_success": False, "gn_success": False,
                                 "ad_lowest_rmsd": np.nan, "gn_top_rmsd": np.nan, "wall_s": np.nan})
                else:
                    r = gridstrat.run_strategy(strat, target.protein_mol2(tpl), target.ligand_mol2(tpl), tpl,
                                               work / strat, box, pdb_id=tpl, **kw)
                    rows.append({"target": t, **r})
                log(f"[{t}/{tpl}] {strat}: {rows[-1].get('gn_top_rmsd')}")
                pd.DataFrame(rows).to_csv(rows_path, index=False)
    table = pd.DataFrame(rows)
    table.to_csv(rows_path, index=False)
    gridstrat.summarize(rows).to_csv(Path(cfg.results_dir) / "grid_strategies_summary.csv", index=False)
    # Unpaired success rates cannot support a claim that one box beats another, so the paired
    # test is written beside them rather than left to be computed by hand.
    paired = pd.concat([gridstrat.paired_strategies(rows, scope=s)
                        for s in ("pairwise", "listwise")], ignore_index=True)
    paired.to_csv(Path(cfg.results_dir) / "grid_strategies_paired.csv", index=False)
    return table


# ---------------------------------------------------------------- performance
def perf_study(cfg: Config, target_name: str, n_values=(8, 16, 32), cpu_threads=(1, 4, 8, 16), cpu_n: int = 16,
               repeats: int = 3, log=print) -> pd.DataFrame:
    _, _, rec, plan = _receptor(cfg, target_name)
    ligs = [l for l in _sample_ligands(cfg, target_name, max(n_values), plan) if l.ok]
    plan_cases = [(c, n, cfg.gnina["cpu"]) for c in ("autodock-gpu", "gnina (CNN, GPU)") for n in n_values if n <= len(ligs)]
    for c in ("gnina (CNN, CPU)", "gnina (Vina, CPU)", "vina (CPU)"):
        plan_cases += [(c, min(cpu_n, len(ligs)), th) for th in cpu_threads]
    runs = perf.run_suite(rec, ligs, Path(cfg.results_dir) / "perf_work", plan_cases, repeats=repeats, log=log,
                          seed=cfg.seed, nrun=cfg.autodock["nrun"], exhaustiveness=cfg.gnina["exhaustiveness"])
    runs.to_csv(Path(cfg.results_dir) / "perf_runs.csv", index=False)
    summ = perf.summarize(runs)
    summ.to_csv(Path(cfg.results_dir) / "perf_summary.csv", index=False)
    return summ


# --------------------------------------------------------------------- figures
def make_figures(cfg: Config, targets=None) -> list[Path]:
    """Render every figure that the available result files support."""
    from . import figures as F
    res, out = Path(cfg.results_dir), Path(cfg.results_dir) / "figures"
    targets = targets or cfg.targets
    made = []
    mpath = res / "metrics_all_targets.csv"
    if mpath.exists():
        m = pd.read_csv(mpath)
        b = pd.concat([pd.read_csv(res / t / "analysis" / "baselines.csv").assign(target=t)
                       for t in targets if (res / t / "analysis" / "baselines.csv").exists()], ignore_index=True)
        made.append(F.plot_auc_summary(m, b if len(b) else None, out / "auc_summary.png"))
        made.append(F.plot_negative_set_auc(m, "consensus", out / "negative_sets_consensus.png"))
    reds = [pd.read_csv(res / t / "redock.csv").assign(target=t) for t in targets if (res / t / "redock.csv").exists()]
    if reds:
        made.append(F.plot_redock(pd.concat(reds, ignore_index=True), out / "redocking.png"))
    pr = [pd.read_csv(res / t / "analysis" / "prioritization.csv").assign(target=t)
          for t in targets if (res / t / "analysis" / "prioritization.csv").exists()]
    if pr:
        made.append(F.plot_prioritization(pd.concat(pr, ignore_index=True), out / "prioritization.png"))
    sc = res / "size_control_macro.csv"
    if sc.exists():
        gaps = res / "size_profile.csv"
        made.append(F.plot_size_control(pd.read_csv(sc), pd.read_csv(gaps) if gaps.exists() else None,
                                        out / "size_control.png"))
    td = res / "template_difficulty.csv"
    if td.exists():
        made.append(F.plot_template_difficulty(pd.read_csv(td), out / "template_difficulty.png"))
    for base in (res, res.parent):
        if (base / "sites.csv").exists():
            made.append(F.plot_site_recovery(pd.read_csv(base / "sites.csv"), out / "site_recovery.png"))
            break
    if (res / "grid_strategies_summary.csv").exists():
        made.append(F.plot_grid_strategies(pd.read_csv(res / "grid_strategies_summary.csv"), out / "grid_strategies.png"))
    if (res / "perf_summary.csv").exists():
        made.append(F.plot_performance(pd.read_csv(res / "perf_summary.csv"), out / "performance.png"))
    return made


# ------------------------------------------- post-screen analyses (Amendment 3)
def _engine_rmsd(cfg: Config, target_name: str) -> tuple[float, float]:
    """In-place redocking RMSD of each engine on the template the campaign actually used."""
    tdir = cfg.target_dir(target_name)
    sel = json.loads((tdir / "selected.json").read_text())["template"]
    red = pd.read_csv(tdir / "redock.csv")
    row = red[red["template"].astype(str) == str(sel)]
    if not len(row):
        return float("nan"), float("nan")
    return float(row.iloc[0]["ad_lowest_rmsd"]), float(row.iloc[0]["gn_top_rmsd"])


def _scored_with_variants(cfg: Config, target_name: str) -> pd.DataFrame:
    """The frozen screening table with the registered scores, the unevaluated scores and the fusion
    variants, plus interaction recovery as a score where it has been computed for every compound."""
    from . import analysis, rescore
    tdir = cfg.target_dir(target_name)
    df = analysis.add_scores(pd.read_csv(tdir / "screen_scores.csv"))
    ad_rmsd, gn_rmsd = _engine_rmsd(cfg, target_name)
    out = rescore.build_variants(df, rescore.redocking_weights(ad_rmsd, gn_rmsd))
    ipath = tdir / "interactions.csv"
    if ipath.exists():
        ifp = pd.read_csv(ipath, dtype={"name": str})
        have = set(ifp["name"])
        if set(out["name"].astype(str)) <= have:                       # only a full pass can rank a whole set
            rec = ifp.set_index("name")["recall"].astype(float)
            vals = out["name"].astype(str).map(rec).fillna(0.0)        # unanalysable pose recovers nothing
            out["score_interaction_recall"] = vals.to_numpy()
    return out


def variant_macro(allv: pd.DataFrame) -> pd.DataFrame:
    """Per-target variant comparisons averaged across targets.

    Counts the targets where a variant beat the primary both ways. The adjusted count is the one
    the promotion rule reads; the pointwise count is kept beside it because the gap between them
    is itself worth seeing, being the number of apparent wins that are multiple comparisons.
    """
    agg = dict(value=("value", "mean"), sd_across_targets=("value", "std"),
               mean_vs_primary=("vs_primary", "mean"),
               n_better=("better", "sum"), n_worse=("worse", "sum"),
               n_targets=("target", "nunique"))
    if "better_adjusted" in allv.columns:
        agg["n_better_adjusted"] = ("better_adjusted", "sum")
        agg["n_worse_adjusted"] = ("worse_adjusted", "sum")
    return allv.groupby(["score", "negative_set", "metric"], as_index=False).agg(**agg)


def rescore_study(cfg: Config, targets=None, n_boot: int | None = None,
                  variant_n_boot: int | None = None, log=print) -> dict:
    """Items A to D of Amendment 3: template difficulty, size-controlled AUC, and the score variants.

    Re-docks nothing. Every number comes from files the screening stage already wrote.

    ``variant_n_boot`` applies to the variant comparison alone. That comparison is the only one
    carrying a multiplicity-adjusted interval, which reads a far quantile and so needs more draws
    than the pointwise intervals elsewhere; raising the global ``n_boot`` to suit it would pay that
    cost across every table that does not need it.

    Both counts fall back to the config, not to each other. ``variant_n_boot`` used to default to
    whatever ``n_boot`` resolved to, so ``cfg.variant_n_boot`` - which exists, defaults to 5000 and
    is read correctly by ``baseline_comparison`` - was accepted, written into every manifest and
    then ignored here. A run that did not pass ``--variant-n-boot`` on the command line silently
    compared the variants at 1000 draws, which does not resolve the adjusted interval for a family
    of this size (``draws_resolve_interval`` is False below 2400 for twelve comparisons). The
    registered arm passed the flag explicitly and so was unaffected, but nothing in the config made
    that necessary and nothing failed when it was omitted.
    """
    from . import analysis, metrics, rescore, sizecontrol, templatestats
    res = Path(cfg.results_dir)
    targets = list(targets or cfg.targets)
    n_boot = cfg.n_boot if n_boot is None else n_boot
    variant_n_boot = cfg.variant_n_boot if variant_n_boot is None else variant_n_boot

    # A. template difficulty, across every target that was redocked
    tt = templatestats.template_table(Path(cfg.data_root), res, targets)
    tt.to_csv(res / "template_difficulty.csv", index=False)
    templatestats.correlations(tt).to_csv(res / "template_difficulty_correlations.csv", index=False)
    templatestats.by_macrocycle(tt).to_csv(res / "template_difficulty_by_macrocycle.csv", index=False)
    templatestats.per_target(tt).to_csv(res / "template_difficulty_per_target.csv", index=False)
    log(f"template difficulty: {len(tt)} templates, "
        f"{int((tt['macrocycle'] == True).sum())} macrocyclic, {int(tt['validated'].sum())} validated")

    var_frames, size_frames, strata_frames, profile_rows = [], [], [], []
    for t in targets:
        if not (cfg.target_dir(t) / "screen_scores.csv").exists():
            continue
        df = _scored_with_variants(cfg, t)
        adir = cfg.target_dir(t) / "analysis"
        adir.mkdir(parents=True, exist_ok=True)
        cols = dict(analysis.SCORE_COLUMNS)
        # A variant label can point at a column a registered score already covers ("gnina alone"
        # is the registered "gnina"). Carrying both would put one score in the size-control table
        # twice under two names, so the registered label wins and the alias is dropped.
        taken = set(cols.values())
        for label, col in rescore.variant_columns(df).items():
            if col not in taken:
                cols[label] = col
                taken.add(col)
        if "score_interaction_recall" in df.columns:
            cols["interaction recall"] = "score_interaction_recall"
        # The family is what compare() adjusts over: the variants, not the registered scores that
        # the size-control table also carries.
        n_variants = len(rescore.variant_columns(df))
        if not metrics.draws_resolve_interval(variant_n_boot, n_comparisons=n_variants):
            log(f"[{t}] note: {variant_n_boot} draws leave the multiplicity-adjusted interval for "
                f"{n_variants} variants at the mercy of a couple of extreme replicates; the pointwise "
                f"interval is unaffected and no variant can be promoted on an unresolved one")
        v, s, st = [], [], []
        for name, sub in analysis.evaluation_sets(df).items():
            v.append(rescore.compare(sub, name, n_boot=variant_n_boot, seed=cfg.seed))
            a, b = sizecontrol.table(sub, name, cols, n_boot=n_boot, seed=cfg.seed)
            s.append(a)
            st.append(b)
            profile_rows.append({"target": t, **sizecontrol.size_profile(sub, name)})
        pd.concat(v, ignore_index=True).to_csv(adir / "score_variants.csv", index=False)
        pd.concat(s, ignore_index=True).to_csv(adir / "size_control.csv", index=False)
        pd.concat(st, ignore_index=True).to_csv(adir / "size_control_strata.csv", index=False)
        var_frames.append(pd.concat(v, ignore_index=True).assign(target=t))
        size_frames.append(pd.concat(s, ignore_index=True).assign(target=t))
        strata_frames.append(pd.concat(st, ignore_index=True).assign(target=t))
        log(f"[{t}] score variants and size-controlled AUC written")

    out = {}
    if var_frames:
        allv = pd.concat(var_frames, ignore_index=True)
        allv.to_csv(res / "score_variants_all_targets.csv", index=False)
        macro = variant_macro(allv)
        macro.to_csv(res / "score_variants_macro.csv", index=False)
        out["variants"] = macro
        prom = rescore.promoted(allv)
        prom.to_csv(res / "score_variants_promoted.csv", index=False)
        log(f"variants beating the registered consensus on both hard negative sets: "
            f"{sorted(set(prom['score'])) or 'none'}")
    if size_frames:
        alls = pd.concat(size_frames, ignore_index=True)
        alls.to_csv(res / "size_control_all_targets.csv", index=False)
        macro = (alls.groupby(["score", "negative_set"], as_index=False)
                 .agg(roc_auc=("roc_auc", "mean"), conditional_auc=("conditional_auc", "mean"),
                      mean_shift=("shift", "mean"), sd_shift=("shift", "std"),
                      n_targets=("target", "nunique")))
        macro.to_csv(res / "size_control_macro.csv", index=False)
        pd.concat(strata_frames, ignore_index=True).to_csv(res / "size_control_strata_all_targets.csv", index=False)
        pd.DataFrame(profile_rows).to_csv(res / "size_profile.csv", index=False)
        out["size_control"] = macro
    out["templates"] = tt
    return out


def interaction_full_study(cfg: Config, targets=None, log=print) -> dict:
    """Item D's prerequisite: interaction recovery for every screened compound, not just the pools.

    Uses the poses the screen already produced; PLIP runs on the CPU, so this adds no GPU load.
    """
    from . import analysis, pipeline
    targets = list(targets or cfg.targets)
    counts = {}
    for t in targets:
        tdir = cfg.target_dir(t)
        if not (tdir / "screen_scores.csv").exists():
            continue
        df = analysis.add_scores(pd.read_csv(tdir / "screen_scores.csv"))
        df = df.assign(name=df["name"].astype(str))
        sets = analysis.evaluation_sets(df)
        ifp = pipeline.run_interactions(cfg, t, df, sets, log=log, all_compounds=True)
        counts[t] = int(len(ifp))
    return counts


def box_diagnostic(cfg: Config, target_name: str, template: str, edges=(16.0, 20.0, 24.0, 30.0),
                   log=print) -> pd.DataFrame:
    """Item E: redock one template's own ligand across box sizes, to separate a box artefact from a
    genuinely undockable pose. Uses no activity label, so the choice it informs stays label-free."""
    from rdkit import Chem

    from .box import Box
    res = Path(cfg.results_dir)
    target = data.load_target(Path(cfg.data_root), target_name)
    kw = dict(seed=cfg.seed, nrun=cfg.autodock["nrun"], exhaustiveness=cfg.gnina["exhaustiveness"],
              num_modes=cfg.gnina["num_modes"], cpu=cfg.gnina["cpu"], metals_cfg=cfg.metals)
    rows = []
    for edge in edges:
        work = res / "box_diagnostic" / f"{target_name}_{template}_{edge:g}"
        work.mkdir(parents=True, exist_ok=True)
        xtal, _ = redock.crystal_ligand(target.ligand_mol2(template), work)
        centre = Chem.RemoveHs(xtal).GetConformer().GetPositions().mean(axis=0)
        box = Box(center=tuple(float(c) for c in centre), size=(float(edge),) * 3)
        r = redock.redock_template(target.protein_mol2(template), target.ligand_mol2(template), template,
                                   work, box_override=box, pdb_id=template, **kw)
        worst = max(r.ad_lowest_rmsd, r.gn_top_rmsd) if r.ok else float("nan")
        rows.append({"target": target_name, "template": template, "box_edge": float(edge),
                     "ad_rmsd": r.ad_lowest_rmsd, "gn_rmsd": r.gn_top_rmsd,
                     "gn_best_rmsd": r.gn_best_rmsd, "ok": bool(r.ok), "error": r.error,
                     "recovered": bool(np.isfinite(worst) and worst <= 2.0)})
        log(f"[{target_name}/{template}] box {edge:g} A: AD {r.ad_lowest_rmsd:.2f} "
            f"gnina {r.gn_top_rmsd:.2f} (best {r.gn_best_rmsd:.2f})")
    out = pd.DataFrame(rows)
    out.to_csv(res / f"box_diagnostic_{target_name}_{template}.csv", index=False)
    return out
