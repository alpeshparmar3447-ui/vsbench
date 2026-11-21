"""Config-driven, resumable benchmark pipeline: templates -> screen -> analyze."""
from __future__ import annotations

import json
import time
from dataclasses import MISSING, asdict, dataclass, field, fields
from pathlib import Path

import pandas as pd
from rdkit import Chem

from . import analysis, data, engines, gates, ligprep, prioritize, redock
from .box import box_from_coords
from .redock import RedockResult


@dataclass
class Config:
    data_root: str
    results_dir: str
    targets: list[str]
    n_actives: int = 40
    n_random: int = 300
    per_active_matched: int = 8
    seed: int = 42
    templates: dict[str, list[str]] | None = None   # optional per-target restriction
    box: dict = field(default_factory=lambda: {"buffer": 6.0, "min_edge": 24.0})
    autodock: dict = field(default_factory=lambda: {"nrun": 20, "autostop": True, "heuristics": True})
    gnina: dict = field(default_factory=lambda: {"exhaustiveness": 8, "num_modes": 3, "cpu": 16})
    n_boot: int = 1000
    # Replicates for the comparisons read against a multiplicity-adjusted interval (score
    # variants, and docking against the best structure-free baseline). A far quantile needs
    # more draws per tail than the 1,000 the other tables use; see PROTOCOL.md Amendments 4-5.
    variant_n_boot: int = 5000
    n_jobs: int = 16
    seeds: list | None = None            # AutoDock-GPU seeds; None means [seed]. The score is the mean over seeds.
    gnina_seed_study: dict | None = None # e.g. {"n": 60, "seeds": [43, 44]}: gnina repeated on a subset
    decoy_per_active: int = 8            # DUD-E-protocol decoys per active (0 disables the arm)
    decoy_pool_per_target: int = 4000    # inactives sampled from each other target to build the decoy pool
    # Receptor-level assessment (structure QC, pocket detection and druggability, metal audit).
    # None uses receptors.DEFAULTS: fpocket and the controls on, DoGSite and P2Rank off.
    receptors: dict | None = None
    metals: dict = field(default_factory=lambda: dict(redock.DEFAULT_METALS))
    # Hit analysis and lead selection: ADMET profiling of every compound, co-crystal interaction recovery
    # for the actives and the top-k of each evaluation set, and the prioritization comparison. None disables.
    lead_selection: dict | None = field(default_factory=lambda: {
        "top_k": 100, "top_n": [10, 20, 50], "admet": True, "interactions": True, "n_jobs": 8})

    @classmethod
    def from_json(cls, path: Path) -> "Config":
        """Read a run configuration, resolving relative paths against the config file itself.

        Relative ``data_root``/``results_dir`` are taken relative to the directory holding the
        config, not to the current working directory, so a checkout can be run from anywhere and
        the shipped configs need no absolute paths. Absolute values are left alone.
        """
        path = Path(path)
        try:
            raw = json.loads(path.read_text())
        except json.JSONDecodeError as exc:
            raise ValueError(f"{path}: not valid JSON ({exc})") from exc
        if not isinstance(raw, dict):
            raise ValueError(f"{path}: expected a JSON object, got {type(raw).__name__}")
        known = {f.name for f in fields(cls)}
        unknown = sorted(set(raw) - known)
        if unknown:
            raise ValueError(f"{path}: unknown configuration key(s) {unknown}; "
                             f"known keys are {sorted(known)}")
        missing = sorted({f.name for f in fields(cls)
                          if f.default is MISSING and f.default_factory is MISSING} - set(raw))
        if missing:
            raise ValueError(f"{path}: missing required configuration key(s) {missing}")
        for key in ("data_root", "results_dir"):
            raw[key] = str((path.parent / raw[key]).resolve()) if not Path(raw[key]).is_absolute() \
                else str(raw[key])
        return cls(**raw)

    def target_dir(self, target: str) -> Path:
        return Path(self.results_dir) / target


class Logger:
    def __init__(self, path: Path | None = None):
        self.path = path
        if path:
            path.parent.mkdir(parents=True, exist_ok=True)

    def __call__(self, msg: str) -> None:
        line = f"{time.strftime('%Y-%m-%d %H:%M:%S')}  {msg}"
        print(line, flush=True)
        if self.path:
            with open(self.path, "a") as fh:
                fh.write(line + "\n")


# ------------------------------------------------------------------ templates
def run_templates(cfg: Config, target_name: str, log=print) -> tuple[list[RedockResult], RedockResult]:
    """Redock every template's co-crystal ligand and select one without using labels."""
    target = data.load_target(Path(cfg.data_root), target_name)
    wanted = (cfg.templates or {}).get(target_name) or list(target.templates)
    out = cfg.target_dir(target_name) / "redock"
    results = []
    for t in wanted:
        cache = out / t / "result.json"
        if cache.exists():
            results.append(RedockResult(**json.loads(cache.read_text())))
            continue
        log(f"[{target_name}] redocking template {t}")
        r = redock.redock_template(target.protein_mol2(t), target.ligand_mol2(t), t, out / t,
                                   seed=cfg.seed, nrun=cfg.autodock["nrun"],
                                   exhaustiveness=cfg.gnina["exhaustiveness"],
                                   num_modes=cfg.gnina["num_modes"], cpu=cfg.gnina["cpu"],
                                   box_kwargs=cfg.box, metals_cfg=cfg.metals, pdb_id=t)
        (out / t).mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(r.as_dict()))
        results.append(r)
        log(f"[{target_name}]   {t}: ok={r.ok} AD {r.ad_lowest_rmsd:.2f}  gnina {r.gn_top_rmsd:.2f}  {r.error[:80]}")
    table = pd.DataFrame([r.as_dict() for r in results])
    table["success"] = [redock.is_success(r) for r in results]
    table.to_csv(cfg.target_dir(target_name) / "redock.csv", index=False)
    chosen = redock.select_template(results)
    (cfg.target_dir(target_name) / "selected.json").write_text(json.dumps(
        {"template": chosen.template, "worst_engine_rmsd": redock.worst_engine_rmsd(chosen),
         "validated": redock.is_success(chosen), "n_templates": len(results),
         "n_validated": int(table["success"].sum())}))
    return results, chosen


# --------------------------------------------------------------------- screen
def build_compound_table(cfg: Config, target) -> tuple[pd.DataFrame, data.NegativeSets]:
    pool = (data.build_decoy_pool(Path(cfg.data_root), target.name, cfg.decoy_pool_per_target, cfg.seed)
            if cfg.decoy_per_active else None)
    sets = data.draw_negative_sets(
        data.sample_rows(target.actives, cfg.n_actives, cfg.seed), target.inactives,
        n_random=cfg.n_random, per_active=cfg.per_active_matched, seed=cfg.seed, n_jobs=cfg.n_jobs,
        decoy_pool=pool, decoy_per_active=cfg.decoy_per_active)
    groups = [(sets.actives, 1, "active"), (sets.random, 0, "random"), (sets.matched, 0, "matched")]
    if sets.decoy is not None:
        groups.append((sets.decoy, 0, "decoy"))
    parts = []
    for frame, label, tag in groups:
        f = frame.copy()
        f["label"], f["neg_set"] = label, tag
        parts.append(f)
    table = pd.concat(parts, ignore_index=True)
    from . import similarity
    act_smiles = list(table.loc[table["neg_set"] == "active", "smiles"])
    act_fps = similarity.morgan_fps(act_smiles)
    all_fps = similarity.morgan_fps(list(table["smiles"]))
    nn = similarity.max_similarity(all_fps, act_fps)
    is_act = (table["neg_set"] == "active").to_numpy()
    nn[is_act] = similarity.max_similarity([f for f, a in zip(all_fps, is_act) if a], act_fps, exclude_self=True)
    table["nn_sim"] = nn                      # similarity to the actives (actives scored leave-one-out)
    if table["id"].duplicated().any():
        raise RuntimeError("a compound appears in more than one set")
    return table.rename(columns={"id": "name"}), sets


def run_screen(cfg: Config, target_name: str, template: str, log=print) -> pd.DataFrame:
    target = data.load_target(Path(cfg.data_root), target_name)
    tdir = cfg.target_dir(target_name)
    xtal, _ = redock.crystal_ligand(target.ligand_mol2(template), tdir / "redock" / template)
    heavy_xyz = Chem.RemoveHs(xtal).GetConformer().GetPositions()
    box = box_from_coords(heavy_xyz, **cfg.box)
    plan = redock.plan_for(target.protein_mol2(template), heavy_xyz, box, template, cfg.metals)
    log(f"[{target_name}] metals: {plan.note}")
    rec = engines.prepare_receptor(target.protein_mol2(template), tdir / "redock" / template / "receptor", box,
                                   keep_ions=plan.ions)

    table, sets = build_compound_table(cfg, target)
    log(f"[{target_name}] {len(table)} compounds ({(table.label == 1).sum()} actives); "
        f"{sets.n_unparsable} unparsable dropped")
    t0 = time.monotonic()
    ligs = ligprep.prepare_many(list(zip(table["smiles"], table["name"])), cfg.seed, cfg.n_jobs,
                                chelate_zn=plan.zn)
    prep_s = time.monotonic() - t0
    table["prep_ok"] = [l.ok for l in ligs]
    table["n_heavy"] = [l.n_heavy for l in ligs]
    table["prep_error"] = [l.error for l in ligs]
    table["chelators_auto"] = [l.chelators_auto for l in ligs]
    table["chelators_manual"] = [l.chelators_manual for l in ligs]
    table["deprotonated"] = [l.deprotonated for l in ligs]
    log(f"[{target_name}] ligand prep {prep_s:.1f}s, {int((~table.prep_ok).sum())} failed")
    table = table.merge(pd.DataFrame(gates.apply_gates(ligs)), on="name")

    seeds = list(cfg.seeds or [cfg.seed])
    ad_runs, ad_infos = {}, {}
    for i, sd in enumerate(seeds):
        wd = tdir / "screen" / ("adgpu" if i == 0 else f"adgpu_s{sd}")       # first seed keeps the historical path
        ad_runs[sd], ad_infos[sd] = engines.run_autodock_gpu(rec, ligs, wd, seed=sd, log=log, **cfg.autodock)
    ad, ad_info = merge_seed_runs(ad_runs), ad_infos[seeds[0]]
    gn, gn_info = engines.run_gnina(rec, ligs, tdir / "screen" / "gnina", seed=cfg.seed, log=log, **cfg.gnina)
    out = table.merge(ad, on="name").merge(gn, on="name")
    study = cfg.gnina_seed_study or {}
    if study.get("n"):
        keep = seed_study_subset(table, int(study["n"]), cfg.seed, study.get("min_negatives"))
        if table.loc[table["name"].isin(keep), "label"].nunique() < 2:
            log(f"  WARNING: the gnina seed-study subset holds one class only "
                f"({len(keep)} compounds); its ROC-AUC statistics cannot be computed")
        sub = [l for l in ligs if l.name in keep]
        for sd in study["seeds"]:
            g2, _ = engines.run_gnina(rec, sub, tdir / "screen" / f"gnina_s{sd}", seed=sd, log=log, **cfg.gnina)
            out = out.merge(g2[["name", "gn_cnn_vs", "gn_vina"]].rename(
                columns={"gn_cnn_vs": f"gn_cnn_vs_s{sd}", "gn_vina": f"gn_vina_s{sd}"}), on="name", how="left")
    out.to_csv(tdir / "screen_scores.csv", index=False)
    info_path = tdir / "screen_info.json"
    prev = json.loads(info_path.read_text()) if info_path.exists() else {}
    for key, info in (("autodock", ad_info), ("gnina", gn_info)):
        before = prev.get(key, {})
        # A resumed run docks nothing new; accumulate so first-run timings are never lost.
        info["wall_s_cumulative"] = before.get("wall_s_cumulative", 0.0) + info["wall_s"]
        info["n_docked_cumulative"] = before.get("n_docked_cumulative", 0) + info["n_docked_now"]
    info_path.write_text(json.dumps(
        {"template": template, "box": asdict(box), "prep_s": prev.get("prep_s", prep_s),
         "metals": {"ions": plan.label, "source": plan.source, "zn_mode": plan.zn, "audited": plan.audited,
                    "remote": list(plan.remote), "note": plan.note},
         "autodock": ad_info, "gnina": gn_info, "seeds": seeds,
         "autodock_seed_wall_s": {str(sd): ad_infos[sd]["wall_s"] for sd in seeds}, "n": len(out)}))
    return out


def seed_study_subset(table: pd.DataFrame, n: int, seed: int, min_negatives: int | None = None) -> set:
    """The compounds gnina is re-docked on to measure its own seed variance.

    Both classes must survive, because the study reports a per-seed ROC-AUC. Taking every active
    first and filling the remainder with negatives does that only while the actives fit inside the
    budget: a target with more actives than `n` left no room at all, so the subset was pure actives,
    the AUC was undefined and `gnina_seed_study` returned nothing. The registered 40-active arm never
    reached that point, but the 1000-active arm did on four of five targets, silently. The active
    count is therefore capped at `n - min_negatives` and the rest of the budget goes to negatives,
    which leaves a subset at or below `n` carrying both classes whenever the target has both.
    """
    min_neg = int(min_negatives if min_negatives is not None else max(1, n // 3))
    acts = table.loc[table.label == 1, "name"]
    negs = table.loc[table.label == 0, "name"]
    n_act = min(len(acts), max(1, n - min_neg))
    n_neg = min(len(negs), max(min_neg, n - n_act))
    kept_acts = acts if n_act >= len(acts) else acts.sample(n_act, random_state=seed)
    kept_negs = negs if n_neg >= len(negs) else negs.sample(n_neg, random_state=seed)
    return set(kept_acts) | set(kept_negs)


def merge_seed_runs(frames: dict) -> pd.DataFrame:
    """Combine per-seed AutoDock-GPU tables: per-seed columns plus the mean (the primary score) and SD."""
    import numpy as np
    seeds = list(frames)
    out = frames[seeds[0]].copy()
    names = out["name"]
    cols = {}
    for c in ("ad_energy", "ad_largest_energy", "ad_ok"):
        cols[c] = np.column_stack([frames[sd].set_index("name")[c].reindex(names).to_numpy(float) for sd in seeds])
        for j, sd in enumerate(seeds):
            out[f"{c}_s{sd}"] = cols[c][:, j]
    with np.errstate(all="ignore"):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            out["ad_energy"] = np.nanmean(cols["ad_energy"], axis=1)
            out["ad_largest_energy"] = np.nanmean(cols["ad_largest_energy"], axis=1)
            out["ad_energy_sd"] = np.nanstd(cols["ad_energy"], axis=1, ddof=1) if len(seeds) > 1 else np.nan
    out["ad_ok"] = np.nanmax(cols["ad_ok"], axis=1) > 0
    return out


# -------------------------------------------------------------------- analyze
def run_analysis(cfg: Config, target_name: str, df: pd.DataFrame | None = None, log=print) -> dict:
    tdir = cfg.target_dir(target_name)
    out = tdir / "analysis"
    out.mkdir(parents=True, exist_ok=True)
    if df is None:
        df = pd.read_csv(tdir / "screen_scores.csv")
    df = analysis.add_scores(df)
    sets = analysis.evaluation_sets(df)

    metrics, baselines, basecmp, pose_rows, engine_rows, balance = [], [], [], [], [], {}
    for name, sub in sets.items():
        metrics.append(analysis.evaluate(sub, name, cfg.n_boot, cfg.seed))
        baselines.append(analysis.property_baselines(sub, name, cfg.n_boot, cfg.seed))
        basecmp.append(analysis.baseline_comparison(sub, name, cfg.variant_n_boot, cfg.seed))
        pose_rows.append({"negative_set": name, **analysis.pose_selection_effect(sub, cfg.n_boot, cfg.seed)})
        engine_rows.append({"negative_set": name, **analysis.engine_difference(sub, cfg.n_boot, cfg.seed)})
        balance[name] = analysis.balance_table(sub[sub.label == 1], sub[sub.label == 0])
    result = {
        "metrics": pd.concat(metrics, ignore_index=True),
        "baselines": pd.concat(baselines, ignore_index=True),
        "pose_selection": pd.DataFrame(pose_rows),
        "engine_difference": pd.DataFrame(engine_rows),
        "balance": pd.DataFrame(balance),
        "failures": analysis.failure_table(df),
    }
    basecmp = [b for b in basecmp if len(b)]
    if basecmp:
        result["baseline_comparison"] = pd.concat(basecmp, ignore_index=True)
    import re
    from . import similarity
    seeds = sorted(int(c.rsplit("_s", 1)[1]) for c in df.columns if re.fullmatch(r"ad_energy_s\d+", c))
    if len(seeds) > 1:
        result["seed_stability"] = pd.DataFrame(
            [{"negative_set": n, **analysis.seed_stability(sub, seeds)} for n, sub in sets.items()])
    gseeds = sorted(int(c.rsplit("_s", 1)[1]) for c in df.columns if re.fullmatch(r"gn_cnn_vs_s\d+", c))
    if gseeds:
        gs = analysis.gnina_seed_study(df, gseeds)
        if gs:
            result["gnina_seed_study"] = pd.DataFrame([gs])
    diffs = []
    for a, b in (("random", "matched"), ("matched", "decoy"), ("random", "decoy")):
        if a in sets and b in sets:
            for label in analysis.PRIMARY:
                diffs.append({"set_a": a, "set_b": b, "score": label,
                              **analysis.negative_set_difference(sets[a], sets[b], analysis.SCORE_COLUMNS[label],
                                                                 cfg.n_boot, cfg.seed)})
    if diffs:
        result["negative_set_differences"] = pd.DataFrame(diffs)
    if cfg.lead_selection:
        result.update(run_lead_selection(cfg, target_name, df, sets, log))
    act_smi = list(df.loc[df["neg_set"] == "active", "smiles"])
    sim = similarity.audit(act_smi, {n: list(df.loc[df["neg_set"] == n, "smiles"]) for n in sets})
    (out / "similarity.json").write_text(json.dumps(sim, indent=2, default=float))
    for k, v in result.items():
        v.to_csv(out / f"{k}.csv", index=k == "balance")
    df.to_csv(out / "scored.csv", index=False)
    return result


# ------------------------------------------------------- hit analysis and lead selection
def run_admet(cfg: Config, target_name: str, df: pd.DataFrame, log=print) -> pd.DataFrame:
    """ADMET profile of every distinct compound (cached), one row per compound name."""
    from . import admet
    path = cfg.target_dir(target_name) / "admet.csv"
    names = df["name"].astype(str)          # compound ids may be numeric; always compare them as text
    if path.exists():
        cached = pd.read_csv(path, dtype={"name": str})
        if set(names) <= set(cached["name"]):
            return cached
    t0 = time.monotonic()
    prof = admet.profile(list(df["smiles"]))
    prof.insert(0, "name", list(names))
    prof.drop(columns=["smiles"]).to_csv(path, index=False)
    log(f"[{target_name}] ADMET profile of {len(prof)} compounds in {time.monotonic() - t0:.0f}s")
    return pd.read_csv(path, dtype={"name": str})


def gnina_top_poses(workdir: Path, names: set) -> dict:
    """Top-ranked gnina pose of each wanted compound, read from the chunked output SDFs."""
    wanted = {str(n) for n in names}
    out = {}
    for f in sorted(Path(workdir).glob("gnina_*.out.sdf")):
        for name, poses in engines.read_gnina_sdf(f).items():
            if name in wanted and poses:
                out[name] = poses[0]
    return out


def run_interactions(cfg: Config, target_name: str, df: pd.DataFrame, sets: dict, log=print,
                     all_compounds: bool = False) -> pd.DataFrame:
    """Recovery of the co-crystal interactions by the top gnina pose of every active and of the top-k
    compounds of each evaluation set (the pool the lead-selection comparison draws from).

    ``all_compounds`` extends that to every screened compound, which is what ranking a whole evaluation
    set by interaction recovery needs (PROTOCOL.md Amendment 3, item D). Results are cached per compound,
    so the wider pass only computes the poses the narrower one has not already done.
    """
    from . import interactions, sites
    tdir = cfg.target_dir(target_name)
    sel = json.loads((tdir / "selected.json").read_text())["template"]
    k = cfg.lead_selection["top_k"]
    wanted = set(df.loc[df["neg_set"] == "active", "name"].astype(str))
    if all_compounds:
        wanted = set(df["name"].astype(str))
    for sub in sets.values():
        wanted |= set(prioritize.select_pool(sub, k)["name"].astype(str))
    path = tdir / "interactions.csv"
    cached = pd.read_csv(path, dtype={"name": str}) if path.exists() else pd.DataFrame(columns=["name"])
    todo = wanted - set(cached["name"])
    if todo:
        t0 = time.monotonic()
        target = data.load_target(Path(cfg.data_root), target_name)
        work = tdir / "interactions"
        work.mkdir(exist_ok=True)
        prot = work / "protein.pdb"
        sites.write_pdb_from_mol2(target.protein_mol2(sel), prot)
        xtal, _ = redock.crystal_ligand(target.ligand_mol2(sel), tdir / "redock" / sel)
        ref = interactions.reference_set(prot, xtal)
        (work / "reference.json").write_text(json.dumps(sorted(map(list, ref)), default=str))
        poses = gnina_top_poses(tdir / "screen" / "gnina", todo)
        missing = pd.DataFrame({"name": sorted(todo - set(poses)), "ok": False, "error": "no docked pose"})
        # Written in batches: PLIP on a large receptor runs for over an hour, and a single write at the
        # end means an interruption loses all of it. Each batch lands in the cache, so a re-run resumes.
        names = sorted(poses)
        batch = max(1, int(cfg.lead_selection.get("interaction_batch", 100)))
        for i in range(0, len(names), batch):
            part = {n: poses[n] for n in names[i:i + batch]}
            new = interactions.score_poses(prot, ref, part, n_jobs=cfg.lead_selection.get("n_jobs", 8))
            cached = pd.concat([cached, new], ignore_index=True)
            cached["n_reference"] = len(ref)
            cached.to_csv(path, index=False)
            log(f"[{target_name}] interaction recovery {min(i + batch, len(names))}/{len(names)} "
                f"({time.monotonic() - t0:.0f}s)")
        if len(missing):
            cached = pd.concat([cached, missing], ignore_index=True)
        cached["n_reference"] = len(ref)
        cached.to_csv(path, index=False)
        log(f"[{target_name}] interaction recovery for {len(todo)} poses vs {len(ref)} crystal contacts "
            f"({time.monotonic() - t0:.0f}s)")
    return cached


def run_lead_selection(cfg: Config, target_name: str, df: pd.DataFrame, sets: dict, log=print) -> dict:
    """ADMET, interaction recovery and the prioritization comparison for one target."""
    ls = cfg.lead_selection
    tables = {}
    df = df.assign(name=df["name"].astype(str))
    merged = df
    if ls.get("admet", True):
        adm = run_admet(cfg, target_name, df, log)
        keep = [c for c in adm.columns if c not in merged.columns or c == "name"]
        merged = merged.merge(adm[keep], on="name", how="left")
        tables["admet"] = adm
    if ls.get("interactions", True):
        ifp = run_interactions(cfg, target_name, df, {k: v.assign(name=v["name"].astype(str)) for k, v in sets.items()}, log)
        sub = ifp.rename(columns={"recall": "ifp_recall", "tanimoto": "ifp_tanimoto",
                                  "residue_recall": "ifp_residue_recall", "ok": "ifp_ok"})
        merged = merged.merge(sub[["name", "ifp_recall", "ifp_tanimoto", "ifp_residue_recall", "ifp_ok"]],
                              on="name", how="left")
        tables["interactions"] = ifp
    msets = analysis.evaluation_sets(merged)
    if "ifp_recall" in merged:
        for s in msets.values():
            # a pose PLIP could not analyse, or no pose at all, reproduces nothing: counted, not dropped
            s["ifp_recall"] = s["ifp_recall"].fillna(0.0) if s["ifp_recall"].notna().any() else s["ifp_recall"]
    tables.update(prioritize.run(msets, ls["top_k"], tuple(ls.get("top_n", (10, 20, 50))), cfg.n_boot, cfg.seed))
    return {k: v for k, v in tables.items() if k in ("prioritization", "admet_filters")} | {
        "interactions_pool": tables.get("interactions", pd.DataFrame())}


def aggregate(cfg: Config) -> pd.DataFrame:
    """Stack every target's metrics and add the macro-average over targets.

    The panel files describe everything analysed under the results directory, not only the targets
    this invocation happened to run. Keying them to ``cfg.targets`` meant that re-running one target
    after an interruption - which ``--targets`` sets to that single name - rewrote
    ``metrics_all_targets.csv`` from it alone. ``report.py`` reads that file as its metrics source, so
    the next report would have described a one-target panel without saying so. The per-target
    ``analysis/metrics.csv`` files are the record and these two are derived from them, so rebuilding
    from every target present on disk is both the correct scope and self-healing.
    """
    root = Path(cfg.results_dir)
    targets = sorted(set(cfg.targets) | {p.parent.parent.name for p in root.glob("*/analysis/metrics.csv")})
    frames = []
    for t in targets:
        p = cfg.target_dir(t) / "analysis" / "metrics.csv"
        if p.exists():
            frames.append(pd.read_csv(p).assign(target=t))
    if not frames:
        raise FileNotFoundError("no per-target metrics found; run the analyze stage first")
    allm = pd.concat(frames, ignore_index=True)
    macro = (allm.groupby(["score", "negative_set", "metric"], as_index=False)
             .agg(value=("value", "mean"), sd_across_targets=("value", "std"),
                  n_targets=("target", "nunique")))
    macro["target"] = "MACRO"
    allm.to_csv(Path(cfg.results_dir) / "metrics_all_targets.csv", index=False)
    macro.to_csv(Path(cfg.results_dir) / "metrics_macro.csv", index=False)
    return macro


def run_all(cfg: Config, stages=("templates", "screen", "analyze"), log=None) -> list[str]:
    """Run the stages for every target. A target that fails is recorded and skipped, so one
    bad receptor cannot abort a multi-hour campaign. Returns the list of failed targets."""
    import traceback
    Path(cfg.results_dir).mkdir(parents=True, exist_ok=True)
    log = log or Logger(Path(cfg.results_dir) / "run.log")
    (Path(cfg.results_dir) / "config.json").write_text(json.dumps(asdict(cfg), indent=2))
    failed = []
    for t in cfg.targets:
        try:
            chosen = None
            if "templates" in stages:
                _, chosen = run_templates(cfg, t, log)
            if "screen" in stages:
                if chosen is None:
                    sel = json.loads((cfg.target_dir(t) / "selected.json").read_text())
                    chosen = RedockResult(template=sel["template"], ok=True)
                run_screen(cfg, t, chosen.template, log)
            if "analyze" in stages:
                run_analysis(cfg, t, log=log)
                log(f"[{t}] analysis written")
        except Exception as exc:
            failed.append(t)
            cfg.target_dir(t).mkdir(parents=True, exist_ok=True)
            (cfg.target_dir(t) / "error.txt").write_text(traceback.format_exc())
            log(f"[{t}] FAILED: {type(exc).__name__}: {exc}")
    if "analyze" in stages and len(failed) < len(cfg.targets):
        aggregate(Config(**{**asdict(cfg), "targets": [t for t in cfg.targets if t not in failed]}))
    if failed:
        log(f"targets that failed: {failed}")
    return failed
