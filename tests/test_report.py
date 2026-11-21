import numpy as np
import pandas as pd

from vsbench import report as R


def test_to_markdown_shape_and_blank_for_nan():
    df = pd.DataFrame({"a": [1, 2], "b": ["x", np.nan], "c": [0.5, float("nan")]})
    lines = R.to_markdown(df).splitlines()
    assert lines[0] == "| a | b | c |" and lines[1] == "|---|---|---|"
    assert lines[2] == "| 1 | x | 0.5 |"
    assert lines[3].startswith("| 2 |  |")       # NaN renders blank, never the string 'nan'
    assert "nan" not in R.to_markdown(df)


def _metrics():
    rows = []
    for t, v in (("T1", 0.61), ("T2", 0.55)):
        for neg in ("random", "matched"):
            for eng in ("autodock-gpu", "gnina", "consensus"):
                rows.append(dict(target=t, negative_set=neg, score=eng, metric="roc_auc",
                                 value=v, ci_low=v - 0.05, ci_high=v + 0.05))
                rows.append(dict(target=t, negative_set=neg, score=eng, metric="pr_auc",
                                 value=0.1, ci_low=0.05, ci_high=0.2))
    return pd.DataFrame(rows)


def test_auc_table_formats_intervals_and_orders_engines():
    t = R.auc_table(_metrics(), "matched")
    assert list(t.columns) == ["target", "autodock-gpu", "gnina", "consensus"]
    assert t.loc[t.target == "T1", "gnina"].item() == "0.61 [0.56, 0.66]"
    assert len(t) == 2


def test_baseline_table_pivots_by_set_and_property():
    rows = [dict(target=t, negative_set=n, property=p, roc_auc=0.5 + 0.1 * (n == "random"))
            for t in ("T1", "T2") for n in ("random", "matched") for p in ("n_heavy", "mw", "tpsa", "hba")]
    t = R.baseline_table(pd.DataFrame(rows))
    assert "random: n_heavy" in t.columns and "matched: tpsa" in t.columns
    assert "random: hba" not in t.columns
    assert t.loc[0, "random: n_heavy"] == 0.6 and t.loc[0, "matched: n_heavy"] == 0.5


def test_metal_audit_table_summarises_pocket_and_remote_ions():
    audit = pd.DataFrame([
        dict(target="FEN1", template="5fv7", frame_overlap=0.999, ions_in_benchmark_file=0, ions_in_deposition=2,
             relevant_ions=2, element="MG", ion="MG1A", dist_to_ligand_A=2.29, nearest_functional_ion_within_30A="MG@2.3"),
        dict(target="TP53", template="1abc", frame_overlap=0.996, ions_in_benchmark_file=0, ions_in_deposition=1,
             relevant_ions=0, element=np.nan, ion=np.nan, dist_to_ligand_A=np.nan, nearest_functional_ion_within_30A="ZN@29.0"),
        dict(target="TP53", template="2abc", frame_overlap=1.0, ions_in_benchmark_file=0, ions_in_deposition=0,
             relevant_ions=0, element=np.nan, ion=np.nan, dist_to_ligand_A=np.nan, nearest_functional_ion_within_30A=""),
        dict(target="ESR1_ant", template="1xp1", frame_overlap=1.0, ions_in_benchmark_file=0, ions_in_deposition=0,
             relevant_ions=0, element=np.nan, ion=np.nan, dist_to_ligand_A=np.nan, nearest_functional_ion_within_30A=""),
    ])
    t = R.metal_audit_table(audit).set_index("target")
    assert t.loc["FEN1", "templates with a pocket ion"] == 1 and t.loc["FEN1", "pocket ion (nearest)"] == "MG 2.3 A"
    assert t.loc["TP53", "templates"] == 2 and t.loc["TP53", "templates with a pocket ion"] == 0
    assert t.loc["TP53", "nearest functional ion within 30 A"] == "ZN 29.0 A"
    assert t.loc["TP53", "min frame overlap"] == "0.996"
    assert t.loc["ESR1_ant", "pocket ion (nearest)"] == "" and t.loc["ESR1_ant", "nearest functional ion within 30 A"] == ""


# ------------------------------------------------------------- tables for the extended pipeline
def _target_dir(tmp_path, t="T1"):
    d = tmp_path / t
    (d / "analysis").mkdir(parents=True)
    return d


def test_negset_table_counts_targets_whose_interval_excludes_zero(tmp_path):
    for t, (lo, hi, d) in {"A": (0.02, 0.10, 0.06), "B": (-0.03, 0.05, 0.01)}.items():
        pd.DataFrame([{"set_a": "random", "set_b": "matched", "score": "gnina", "delta_auc": d,
                       "ci_low": lo, "ci_high": hi, "p_two_sided": 0.1}]).to_csv(
            _target_dir(tmp_path, t) / "analysis" / "negative_set_differences.csv", index=False)
    out = R.negset_table(tmp_path, ["A", "B"])
    assert out.loc[0, "mean_delta_auc"] == "+0.035" and out.loc[0, "targets_with_ci_excluding_zero"] == 1
    assert R.negset_table(tmp_path, ["MISSING"]) is None


def test_admet_and_interaction_tables_join_on_text_ids(tmp_path):
    d = _target_dir(tmp_path)
    pd.DataFrame({"name": [1, 2, 3, 4], "label": [1, 1, 0, 0]}).to_csv(d / "screen_scores.csv", index=False)
    pd.DataFrame({"name": ["1", "2", "3", "4"], "tier": ["PASS", "EXCLUDE", "PASS", "PASS"]}).to_csv(d / "admet.csv", index=False)
    pd.DataFrame({"name": ["1", "2", "3", "4"], "ok": [True, True, True, False], "recall": [0.8, 0.6, 0.2, np.nan],
                  "n_reference": 9}).to_csv(d / "interactions.csv", index=False)
    a = R.admet_class_table(tmp_path, ["T1"]).iloc[0]
    assert a["actives PASS"] == "50%" and a["actives EXCLUDE"] == "50%" and a["negatives PASS"] == "100%"
    i = R.interaction_table(tmp_path, ["T1"]).iloc[0]
    assert i["recall, actives"] == "0.70" and i["recall, top-ranked negatives"] == "0.20"
    assert i["n, top-ranked negatives"] == 1 and i["crystal contacts"] == 9          # the failed pose is not averaged in


def test_prioritization_and_filter_macro_tables(tmp_path):
    for t, hr in (("A", 0.4), ("B", 0.2)):
        d = _target_dir(tmp_path, t)
        pd.DataFrame([{"negative_set": "matched", "strategy": "docking", "top_n": 10, "hit_rate": hr,
                       "pool_active_fraction": 0.3}]).to_csv(d / "analysis" / "prioritization.csv", index=False)
        pd.DataFrame([{"negative_set": "matched", "filter": "Veber pass", "active_retained": 0.9, "inactive_retained": 0.8,
                       "actives_lost": 2, "roc_auc_before": 0.6, "roc_auc_after": 0.62}]).to_csv(
            d / "analysis" / "admet_filters.csv", index=False)
    p = R.prioritization_table(tmp_path, ["A", "B"]).iloc[0]
    assert p["hit_rate"] == "0.300" and p["targets"] == 2
    f = R.filter_macro_table(tmp_path, ["A", "B"]).iloc[0]
    assert f["actives_lost"] == 4 and f["actives_retained"] == "90%" and f["auc_after"] == "0.620"


def test_site_tables_use_wilson_intervals_and_keep_controls_separate():
    sites = pd.DataFrame({"method": ["fpocket"] * 4 + ["dogsite"] * 2, "top1_hit": [True, True, False, False, True, True],
                          "top3_hit": [True] * 4 + [True, True], "best_recovery": [0.8] * 6, "n_pockets": [10] * 6})
    t = R.site_summary_table(sites).set_index("detector")
    assert t.loc["fpocket", "top-1 recovery"].startswith("2/4 (50%) [15%, 85%]")
    assert t.loc["dogsite", "top-3 recovery"].startswith("2/2 (100%)")
    c = R.site_control_table(pd.DataFrame([{"template": "1HXB", "control_role": "positive", "description": "x",
                                           "method": "fpocket", "rank_of_site": 4, "site_pocket_druggability": 0.055,
                                           "ligand_buried_fraction": 0.5}]))
    assert c.loc[0, "site_pocket_druggability"] == "0.06" and c.loc[0, "ligand_buried_fraction"] == "0.50"


def test_structure_qc_table_counts_flags():
    qc = pd.DataFrame({"target": ["A", "A", "B"], "template": ["1", "2", "3"], "resolution": [2.0, 3.0, 1.5],
                       "n_unexplained_breaks": [0, 1, 0], "flags": ["multiple_chains", "", None],
                       "engineered_mutation": [0, 2, 1], "rama_core": [0.9, 0.8, 0.95]})
    t = R.structure_qc_table(qc).set_index("target")
    assert t.loc["A", "templates"] == 2 and t.loc["A", "with_chain_break"] == 1 and t.loc["A", "multiple_chains"] == 1
    assert t.loc["B", "engineered_mutation"] == 1 and t.loc["A", "median_resolution"] == 2.5


def test_table_numbering_is_gapless_and_split_by_placement():
    """Numbers follow position, so an absent optional table cannot leave a gap."""
    d = pd.DataFrame({"x": [1]})
    spec = [(R.MAIN, "a", d), (R.SUPP, "b", d), (R.MAIN, "c", d), (R.SUPP, "e", d), (R.SUPP, "f", d)]
    got = [title for _, title, _ in R.number_tables(spec)]
    assert got == ["Table 1. a", "Table S1. b", "Table 2. c", "Table S2. e", "Table S3. f"]
    # Dropping an entry renumbers rather than leaving a hole: the Table-28 gap cannot recur.
    kept = [e for e in spec if e[1] != "c"]
    assert [t for _, t, _ in R.number_tables(kept)] == [
        "Table 1. a", "Table S1. b", "Table S2. e", "Table S3. f"]


def test_merging_across_negative_sets_keeps_every_row():
    """One analysis per negative set is one table; no row may be lost in the stacking."""
    parts = {"random": pd.DataFrame({"score": ["a", "b"], "v": [1, 2]}),
             "matched": pd.DataFrame({"score": ["a"], "v": [3]}),
             "decoy": None}
    out = R._stack_by_negative_set(parts.get)
    assert list(out.columns) == ["negative set", "score", "v"]
    assert len(out) == 3
    assert set(out["negative set"]) == {"random inactives", "property-matched inactives"}
    assert R._stack_by_negative_set(lambda ns: None) is None


def test_size_control_table_lists_each_score_once(tmp_path):
    """"gnina alone" is the registered "gnina": the same column must not appear twice."""
    pd.DataFrame([{"score": s, "negative_set": "matched", "roc_auc": 0.6,
                   "conditional_auc": 0.59, "mean_shift": -0.01, "sd_shift": 0.03, "n_targets": 5}
                  for s in ("gnina", "gnina alone", "autodock-gpu", "AutoDock-GPU alone",
                            "fusion: rank mean")]).to_csv(
        tmp_path / "size_control_macro.csv", index=False)
    t = R.size_control_table(tmp_path, "matched")
    assert sorted(t["score"]) == ["autodock-gpu", "fusion: rank mean", "gnina"]


def test_single_target_standard_deviation_renders_blank_not_nan(tmp_path):
    d = _target_dir(tmp_path)
    pd.DataFrame([{"negative_set": "matched", "strategy": "docking", "top_n": 10, "hit_rate": 0.4,
                   "pool_active_fraction": 0.3}]).to_csv(d / "analysis" / "prioritization.csv", index=False)
    t = R.prioritization_table(tmp_path, ["T1"])
    assert t.loc[0, "sd_across_targets"] == "" and "nan" not in R.to_markdown(t)


# ------------------------------------------------- Amendment 3 tables
def _write_amendment3(tmp_path):
    pd.DataFrame([
        {"score": "consensus", "negative_set": "matched", "roc_auc": 0.541,
         "conditional_auc": 0.540, "mean_shift": -0.001, "sd_shift": 0.023, "n_targets": 5},
        {"score": "autodock-gpu", "negative_set": "matched", "roc_auc": 0.529,
         "conditional_auc": 0.517, "mean_shift": -0.012, "sd_shift": 0.030, "n_targets": 5},
    ]).to_csv(tmp_path / "size_control_macro.csv", index=False)
    pd.DataFrame([
        {"target": "T1", "negative_set": "random", "median_gap": 5.5},
        {"target": "T1", "negative_set": "matched", "median_gap": 0.5},
        {"target": "T2", "negative_set": "random", "median_gap": -1.0},
        {"target": "T2", "negative_set": "matched", "median_gap": 0.0},
    ]).to_csv(tmp_path / "size_profile.csv", index=False)
    pd.DataFrame([
        {"score": "gnina (CNN score)", "negative_set": "decoy", "metric": "ef1", "value": 3.45,
         "sd_across_targets": 1.0, "mean_vs_primary": 1.95, "n_better": 0, "n_worse": 0, "n_targets": 5},
        {"score": "gnina alone", "negative_set": "decoy", "metric": "ef1", "value": 2.70,
         "sd_across_targets": 1.0, "mean_vs_primary": 1.20, "n_better": 1, "n_worse": 0, "n_targets": 5},
    ]).to_csv(tmp_path / "score_variants_macro.csv", index=False)
    pd.DataFrame([
        {"target": "MTORC1", "templates": 11, "ran": 10, "validated": 0, "macrocyclic": 9,
         "median_heavy_atoms": 65.0, "median_mol_weight": 914.2, "median_rotatable_bonds": 6.0,
         "best_rmsd": 6.2, "dockable_structures": True},
    ]).to_csv(tmp_path / "template_difficulty_per_target.csv", index=False)
    pd.DataFrame([
        {"property": "max_ring_size", "spearman_rho": 0.490, "p_permutation": 0.002, "n_templates": 42},
        {"property": "rotatable_bonds", "spearman_rho": 0.207, "p_permutation": 0.190, "n_templates": 42},
    ]).to_csv(tmp_path / "template_difficulty_correlations.csv", index=False)
    pd.DataFrame([
        {"group": "non-macrocyclic", "n_templates": 42, "validated": 16, "rate": 0.381,
         "ci_low": 0.25, "ci_high": 0.532},
        {"group": "macrocyclic", "n_templates": 10, "validated": 0, "rate": 0.0,
         "ci_low": 0.0, "ci_high": 0.278},
        {"group": "ligand not sanitizable", "n_templates": 10, "validated": 1, "rate": 0.1,
         "ci_low": float("nan"), "ci_high": float("nan")},
    ]).to_csv(tmp_path / "template_difficulty_by_macrocycle.csv", index=False)
    pd.DataFrame([
        {"target": "MTORC1", "template": "4jsx", "box_edge": 16.0, "ad_rmsd": 1.1, "gn_rmsd": 1.4,
         "gn_best_rmsd": 0.9, "ok": True, "error": "", "recovered": True},
        {"target": "MTORC1", "template": "4jsx", "box_edge": 24.0, "ad_rmsd": 6.13, "gn_rmsd": 6.21,
         "gn_best_rmsd": 5.79, "ok": True, "error": "", "recovered": False},
    ]).to_csv(tmp_path / "box_diagnostic_MTORC1_4jsx.csv", index=False)


def test_size_control_table_shows_the_change_and_sorts_by_the_controlled_value(tmp_path):
    _write_amendment3(tmp_path)
    t = R.size_control_table(tmp_path, "matched")
    assert list(t.score) == ["consensus", "autodock-gpu"]
    assert t.loc[1, "change"] == "-0.012" and t.loc[0, "change"] == "-0.001"
    assert t.loc[0, "size-controlled AUC"] == "0.540"
    assert R.size_control_table(tmp_path, "decoy") is None      # absent set, not an empty table


def test_size_gap_table_signs_the_gap(tmp_path):
    _write_amendment3(tmp_path)
    t = R.size_gap_table(tmp_path).set_index("target")
    assert t.loc["T1", "random"] == "+5.5" and t.loc["T2", "random"] == "-1.0"
    assert t.loc["T2", "matched"] == "+0.0"


def test_score_variant_table_keeps_ef1_precision_and_the_target_counts(tmp_path):
    _write_amendment3(tmp_path)
    t = R.score_variant_table(tmp_path, "decoy", "ef1")
    assert list(t.score) == ["gnina (CNN score)", "gnina alone"]
    assert t.loc[0, "ef1"] == "3.45" and t.loc[0, "vs registered consensus"] == "+1.95"
    assert t.loc[1, "targets better"] == 1 and t.loc[0, "targets"] == 5


def test_template_tables_render_the_mtorc1_finding(tmp_path):
    _write_amendment3(tmp_path)
    d = R.template_difficulty_table(tmp_path)
    assert d.loc[0, "macrocyclic ligands"] == 9 and d.loc[0, "median MW"] == "914"
    assert d.loc[0, "validated (<=2 A)"] == 0
    p = R.template_property_table(tmp_path)
    assert p.loc[0, "ligand property"] == "largest ring size"
    assert p.loc[0, "Spearman rho vs RMSD"] == "+0.490" and p.loc[0, "p (permutation)"] == "0.002"
    m = R.macrocycle_table(tmp_path)
    assert m.loc[1, "validation rate"] == "0.0%" and m.loc[1, "95% CI"] == "0.0% to 27.8%"
    assert m.loc[2, "95% CI"] == ""                              # unknown group has no interval


def test_box_diagnostic_table_marks_recovery(tmp_path):
    _write_amendment3(tmp_path)
    t = R.box_diagnostic_table(tmp_path)
    assert list(t["recovered"]) == ["yes", "no"]
    assert t.loc[0, "box edge (A)"] == "16" and t.loc[1, "gnina RMSD (A)"] == "6.21"


def test_amendment3_tables_are_absent_rather_than_empty_when_nothing_ran(tmp_path):
    for fn in (R.size_gap_table, R.template_difficulty_table, R.template_property_table,
               R.macrocycle_table, R.box_diagnostic_table):
        assert fn(tmp_path) is None
    assert R.score_variant_table(tmp_path, "decoy") is None
