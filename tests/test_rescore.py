import numpy as np
import pandas as pd
import pytest

from vsbench import rescore as R


def test_rank_fusion_is_invariant_to_a_monotone_rescaling_of_one_engine():
    """The registered z-score consensus can be dominated by whichever engine has the wider spread.
    Rank fusion cannot, and that is the reason for including it."""
    rng = np.random.default_rng(0)
    a, b = rng.normal(size=50), rng.normal(size=50)
    plain = R.rank_mean({"ad": a, "gn": b})
    blown = R.rank_mean({"ad": a * 1000 + 5, "gn": b})
    assert np.argsort(plain).tolist() == np.argsort(blown).tolist()


def test_zscore_fusion_is_not_invariant_which_is_what_rank_fusion_fixes():
    from vsbench.metrics import consensus_zscore
    rng = np.random.default_rng(1)
    a, b = rng.normal(size=40), rng.normal(size=40)
    z1 = consensus_zscore({"ad": a, "gn": b})
    z2 = consensus_zscore({"ad": a * 1000, "gn": b})
    assert np.allclose(np.argsort(z1), np.argsort(z2))          # z-scoring already removes pure scaling
    z3 = consensus_zscore({"ad": np.sign(a) * a ** 2, "gn": b})  # a non-linear stretch does change it
    assert not np.allclose(np.argsort(z1), np.argsort(z3))


def test_reciprocal_rank_rewards_agreement_at_the_top():
    both = {"ad": np.array([10.0, 9.0, 1.0]), "gn": np.array([10.0, 1.0, 9.0])}
    s = R.reciprocal_rank(both)
    assert s[0] > s[1] and s[0] > s[2]                          # ranked first by both engines
    assert s[1] == pytest.approx(s[2])                          # symmetric disagreement


def test_reciprocal_rank_rejects_a_non_positive_k():
    with pytest.raises(ValueError):
        R.reciprocal_rank({"a": np.array([1.0, 2.0])}, k=0)


def test_weights_track_redocking_quality_and_fall_back_to_equal():
    w = R.redocking_weights(0.5, 6.0)
    assert w["ad"] > w["gn"] and w["ad"] == pytest.approx(1 / 1.5)
    assert R.redocking_weights(np.nan, np.nan) == {"ad": 1.0, "gn": 1.0}
    assert R.redocking_weights(0.0, np.nan)["gn"] == 0.0


def test_weighted_fusion_follows_the_trusted_engine():
    good = np.array([3.0, 2.0, 1.0])
    bad = np.array([1.0, 2.0, 3.0])
    s = R.weighted_zscore({"ad": good, "gn": bad}, {"ad": 10.0, "gn": 0.0})
    assert s[0] > s[1] > s[2]
    with pytest.raises(ValueError):
        R.weighted_zscore({"ad": good}, {})
    with pytest.raises(ValueError):
        R.weighted_zscore({"ad": good}, {"ad": 0.0})


def test_a_constant_engine_contributes_no_spread_rather_than_dividing_by_zero():
    s = R.weighted_zscore({"ad": np.array([1.0, 2.0, 3.0]), "gn": np.array([5.0, 5.0, 5.0])},
                          {"ad": 1.0, "gn": 1.0})
    assert np.all(np.isfinite(s)) and s[2] > s[0]


def _frame(n=60, seed=0):
    rng = np.random.default_rng(seed)
    y = np.r_[np.ones(15, int), np.zeros(n - 15, int)]
    return pd.DataFrame({
        "label": y, "ad_lowest": rng.normal(size=n) + y * 1.2, "gn_cnn_vs": rng.normal(size=n) + y * 0.9,
        "gn_cnn_score": rng.normal(size=n), "gn_cnn_affinity": rng.normal(size=n) + y * 0.4,
        "consensus": rng.normal(size=n) + y * 1.0})


def test_build_variants_adds_every_column_and_needs_both_engines():
    d = R.build_variants(_frame())
    for c in ("fusion_rank_mean", "fusion_rrf", "fusion_weighted_z",
              "score_gnina_cnn_score", "score_gnina_cnn_affinity"):
        assert c in d.columns and np.isfinite(d[c]).all()
    assert "gnina (CNN affinity)" in R.variant_columns(d)
    with pytest.raises(ValueError):
        R.build_variants(_frame().drop(columns=["gn_cnn_vs"]))


def test_non_finite_recorded_scores_are_imputed_to_the_worst_value():
    f = _frame()
    f.loc[0, "gn_cnn_affinity"] = np.nan
    d = R.build_variants(f)
    col = d["score_gnina_cnn_affinity"]
    assert np.isfinite(col).all() and col.iloc[0] == col.min()


def test_compare_reports_paired_differences_against_the_primary(monkeypatch):
    d = R.build_variants(_frame())
    t = R.compare(d, "matched", n_boot=120, seed=0)
    assert set(t.metric) == {"roc_auc", "ef1", "bedroc20"}
    assert (t.ci_low <= t.value + 1e-9).all() and (t.value <= t.ci_high + 1e-9).all()
    assert (t.diff_ci_low <= t.vs_primary + 1e-9).all()
    assert ((t.p >= 0) & (t.p <= 1)).all()
    assert not (t["better"] & t["worse"]).any()
    with pytest.raises(ValueError):
        R.compare(d.drop(columns=["consensus"]), "matched", n_boot=10)


def test_a_variant_identical_to_the_primary_shows_no_difference():
    d = R.build_variants(_frame())
    d["consensus"] = d["ad_lowest"]
    t = R.compare(d, "random", n_boot=80, seed=0)
    row = t[(t.score == "AutoDock-GPU alone") & (t.metric == "roc_auc")].iloc[0]
    assert row.vs_primary == pytest.approx(0.0) and not row.better and not row.worse


def test_promotion_requires_both_hard_negative_sets():
    rows = []
    for s, sets in (("A", {"matched": True, "decoy": True}), ("B", {"matched": True, "decoy": False}),
                    ("C", {"matched": False, "decoy": False})):
        for neg, better in sets.items():
            rows.append({"score": s, "negative_set": neg, "metric": "roc_auc", "better": better})
    t = pd.DataFrame(rows)
    assert set(R.promoted(t).score) == {"A"}


def test_promotion_ignores_a_win_on_the_random_set_alone():
    t = pd.DataFrame([{"score": "A", "negative_set": "random", "metric": "roc_auc", "better": True}])
    assert len(R.promoted(t)) == 0


def test_size_normalisation_divides_by_ligand_size():
    s = R.size_normalise(np.array([10.0, 10.0]), np.array([10.0, 20.0]))
    assert s[0] == pytest.approx(1.0) and s[1] == pytest.approx(0.5)
    c = R.size_normalise(np.array([8.0]), np.array([8.0]), power=1 / 3)
    assert c[0] == pytest.approx(4.0)                            # 8 / 8**(1/3) = 8 / 2


def test_size_normalisation_removes_a_pure_size_signal():
    """A score that is exactly proportional to heavy atoms carries no information once divided by them."""
    n = np.arange(10.0, 30.0)
    s = R.size_normalise(n * 0.7, n)
    assert np.allclose(s, s[0])


def test_a_missing_or_zero_heavy_atom_count_cannot_make_the_score_infinite():
    s = R.size_normalise(np.array([5.0, 5.0, 5.0]), np.array([10.0, 0.0, np.nan]))
    assert np.isfinite(s).all() and s[1] == s.min() and s[2] == s.min()


def test_build_variants_adds_the_size_normalised_scores_when_sizes_are_present():
    f = _frame()
    f["n_heavy"] = np.linspace(15, 40, len(f))
    f["gn_vina_hb"] = np.linspace(1, 5, len(f))
    d = R.build_variants(f)
    cols = R.variant_columns(d)
    assert "AutoDock-GPU per heavy atom" in cols
    assert "gnina (Vina-type) per heavy atom^(1/3)" in cols
    for c in cols.values():
        assert np.isfinite(d[c]).all()
    assert len(set(cols.values())) == len(cols)                  # no two variants share a column


def test_size_normalised_scores_are_skipped_when_no_sizes_are_recorded():
    d = R.build_variants(_frame())                               # no n_heavy column
    assert not any("per heavy atom" in k for k in R.variant_columns(d))


def test_interaction_recall_is_compared_against_the_primary_when_present():
    """Item D ranks a whole evaluation set by how much of the crystal interaction pattern each pose
    reproduces. It is supplied by the caller, so variant_columns must pick it up when it is there."""
    d = R.build_variants(_frame())
    assert "interaction recall" not in R.variant_columns(d)
    d["score_interaction_recall"] = np.linspace(0.0, 1.0, len(d))
    cols = R.variant_columns(d)
    assert cols["interaction recall"] == "score_interaction_recall"
    t = R.compare(d, "matched", n_boot=80, seed=0)
    assert "interaction recall" in set(t.score)
    row = t[(t.score == "interaction recall") & (t.metric == "roc_auc")].iloc[0]
    assert np.isfinite(row.value) and np.isfinite(row.vs_primary)


# ------------------------------------------- multiplicity across the variants tried at once
def test_compare_reports_an_adjusted_interval_wider_than_the_pointwise_one():
    df = R.build_variants(_frame())
    t = R.compare(df, "matched", n_boot=200, seed=0)
    assert (t["n_comparisons"] == len(R.variant_columns(df))).all()
    assert (t["diff_adj_low"] <= t["diff_ci_low"]).all()
    assert (t["diff_adj_high"] >= t["diff_ci_high"]).all()
    # Widening can only remove a win, never create one.
    assert not (t["better_adjusted"] & ~t["better"]).any()
    assert not (t["worse_adjusted"] & ~t["worse"]).any()


def test_compare_records_whether_the_adjusted_interval_is_resolved():
    """At the configured 1000 draws over a dozen variants it is not, and the table says so."""
    t = R.compare(R.build_variants(_frame()), "matched", n_boot=200, seed=0)
    assert not t["adj_resolved"].any()


def test_promotion_uses_the_adjusted_flag_when_the_table_has_one():
    rows = []
    for neg in ("matched", "decoy"):
        rows.append({"score": "A", "negative_set": neg, "metric": "roc_auc",
                     "better": True, "better_adjusted": False})
        rows.append({"score": "B", "negative_set": neg, "metric": "roc_auc",
                     "better": True, "better_adjusted": True})
    # A wins on every pointwise interval and would have been promoted before the adjustment.
    assert set(R.promoted(pd.DataFrame(rows)).score) == {"B"}


def _studies_cfg(tmp_path, **kw):
    from vsbench.pipeline import Config
    (tmp_path / "T").mkdir(exist_ok=True)
    (tmp_path / "T" / "screen_scores.csv").write_text("name\nc1\n")
    return Config(data_root=str(tmp_path), results_dir=str(tmp_path), targets=["T"], **kw)


def _stub_studies(monkeypatch, captured):
    """Everything rescore_study touches except the draw-count resolution under test."""
    from vsbench import analysis, rescore, sizecontrol, studies, templatestats
    monkeypatch.setattr(templatestats, "template_table",
                        lambda *a, **k: pd.DataFrame({"macrocycle": [False], "validated": [True]}))
    for fn in ("correlations", "by_macrocycle", "per_target"):
        monkeypatch.setattr(templatestats, fn, lambda *a, **k: pd.DataFrame({"x": [1]}))
    monkeypatch.setattr(studies, "_scored_with_variants",
                        lambda cfg, t: pd.DataFrame({"name": ["c1"], "label": [1]}))
    monkeypatch.setattr(rescore, "variant_columns",
                        lambda d: {f"v{i}": f"col{i}" for i in range(12)})
    monkeypatch.setattr(analysis, "evaluation_sets", lambda d: {"matched": d})

    def fake_compare(sub, name, n_boot=1000, seed=0, **k):
        captured.append(n_boot)
        return pd.DataFrame({"score": ["v0"], "negative_set": [name]})

    monkeypatch.setattr(rescore, "compare", fake_compare)
    monkeypatch.setattr(sizecontrol, "table", lambda *a, **k: (
        pd.DataFrame({"score": ["s"], "negative_set": ["matched"], "roc_auc": [0.5],
                      "conditional_auc": [0.5], "shift": [0.0]}), pd.DataFrame({"x": [1]})))
    monkeypatch.setattr(sizecontrol, "size_profile", lambda *a, **k: {"n": 1})
    monkeypatch.setattr(studies, "variant_macro", lambda allv: pd.DataFrame({"x": [1]}))
    monkeypatch.setattr(rescore, "promoted", lambda allv: pd.DataFrame({"score": []}))


def test_the_variant_comparison_takes_its_draw_count_from_the_config(tmp_path, monkeypatch):
    """`variant_n_boot` used to fall back to `n_boot` rather than to the config, so a run that did
    not pass --variant-n-boot compared the variants at 1000 draws while its own manifest recorded
    5000. Nothing failed; the adjusted interval was simply unresolved."""
    from vsbench import studies
    captured: list[int] = []
    _stub_studies(monkeypatch, captured)
    studies.rescore_study(_studies_cfg(tmp_path, n_boot=1000, variant_n_boot=5000), log=lambda *a: None)
    assert captured == [5000]


def test_an_explicit_variant_draw_count_still_overrides_the_config(tmp_path, monkeypatch):
    from vsbench import studies
    captured: list[int] = []
    _stub_studies(monkeypatch, captured)
    studies.rescore_study(_studies_cfg(tmp_path, n_boot=1000, variant_n_boot=5000),
                          variant_n_boot=250, log=lambda *a: None)
    assert captured == [250]


def test_the_configured_defaults_are_why_the_fallback_mattered():
    """The adjusted interval reads a far quantile: the config default resolves it for the reported
    twelve-variant family and the pointwise default does not, which is the whole point of the
    separate setting."""
    from vsbench import metrics as M
    from vsbench.pipeline import Config
    d = Config(data_root=".", results_dir=".", targets=[])
    assert M.draws_resolve_interval(d.variant_n_boot, n_comparisons=12)
    assert not M.draws_resolve_interval(d.n_boot, n_comparisons=12)
