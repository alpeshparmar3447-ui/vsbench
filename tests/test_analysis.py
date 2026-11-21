import numpy as np
import pandas as pd
import pytest

from vsbench import analysis as A
from vsbench import metrics as M


def _table(n_act=40, n_inact=800, seed=0, signal=0.0, size_effect=0.0):
    """Synthetic screening table. ``signal`` adds real discrimination to the engine
    scores; ``size_effect`` makes scores depend on heavy-atom count only, while
    actives are larger than inactives (a pure property confound)."""
    rng = np.random.default_rng(seed)
    y = np.r_[np.ones(n_act), np.zeros(n_inact)].astype(int)
    n = len(y)
    heavy = np.where(y == 1, rng.normal(32, 5, n), rng.normal(26, 5, n)).round()
    mw = heavy * 13.5 + rng.normal(0, 8, n)
    df = pd.DataFrame({
        "name": [str(i) for i in range(n)], "label": y, "prep_ok": True,
        "n_heavy": heavy, "mw": mw, "clogp": rng.normal(3, 1, n) + 0.02 * heavy,
        "hbd": rng.integers(0, 4, n), "hba": rng.integers(2, 9, n), "rotb": rng.integers(1, 9, n),
        "tpsa": rng.normal(80, 20, n), "charge": 0,
    })
    base = rng.normal(size=n)
    df["ad_energy"] = -(8 + base + signal * y + size_effect * (heavy - 26) / 5)
    df["ad_largest_energy"] = df["ad_energy"] + rng.normal(0, 0.1, n)
    df["ad_ok"] = True
    df["gn_vina"] = -(7 + rng.normal(size=n) + signal * y)
    df["gn_cnn_vs"] = 3 + rng.normal(size=n) + signal * y + size_effect * (heavy - 26) / 5
    df["gn_ok"] = True
    return df


def test_impute_worst_and_scores_are_higher_is_better():
    v = A.impute_worst([1.0, np.nan, 5.0])
    assert v.tolist() == [1.0, 1.0, 5.0]
    with pytest.raises(ValueError):
        A.impute_worst([np.nan, np.nan])
    df = A.add_scores(_table(signal=1.0))
    assert M.roc_auc(df["label"], df["ad_lowest"]) > 0.6         # more negative energy = better
    assert M.roc_auc(df["label"], df["consensus"]) > 0.6


def test_a_failed_docking_counts_against_the_compound():
    df = _table(signal=2.0)
    df.loc[df.index[0], ["ad_energy", "ad_ok"]] = [np.nan, False]   # a failed ACTIVE
    out = A.add_scores(df)
    assert out.loc[out.index[0], "ad_lowest"] == out["ad_lowest"].min()
    ft = A.failure_table(df)
    assert ft.loc[ft["class"] == "active", "autodock_failed"].item() == 1
    assert ft.loc[ft["class"] == "inactive", "autodock_failed"].item() == 0


def _sets_from_pool(size_effect=0.0, signal=0.0, seed=0, per_active=6):
    """Build screening tables the way the pipeline does: draw random and matched
    negatives from a large inactive pool, then attach engine scores."""
    from vsbench import data as D
    rng = np.random.default_rng(seed)
    n_act, n_pool = 40, 4000
    act = pd.DataFrame({"smiles": ["C" * int(k) + "O" for k in rng.integers(16, 40, n_act)],
                        "id": [f"a{i}" for i in range(n_act)]})
    ina = pd.DataFrame({"smiles": ["C" * int(k) + "O" for k in rng.integers(2, 40, n_pool)],
                        "id": [f"i{i}" for i in range(n_pool)]})
    sets = D.draw_negative_sets(act, ina, n_random=240, per_active=per_active, seed=seed)

    def frame(d, label, tag):
        out = d.copy()
        out["label"], out["neg_set"], out["prep_ok"] = label, tag, True
        out["n_heavy"] = [__import__("rdkit.Chem", fromlist=["x"]).MolFromSmiles(s).GetNumHeavyAtoms()
                          for s in out.smiles]
        r = np.random.default_rng(seed + 1)
        eff = size_effect * (out["n_heavy"] - 20) / 5
        out["ad_energy"] = -(8 + r.normal(size=len(out)) + signal * label + eff)
        out["ad_largest_energy"] = out["ad_energy"] + r.normal(0, 0.1, len(out))
        out["gn_vina"] = -(7 + r.normal(size=len(out)) + signal * label)
        out["gn_cnn_vs"] = 3 + r.normal(size=len(out)) + signal * label + eff
        out["ad_ok"] = out["gn_ok"] = True
        return out

    return pd.concat([frame(sets.actives, 1, "active"), frame(sets.random, 0, "random"),
                      frame(sets.matched, 0, "matched")], ignore_index=True)


def test_matching_removes_a_pure_size_confound_that_inflates_random_negatives():
    # Engine "scores" depend on size only. Against random inactives this looks like
    # skill; against property-matched inactives it must fall to near chance.
    df = A.add_scores(_sets_from_pool(size_effect=3.0, seed=4))
    sets = A.evaluation_sets(df)
    auc = {k: M.roc_auc(v["label"], v["gn_cnn_vs"]) for k, v in sets.items()}
    assert auc["random"] > 0.70
    assert auc["matched"] < auc["random"] - 0.15
    assert abs(auc["matched"] - 0.5) < 0.12
    base_r = A.property_baselines(sets["random"], "random", n_boot=100).set_index("property")
    base_m = A.property_baselines(sets["matched"], "matched", n_boot=100).set_index("property")
    assert base_r.loc["n_heavy", "roc_auc"] > 0.65      # the confound is visible without docking
    assert abs(base_m.loc["n_heavy", "roc_auc"] - 0.5) < 0.12


def test_matching_preserves_real_signal():
    df = A.add_scores(_sets_from_pool(signal=1.5, seed=7))
    sets = A.evaluation_sets(df)
    assert M.roc_auc(sets["matched"]["label"], sets["matched"]["gn_cnn_vs"]) > 0.75
    assert M.roc_auc(sets["random"]["label"], sets["random"]["gn_cnn_vs"]) > 0.75


def test_evaluation_sets_share_the_same_actives_and_keep_negatives_apart():
    df = A.add_scores(_sets_from_pool(seed=2))
    sets = A.evaluation_sets(df)
    ra, ma = sets["random"], sets["matched"]
    assert (ra.label == 1).sum() == (ma.label == 1).sum() == 40
    assert set(ra[ra.label == 0].id).isdisjoint(ma[ma.label == 0].id)
    assert set(ra[ra.label == 0].neg_set) == {"random"} and set(ma[ma.label == 0].neg_set) == {"matched"}


def test_evaluate_reports_every_metric_and_matches_direct_calls():
    df = A.add_scores(_table(signal=1.0, seed=2))
    ev = A.evaluate(df, "random", n_boot=100, seed=0, columns={"autodock-gpu": "ad_lowest"})
    assert set(ev.metric) == set(M.METRICS)
    y, s = df["label"].to_numpy(), df["ad_lowest"].to_numpy()
    row = ev[ev.metric == "roc_auc"].iloc[0]
    assert row.value == pytest.approx(M.roc_auc(y, s))
    assert row.ci_low < row.value < row.ci_high
    assert row.n_active == 40 and row.n_inactive == 800


def test_pose_selection_effect_is_zero_when_poses_agree_and_detects_a_real_gap():
    df = A.add_scores(_table(signal=1.0, seed=3))
    same = df.assign(ad_largest=df["ad_lowest"])
    r = A.pose_selection_effect(same, n_boot=100)
    assert r["delta_auc_largest_minus_lowest"] == 0.0 and r["p_two_sided"] == 1.0
    rng = np.random.default_rng(0)
    worse = df.assign(ad_largest=rng.normal(size=len(df)))
    r2 = A.pose_selection_effect(worse, n_boot=200)
    assert r2["delta_auc_largest_minus_lowest"] < -0.1 and r2["p_two_sided"] < 0.05


def test_engine_difference_detects_the_better_engine():
    df = A.add_scores(_table(seed=6))
    df["gn_cnn_vs"] = df["gn_cnn_vs"] + 2.5 * df["label"]          # make gnina clearly better
    r = A.engine_difference(df, n_boot=200)
    assert r["delta_auc_gnina_minus_autodock"] > 0.2 and r["ci_low"] > 0


def test_negative_set_difference_detects_an_easy_negative_set_and_is_null_for_identical_sets():
    rng = np.random.default_rng(0)
    n_act, n_neg = 40, 300
    names = [f"a{i}" for i in range(n_act)]
    act = pd.DataFrame({"name": names, "label": 1, "s": rng.normal(1.0, 1.0, n_act)})
    hard = pd.DataFrame({"name": [f"h{i}" for i in range(n_neg)], "label": 0, "s": rng.normal(0.8, 1.0, n_neg)})
    easy = pd.DataFrame({"name": [f"e{i}" for i in range(n_neg)], "label": 0, "s": rng.normal(-1.0, 1.0, n_neg)})
    hard_set, easy_set = pd.concat([act, hard]), pd.concat([act, easy])
    r = A.negative_set_difference(easy_set, hard_set, "s", n_boot=300, seed=1)
    assert r["delta_auc"] > 0.3 and r["ci_low"] > 0 and r["p_two_sided"] < 0.05
    same = A.negative_set_difference(hard_set, hard_set, "s", n_boot=300, seed=1)
    assert same["delta_auc"] == 0.0 and same["ci_low"] <= 0 <= same["ci_high"]


def test_evaluation_sets_include_the_decoy_arm_only_when_present():
    base = _sets_from_pool(seed=1)
    assert set(A.evaluation_sets(A.add_scores(base))) == {"random", "matched"}
    decoy = base[base.neg_set == "random"].copy()
    decoy["neg_set"] = "decoy"
    decoy["name"] = decoy["id"] + "_d"
    both = pd.concat([base, decoy], ignore_index=True)
    assert set(A.evaluation_sets(A.add_scores(both))) == {"random", "matched", "decoy"}


def test_property_baselines_use_the_analogue_column_when_present():
    df = A.add_scores(_sets_from_pool(seed=2))
    df["nn_sim"] = np.where(df.label == 1, 0.8, 0.2)
    b = A.property_baselines(A.evaluation_sets(df)["random"], "random", n_boot=50).set_index("property")
    assert b.loc["nn_sim", "roc_auc"] == 1.0
