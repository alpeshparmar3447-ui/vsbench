import numpy as np
import pandas as pd
import pytest

from vsbench import metrics as M
from vsbench import prioritize as P


def test_pareto_fronts_match_a_hand_worked_example():
    # two objectives, both maximised
    X = np.array([[3, 3],      # dominates everything below
                  [2, 2],      # dominated by 0 only          -> front 1
                  [3, 1],      # dominated by 0               -> front 1
                  [1, 1],      # dominated by 1, 2, 0         -> front 2
                  [4, 0]])     # nobody dominates it          -> front 0
    f = P.pareto_fronts(X)
    assert list(f) == [0, 1, 1, 2, 0]


def test_pareto_minimised_objective_is_flipped():
    maxi = np.array([[1.0], [1.0], [1.0]])
    mini = np.array([[0.0], [1.0], [2.0]])        # lower is better
    assert list(P.pareto_fronts(maxi, mini)) == [0, 1, 2]


def test_pareto_ties_share_a_front():
    f = P.pareto_fronts(np.array([[1.0, 1.0], [1.0, 1.0]]))
    assert list(f) == [0, 0]                       # equal points do not dominate each other


def test_zscore_handles_a_constant_column():
    assert np.allclose(P.zscore([2.0, 2.0, 2.0]), 0.0)
    z = P.zscore([1.0, 2.0, 3.0])
    assert z.mean() == pytest.approx(0.0) and z.std() == pytest.approx(1.0)


def _pool(n=60, seed=0):
    rng = np.random.default_rng(seed)
    y = np.r_[np.ones(n // 3, int), np.zeros(n - n // 3, int)]
    dock = rng.normal(size=n) + 0.8 * y
    ifp = np.clip(rng.normal(0.4, 0.15, n) + 0.3 * y, 0, 1)
    tier = rng.choice(["PASS", "MODERATE", "HIGH_RISK", "EXCLUDE"], size=n, p=[0.5, 0.2, 0.15, 0.15])
    risk = pd.Series(tier).map({"PASS": 0, "MODERATE": 1, "HIGH_RISK": 2, "EXCLUDE": 3}).to_numpy()
    return pd.DataFrame({"label": y, "consensus": dock, "ifp_recall": ifp, "tier": tier, "risk": risk})


def test_select_pool_uses_only_the_docking_score():
    df = _pool()
    pool = P.select_pool(df, 20)
    assert len(pool) == 20
    assert pool["consensus"].min() >= df["consensus"].sort_values(ascending=False).iloc[19] - 1e-12
    assert list(pool["consensus"]) == sorted(pool["consensus"], reverse=True)


def test_strategies_orderings():
    pool = P.select_pool(_pool(), 40)
    s = P.strategies(pool)
    assert set(s) == {"docking", "docking+ADMET", "docking+IFP", "pareto"}
    assert list(s["docking"]) == list(np.argsort(-pool["consensus"].to_numpy(), kind="mergesort"))
    admet = pool["tier"].to_numpy()[s["docking+ADMET"]]
    assert set(admet) <= {"PASS", "MODERATE"}                         # the filter removed the risky tiers
    assert sorted(s["pareto"]) == list(range(len(pool)))              # a permutation, nobody lost
    top = s["pareto"][0]
    assert P.pareto_fronts(np.column_stack([pool["consensus"], pool["ifp_recall"]]),
                           pool["risk"].to_numpy(float)[:, None])[top] == 0


def test_strategies_without_interaction_data_skip_ifp_and_pareto():
    pool = P.select_pool(_pool().drop(columns=["ifp_recall"]), 30)
    assert set(P.strategies(pool)) == {"docking", "docking+ADMET"}


def test_hit_rate_table_counts_and_interval():
    pool = P.select_pool(_pool(90), 60)
    t = P.hit_rate_table(pool, n_list=(5, 10), n_boot=200, seed=1)
    assert set(t.strategy) == {"docking", "docking+ADMET", "docking+IFP", "pareto"}
    row = t[(t.strategy == "docking") & (t.top_n == 10)].iloc[0]
    order = P.strategies(pool)["docking"][:10]
    assert row.actives_in_top_n == int(pool["label"].to_numpy()[order].sum())
    assert row.hit_rate == pytest.approx(row.actives_in_top_n / 10)
    assert row.ci_low <= row.hit_rate + 0.35 and row.ci_high >= row.hit_rate - 0.35
    assert (t.ci_low <= t.ci_high).all()
    assert row.pool_active_fraction == pytest.approx(pool["label"].mean())
    # reproducible
    t2 = P.hit_rate_table(pool, n_list=(5, 10), n_boot=200, seed=1)
    pd.testing.assert_frame_equal(t, t2)


def test_a_filter_that_removes_actives_lowers_the_hit_rate():
    pool = P.select_pool(_pool(90), 60)
    pool = pool.assign(tier=np.where(pool.label == 1, "EXCLUDE", "PASS"))   # the filter hates every active
    t = P.hit_rate_table(pool, n_list=(10,), n_boot=50)
    assert t[t.strategy == "docking+ADMET"].iloc[0].actives_in_top_n == 0


def test_filter_effect_agrees_with_independent_metrics():
    df = _pool(120, seed=3)
    keep = df["tier"].isin(["PASS", "MODERATE"])
    e = P.filter_effect(df, keep, "x")
    y, s, k = df["label"].to_numpy(), df["consensus"].to_numpy(), keep.to_numpy()
    assert e["active_retained"] == pytest.approx(k[y == 1].mean())
    assert e["inactive_retained"] == pytest.approx(k[y == 0].mean())
    assert e["actives_lost"] == int(((~k) & (y == 1)).sum())
    from sklearn.metrics import roc_auc_score
    assert e["roc_auc_before"] == pytest.approx(roc_auc_score(y, s))
    assert e["roc_auc_after"] == pytest.approx(roc_auc_score(y[k], s[k]))
    assert e["ef5_before"] == pytest.approx(M.enrichment_factor(y, s, 0.05))


def test_filter_effect_with_one_class_left_returns_nan_not_an_error():
    df = _pool(60)
    keep = pd.Series(df["label"] == 0)               # every active removed
    e = P.filter_effect(df, keep, "kill actives")
    assert e["actives_lost"] == int(df["label"].sum()) and np.isnan(e["roc_auc_after"])


def test_filter_table_has_one_row_per_rule():
    df = _pool(80).assign(veber_pass=lambda d: d.risk < 3, n_pains=lambda d: (d.risk == 3).astype(int))
    t = P.filter_table(df)
    assert len(t) == 4 and t["filter"].is_unique
    assert ((t.active_retained >= 0) & (t.active_retained <= 1)).all()
