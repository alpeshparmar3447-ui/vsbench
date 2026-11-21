import numpy as np
import pandas as pd
import pytest

from vsbench import data as D
from vsbench import similarity as Sim

BENZ = ["c1ccccc1C", "c1ccccc1CC", "c1ccccc1CCC", "c1ccccc1C(C)C", "c1ccccc1CCCC", "c1ccccc1OC", "c1ccccc1CO"]
ALK = ["CCCCCCCC", "CCCCCCCCC", "CCCCCCCCCC", "CCCCCCCCCCC", "C1CCCCCCC1", "CC(C)CCCCC"]
OTHER = ["c1ccncc1CC", "C1CCNCC1CC", "CC(=O)OC1CCCCC1", "OC1CCCCC1", "CN1CCCC1", "O=C(O)CCC(=O)O"]


def test_identical_molecules_have_similarity_one_and_unrelated_ones_less():
    a = Sim.morgan_fps(["c1ccccc1CCO", "c1ccccc1CCO", "CCCCCCCCCC"])
    s = Sim.max_similarity(a, a, exclude_self=True)
    assert s[0] == pytest.approx(1.0) and s[2] < 0.3                       # the twin is found, the alkane has no analogue
    assert np.isnan(Sim.max_similarity(Sim.morgan_fps(["not a smiles"]), a)[0])


def test_exclude_self_does_not_match_a_compound_with_itself():
    a = Sim.morgan_fps(["c1ccccc1CCO", "CCCCCCCCCC"])
    assert Sim.max_similarity(a, a, exclude_self=True)[0] < 1.0


def test_butina_clusters_separate_two_chemical_series():
    sizes = Sim.butina_cluster_sizes(Sim.morgan_fps(BENZ + ALK), similarity_cutoff=0.25)
    assert len(sizes) >= 2 and sizes[0] >= 3 and sum(sizes) == len(BENZ) + len(ALK)
    assert Sim.butina_cluster_sizes([]) == []


def test_scaffold_sharing():
    assert Sim.scaffold_sharing(BENZ, ["c1ccccc1CCN", "c1ccccc1Cl"]) == 1.0      # same benzene scaffold
    assert Sim.scaffold_sharing(BENZ, ["c1ccncc1CC", "C1CCNCC1CC"]) == 0.0
    assert np.isnan(Sim.scaffold_sharing(BENZ, ["CCCC", "CCO"]))                  # acyclic: no scaffold to compare


def test_analogue_baseline_detects_analogue_bias_and_not_its_absence():
    point, lo, hi = Sim.analogue_baseline_auc(BENZ, ALK, n_boot=200)
    assert point > 0.8 and lo > 0.6                                                # actives are analogues, negatives are not
    rng = np.random.default_rng(0)
    pool = (BENZ + ALK + OTHER) * 2
    rng.shuffle(pool)
    flat, _, _ = Sim.analogue_baseline_auc(pool[:12], pool[12:], n_boot=200)
    assert 0.2 < flat < 0.8


def test_audit_reports_every_negative_set():
    rep = Sim.audit(BENZ, {"random": ALK, "analogues": ["c1ccccc1CCN", "c1ccccc1CCCl"]})
    assert rep["actives"]["n"] == 7 and rep["actives"]["median_nn_similarity"] > 0.2
    assert rep["analogues"]["median_max_similarity"] > rep["random"]["median_max_similarity"]
    assert rep["analogues"]["scaffold_sharing"] == 1.0 and rep["random"]["analogue_baseline_auc"] > 0.8


# ----------------------------------------------------------------- DUD-E-protocol decoys
def _pool_and_actives():
    rng = np.random.default_rng(3)
    act_smi = ["c1ccccc1" + "C" * k + "N" for k in range(2, 12)]                        # one benzylamine-like series
    pool = [("C" * k + "O", f"p{i}") for i, k in enumerate(rng.integers(2, 30, 1500))]   # alkanols
    pool += [("c1ccccc1" + "C" * k + "N", f"s{i}") for i, k in enumerate(rng.integers(2, 12, 400))]   # close analogues of the actives
    return (pd.DataFrame({"smiles": act_smi, "id": [f"a{i}" for i in range(len(act_smi))]}),
            pd.DataFrame(pool, columns=["smiles", "id"]).drop_duplicates("smiles").reset_index(drop=True))


def test_dude_decoys_are_property_matched_disjoint_and_less_similar_than_property_matching_alone():
    act, pool = _pool_and_actives()
    adesc = D.compute_descriptors(act["smiles"])
    d = D.draw_dude_decoys(adesc, act["smiles"], pool, per_active=3, seed=1, candidates=40)
    assert len(d) == 30 and d["id"].is_unique
    sim_dude = np.nanmean(Sim.max_similarity(Sim.morgan_fps(d["smiles"]), Sim.morgan_fps(act["smiles"])))
    idx = D.match_inactives(adesc, D.compute_descriptors(pool["smiles"]), per_active=3, seed=1)
    sim_prop = np.nanmean(Sim.max_similarity(Sim.morgan_fps(pool.iloc[idx]["smiles"]), Sim.morgan_fps(act["smiles"])))
    assert sim_dude < sim_prop                                                       # the dissimilarity filter removed analogues
    assert D.draw_dude_decoys(adesc, act["smiles"], pool, 3, 1, candidates=40)["id"].tolist() == d["id"].tolist()


def test_decoy_pool_excludes_everything_tested_against_the_target(tmp_path):
    for name, act, ina in (("X", ["CCO 1"], ["CCN 2", "CCC 3"]), ("Y", ["CCCO 10"], ["CCN 2", "CCCC 11", "CCCCC 12"]),
                           ("Z", ["CO 20"], ["CCCCCC 13", "CCCCCCC 3"])):
        (tmp_path / name).mkdir()
        (tmp_path / name / "actives.smi").write_text("\n".join(act) + "\n")
        (tmp_path / name / "inactives.smi").write_text("\n".join(ina) + "\n")
    pool = D.build_decoy_pool(tmp_path, "X", per_target=10, seed=0)
    assert set(pool["id"]) == {"11", "12", "13"}                                      # ids 2 and 3 were tested against X
    assert not set(pool["id"]) & {"1", "2", "3"}


def test_decoy_arm_validates_pool_size():
    act, pool = _pool_and_actives()
    with pytest.raises(ValueError, match="too small"):
        D.draw_dude_decoys(D.compute_descriptors(act["smiles"]), act["smiles"], pool.head(15), per_active=3, seed=0)
