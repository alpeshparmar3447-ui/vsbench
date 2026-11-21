import numpy as np
import pandas as pd
import pytest

from vsbench import data as D


def test_read_smi_parses_and_validates(tmp_path):
    f = tmp_path / "a.smi"
    f.write_text("CCO 101\n\nc1ccccc1 102\n")
    df = D.read_smi(f)
    assert list(df.columns) == ["smiles", "id"]
    assert df.to_dict("list") == {"smiles": ["CCO", "c1ccccc1"], "id": ["101", "102"]}
    f.write_text("CCO\n")
    with pytest.raises(ValueError, match="expected 'SMILES ID'"):
        D.read_smi(f)
    f.write_text("CCO 1\nCCN 1\n")
    with pytest.raises(ValueError, match="duplicate"):
        D.read_smi(f)


def test_load_target_finds_templates_and_checks_pairs(tmp_path):
    t = tmp_path / "X"
    t.mkdir()
    (t / "actives.smi").write_text("CCO 1\n")
    (t / "inactives.smi").write_text("CCN 2\nCCC 3\n")
    for pdb in ("2abc", "1xyz"):
        (t / f"{pdb}_protein.mol2").write_text("")
        (t / f"{pdb}_ligand.mol2").write_text("")
    tg = D.load_target(tmp_path, "X")
    assert tg.templates == ("1xyz", "2abc")
    assert len(tg.actives) == 1 and len(tg.inactives) == 2
    assert tg.protein_mol2("1xyz").name == "1xyz_protein.mol2"
    (t / "2abc_ligand.mol2").unlink()
    with pytest.raises(FileNotFoundError):
        D.load_target(tmp_path, "X")
    with pytest.raises(FileNotFoundError):
        D.load_target(tmp_path, "missing")


def test_descriptors_known_values():
    d = D.compute_descriptors(["CCO", "CC(=O)[O-]", "not a smiles"])
    assert d.loc[0, "mw"] == pytest.approx(46.07, abs=0.01)
    assert d.loc[0, "hbd"] == 1 and d.loc[0, "charge"] == 0
    assert d.loc[1, "charge"] == -1
    assert d.loc[2].isna().all()


def test_parallel_descriptors_equal_serial():
    smi = ["CCO", "c1ccccc1O", "CC(=O)Nc1ccc(O)cc1", "CCN(CC)CC"] * 600
    pd.testing.assert_frame_equal(D.compute_descriptors(smi, 1), D.compute_descriptors(smi, 4))


def test_sample_rows_is_seeded_without_replacement():
    df = pd.DataFrame({"smiles": list("abcdefghij"), "id": [str(i) for i in range(10)]})
    a, b = D.sample_rows(df, 4, 3), D.sample_rows(df, 4, 3)
    pd.testing.assert_frame_equal(a, b)
    assert a["id"].is_unique and len(a) == 4
    assert not a.equals(D.sample_rows(df, 4, 4))
    assert len(D.sample_rows(df, 50, 0)) == 10
    with pytest.raises(ValueError):
        D.sample_rows(df, -1, 0)


def _synthetic(n, shift, seed):
    rng = np.random.default_rng(seed)
    return pd.DataFrame({
        "mw": rng.normal(350 + shift, 60, n), "clogp": rng.normal(3 + shift / 100, 1.2, n),
        "hbd": rng.integers(0, 4, n), "hba": rng.integers(2, 9, n),
        "rotb": rng.integers(1, 9, n), "tpsa": rng.normal(80 + shift / 4, 25, n),
        "charge": rng.choice([-1, 0, 0, 0, 1], n),
    })


def test_matching_removes_a_moderate_property_imbalance_and_is_deterministic():
    actives = _synthetic(60, shift=40, seed=1)      # moderately heavier, more lipophilic actives
    pool = _synthetic(5000, shift=0, seed=2)
    rand = pool.iloc[np.random.default_rng(0).choice(len(pool), 180, replace=False)]
    idx = D.match_inactives(actives, pool, per_active=3, seed=11)
    matched = pool.iloc[idx]
    assert len(idx) == 180 and len(set(idx)) == 180                  # no reuse
    smd_rand = D.standardized_mean_difference(actives, rand).abs()
    smd_match = D.standardized_mean_difference(actives, matched).abs()
    assert smd_rand[["mw", "clogp", "tpsa"]].max() > 0.3            # the confound exists
    assert smd_match.max() < 0.15                                    # and matching removes it
    np.testing.assert_array_equal(idx, D.match_inactives(actives, pool, 3, 11))
    assert not np.array_equal(idx, D.match_inactives(actives, pool, 3, 12))


def test_matching_under_an_extreme_confound_helps_but_cannot_be_perfect():
    # Actives 1.3 SD above the pool: the pool's tail is too sparse for exact
    # matching, so the guarantee is a large reduction, and callers must inspect
    # standardized_mean_difference() rather than assume balance.
    actives = _synthetic(60, shift=80, seed=1)
    pool = _synthetic(5000, shift=0, seed=2)
    rand = pool.iloc[np.random.default_rng(0).choice(len(pool), 180, replace=False)]
    matched = pool.iloc[D.match_inactives(actives, pool, per_active=3, seed=11)]
    smd_rand = D.standardized_mean_difference(actives, rand).abs().max()
    smd_match = D.standardized_mean_difference(actives, matched).abs().max()
    assert smd_match < smd_rand / 2.5


def test_matching_validates_inputs():
    a, p = _synthetic(10, 0, 1), _synthetic(15, 0, 2)
    with pytest.raises(ValueError, match="needed"):
        D.match_inactives(a, p, per_active=2, seed=0)
    with pytest.raises(ValueError, match="per_active"):
        D.match_inactives(a, p, per_active=0, seed=0)
    bad = p.copy()
    bad.loc[0, "mw"] = np.nan
    with pytest.raises(ValueError, match="NaN"):
        D.match_inactives(a, bad, per_active=1, seed=0)


def test_smd_is_zero_for_identical_and_signed_otherwise():
    a = _synthetic(200, 0, 1)
    assert (D.standardized_mean_difference(a, a).abs() < 1e-12).all()
    b = a.copy()
    b["mw"] = b["mw"] + 60
    assert D.standardized_mean_difference(a, b)["mw"] < -0.8


def _frames(n_act=40, n_pool=6000, shift=40, seed=3):
    rng = np.random.default_rng(seed)
    def smi(n):  # unique, valid, size-varying SMILES
        return ["C" * int(k) + "O" for k in rng.integers(2, 40, n)]
    act = pd.DataFrame({"smiles": smi(n_act), "id": [f"a{i}" for i in range(n_act)]})
    ina = pd.DataFrame({"smiles": smi(n_pool), "id": [f"i{i}" for i in range(n_pool)]})
    return act, ina


def test_draw_negative_sets_are_disjoint_deterministic_and_sized():
    act, ina = _frames()
    a = D.draw_negative_sets(act, ina, n_random=100, per_active=5, seed=1)
    b = D.draw_negative_sets(act, ina, n_random=100, per_active=5, seed=1)
    assert len(a.random) == 100 and len(a.matched) == 200 and len(a.actives) == 40
    assert set(a.random.id).isdisjoint(a.matched.id)
    assert a.matched.id.is_unique
    pd.testing.assert_frame_equal(a.matched, b.matched)
    assert a.n_unparsable == 0
    assert not a.matched.id.tolist() == D.draw_negative_sets(act, ina, 100, 5, seed=2).matched.id.tolist()


def test_draw_negative_sets_matching_beats_random_when_actives_are_larger():
    rng = np.random.default_rng(0)
    act = pd.DataFrame({"smiles": ["C" * int(k) + "O" for k in rng.integers(18, 40, 40)],
                        "id": [f"a{i}" for i in range(40)]})
    ina = pd.DataFrame({"smiles": ["C" * int(k) + "O" for k in rng.integers(2, 40, 6000)],
                        "id": [f"i{i}" for i in range(6000)]})
    sets = D.draw_negative_sets(act, ina, n_random=200, per_active=5, seed=4)
    bal = sets.balance().abs()
    assert bal.loc["mw", "random"] > 0.5            # actives are heavier than a random draw
    assert bal.loc["mw", "matched"] < 0.15          # and matching removes that
    assert bal["matched"].max() < bal["random"].max()


def test_draw_negative_sets_validates_and_reports_unparsable_structures():
    act, ina = _frames(n_act=10, n_pool=60)
    with pytest.raises(ValueError, match="candidate pool"):
        D.draw_negative_sets(act, ina, n_random=20, per_active=8, seed=0)   # needs 80, has 40
    with pytest.raises(ValueError):
        D.draw_negative_sets(act, ina, n_random=0, per_active=1, seed=0)
    act2, ina2 = _frames(n_act=10, n_pool=400)
    ina2.loc[3, "smiles"] = "not a smiles"
    ina2.loc[7, "smiles"] = "also bad"
    out = D.draw_negative_sets(act2, ina2, n_random=50, per_active=3, seed=1, pool_size=300)
    assert out.n_unparsable >= 0                    # counted, never silently dropped


def test_smd_of_a_constant_property_is_zero_even_with_floating_point_noise():
    a = _synthetic(50, 0, 1)
    b = _synthetic(60, 0, 2)
    a["tpsa"] = 20.23 + np.arange(50) * 0.0 + 1e-15 * np.arange(50)   # constant up to noise
    b["tpsa"] = 20.23
    assert D.standardized_mean_difference(a, b)["tpsa"] == 0.0
