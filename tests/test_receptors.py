"""Receptor-level studies: options, path resolution, generality and failure handling."""
from __future__ import annotations

import dataclasses
import json
from pathlib import Path

import pandas as pd
import pytest

from vsbench import receptors
from vsbench.pipeline import Config


def _cfg(tmp_path, **receptor_opts) -> Config:
    data = tmp_path / "data" / "raw" / "set"
    (data / "TGT").mkdir(parents=True)
    (data / "TGT" / "actives.smi").write_text("CCO act1\n")
    (data / "TGT" / "inactives.smi").write_text("CCC ina1\n")
    cfg = {"data_root": str(data), "results_dir": str(tmp_path / "results" / "q1"),
           "targets": ["TGT"]}
    if receptor_opts:
        cfg["receptors"] = receptor_opts
    p = tmp_path / "cfg.json"
    p.write_text(json.dumps(cfg))
    return Config.from_json(p)


# ------------------------------------------------------------------ options
def test_defaults_prefer_local_detectors_and_never_the_web_service(tmp_path):
    opt = receptors.options(_cfg(tmp_path))
    assert opt["fpocket"] and opt["controls"]
    # the one step that leaves the machine stays off unless explicitly asked for
    assert opt["dogsite"] is False
    # the local ranker is resolved against this machine, never left as the literal "auto"
    assert opt["p2rank"] in (True, False)
    assert not opt["offline"], "offline is opt-in; the default still allows a cache fill"


def test_config_overrides_merge_over_the_defaults(tmp_path):
    opt = receptors.options(_cfg(tmp_path, dogsite=True, n_jobs=7))
    assert opt["dogsite"] and opt["n_jobs"] == 7
    assert opt["fpocket"], "an override must not drop the keys it did not mention"


# ------------------------------------------------------------- path resolution
def test_tables_default_beside_the_run_not_inside_it(tmp_path):
    cfg = _cfg(tmp_path)
    assert receptors.receptor_dir(cfg) == Path(cfg.results_dir).parent


def test_receptor_dir_is_configurable(tmp_path):
    cfg = _cfg(tmp_path, dir=str(tmp_path / "elsewhere"))
    assert receptors.receptor_dir(cfg) == tmp_path / "elsewhere"
    assert (tmp_path / "elsewhere").is_dir()


def test_pdb_cache_reuses_an_existing_one_above_the_data_root(tmp_path):
    cfg = _cfg(tmp_path)
    shared = tmp_path / "data" / "pdb_cache"       # two levels above data_root
    shared.mkdir(parents=True)
    assert receptors.pdb_cache(cfg) == shared


def test_pdb_cache_is_created_beside_the_data_root_when_none_exists(tmp_path):
    cfg = _cfg(tmp_path)
    got = receptors.pdb_cache(cfg)
    assert got == Path(cfg.data_root).resolve().parent / "pdb_cache" and got.is_dir()


def test_explicit_pdb_cache_wins(tmp_path):
    (tmp_path / "data" / "pdb_cache").mkdir(parents=True)
    cfg = _cfg(tmp_path, pdb_cache=str(tmp_path / "mine"))
    assert receptors.pdb_cache(cfg) == tmp_path / "mine"


# ------------------------------------------------------------------ generality
@pytest.mark.parametrize("name,is_code", [
    ("4agq", True), ("1HXB", True), ("5o1i", True),
    ("my_receptor", False), ("abcd", False), ("12ab3", False), ("apo", False),
])
def test_only_real_pdb_codes_are_sent_to_the_public_service(name, is_code):
    assert bool(receptors.PDB_CODE.match(name)) is is_code


def test_all_targets_lists_every_target_directory(tmp_path):
    cfg = _cfg(tmp_path)
    root = Path(cfg.data_root)
    for extra in ("AAA", "ZZZ"):
        (root / extra).mkdir()
        (root / extra / "actives.smi").write_text("CCO a\n")
    (root / "not_a_target").mkdir()          # no actives.smi, so not a target
    assert receptors.all_targets(cfg) == ["AAA", "TGT", "ZZZ"]


def test_mol2_heavy_xyz_skips_hydrogens_and_survives_an_empty_file(tmp_path):
    p = tmp_path / "l.mol2"
    p.write_text("@<TRIPOS>ATOM\n"
                 "1 C1 1.0 2.0 3.0 C.3 1 LIG 0.0\n"
                 "2 H1 9.0 9.0 9.0 H 1 LIG 0.0\n"
                 "@<TRIPOS>BOND\n")
    xyz = receptors.mol2_heavy_xyz(p)
    assert xyz.shape == (1, 3) and list(xyz[0]) == [1.0, 2.0, 3.0]
    empty = tmp_path / "e.mol2"
    empty.write_text("@<TRIPOS>ATOM\n@<TRIPOS>BOND\n")
    assert receptors.mol2_heavy_xyz(empty).shape == (0, 3)


# ------------------------------------------------------------ failure handling
def test_an_unreadable_template_becomes_a_row_and_does_not_raise(tmp_path):
    """One bad template must not discard the detections already made for the others."""
    cfg = _cfg(tmp_path)
    rows = receptors._site_rows(cfg, receptors.DEFAULTS, tmp_path, tmp_path,
                                ("TGT", "does_not_exist"))
    assert len(rows) == 1
    assert rows[0]["method"] == "setup" and rows[0]["error"]
    assert rows[0]["template"] == "does_not_exist"


def test_metal_rows_report_the_error_rather_than_raising(tmp_path):
    cfg = _cfg(tmp_path)
    rows = receptors._metal_rows(cfg, receptors.DEFAULTS, ("TGT", "nope"))
    assert len(rows) == 1 and "error" in rows[0]


def test_studies_on_no_templates_return_empty_without_writing(tmp_path):
    cfg = _cfg(tmp_path)
    msgs = []
    for fn in (receptors.qc_study, receptors.sites_study, receptors.metal_study):
        assert fn(cfg, [], log=msgs.append).empty
    assert all("no templates" in m for m in msgs)


# ----------------------------------------------------------- druggability rule
def test_druggability_is_informative_only_when_the_negative_ranks_below_every_positive():
    good = pd.DataFrame({
        "method": ["d", "d", "d"], "control_role": ["positive", "positive", "negative"],
        "site_pocket_druggability": [0.80, 0.64, 0.30]})
    assert receptors.druggability_verdict(good)["d"]["informative"] is True

    # the real case: the negative control outscored a positive, so the score says nothing here
    bad = good.copy()
    bad.loc[2, "site_pocket_druggability"] = 0.74
    verdict = receptors.druggability_verdict(bad)["d"]
    assert verdict["informative"] is False
    assert verdict["negative_control"] == 0.74 and verdict["lowest_positive"] == 0.64


def test_druggability_verdict_is_empty_without_controls():
    assert receptors.druggability_verdict(pd.DataFrame()) == {}
    assert receptors.druggability_verdict(None) == {}


# ------------------------------------------------------------------------- CLI
def test_cli_registers_the_three_receptor_commands():
    from vsbench import cli
    names = {n for n, _ in cli.STUDIES}
    assert {"qc", "sites", "metals"} <= names
    assert set(cli.RECEPTOR_STUDIES) == {"qc", "sites", "metals"}


def test_cli_all_targets_flag_widens_beyond_the_configured_targets(tmp_path, monkeypatch):
    from vsbench import cli
    cfg = _cfg(tmp_path)
    root = Path(cfg.data_root)
    (root / "OTHER").mkdir()
    (root / "OTHER" / "actives.smi").write_text("CCO a\n")
    seen = {}
    monkeypatch.setattr(receptors, "qc_study",
                        lambda c, t, *a, **k: seen.setdefault("targets", list(t)) or pd.DataFrame())
    cli.main(["qc", "--config", str(tmp_path / "cfg.json"), "--all-targets"])
    assert seen["targets"] == ["OTHER", "TGT"]
    seen.clear()
    cli.main(["qc", "--config", str(tmp_path / "cfg.json")])
    assert seen["targets"] == ["TGT"]


def test_receptors_config_key_round_trips_through_json(tmp_path):
    cfg = _cfg(tmp_path, dogsite=True, dir=str(tmp_path / "r"))
    assert cfg.receptors["dogsite"] is True
    assert dataclasses.asdict(cfg)["receptors"]["dir"].endswith("r")


# ------------------------------------------------------ local-only operation
def test_p2rank_auto_enables_only_when_both_binaries_are_found(tmp_path, monkeypatch):
    from vsbench import tools
    cfg = _cfg(tmp_path, p2rank="auto")

    monkeypatch.setattr(tools, "find_p2rank", lambda: Path("/x/prank"))
    monkeypatch.setattr(tools, "find_java", lambda: Path("/x/java"))
    opt = receptors.options(cfg)
    assert opt["p2rank"] is True
    assert opt["p2rank_bin"] == "/x/prank" and opt["java_bin"] == "/x/java"

    monkeypatch.setattr(tools, "find_java", lambda: None)      # P2Rank without a JRE is unusable
    assert receptors.options(cfg)["p2rank"] is False


def test_an_explicit_p2rank_request_is_not_silently_disabled(tmp_path, monkeypatch):
    """'auto' may resolve to off; an explicit True must stay on so the run reports the failure
    instead of quietly producing a table with one detector missing."""
    from vsbench import tools
    monkeypatch.setattr(tools, "find_p2rank", lambda: None)
    monkeypatch.setattr(tools, "find_java", lambda: None)
    assert receptors.options(_cfg(tmp_path, p2rank=True))["p2rank"] is True


def test_offline_forces_off_everything_that_would_reach_the_network(tmp_path):
    opt = receptors.options(_cfg(tmp_path, offline=True, dogsite=True))
    assert opt["offline"] is True
    assert opt["dogsite"] is False, "offline must override an explicit dogsite request"


def test_offline_deposition_fetch_refuses_rather_than_reaching_out(tmp_path, monkeypatch):
    import urllib.request

    from vsbench import metals
    monkeypatch.setattr(urllib.request, "urlopen",
                        lambda *a, **k: pytest.fail("offline mode made a network request"))
    with pytest.raises(FileNotFoundError) as e:
        metals.fetch_deposition("1ABC", tmp_path / "cache", offline=True)
    # the message has to say what to do about it, not just that it failed
    assert "1ABC.pdb" in str(e.value) and str(tmp_path / "cache") in str(e.value)


def test_a_cached_deposition_is_used_offline_without_any_request(tmp_path, monkeypatch):
    import urllib.request

    from vsbench import metals
    cache = tmp_path / "cache"
    cache.mkdir()
    (cache / "1ABC.pdb").write_text("x" * 2000)
    monkeypatch.setattr(urllib.request, "urlopen",
                        lambda *a, **k: pytest.fail("a cached hit must not reach the network"))
    assert metals.fetch_deposition("1abc", cache, offline=True) == cache / "1ABC.pdb"


def test_structure_qc_reads_depositions_from_cache_only(tmp_path):
    """qc never downloads: it uses what is cached and leaves the rest blank."""
    cache = tmp_path / "c"
    cache.mkdir()
    assert receptors._cached_deposition("9zzz", cache) is None
    (cache / "9ZZZ.cif").write_text("data_9ZZZ\n")
    assert receptors._cached_deposition("9zzz", cache) == cache / "9ZZZ.cif"


def test_tools_finds_java_through_java_home(tmp_path, monkeypatch):
    from vsbench import tools
    jbin = tmp_path / "jdk" / "bin"
    jbin.mkdir(parents=True)
    (jbin / "java").write_text("#!/bin/sh\n")
    (jbin / "java").chmod(0o755)
    monkeypatch.setenv("JAVA_HOME", str(tmp_path / "jdk"))
    assert tools.find_java() == jbin / "java"


def test_tools_finds_p2rank_through_its_env_var(tmp_path, monkeypatch):
    from vsbench import tools
    d = tmp_path / "p2rank_9.9"
    d.mkdir()
    (d / "prank").write_text("#!/bin/sh\n")
    (d / "prank").chmod(0o755)
    monkeypatch.setenv("VSBENCH_P2RANK", str(d))          # a directory resolves to its launcher
    assert tools.find_p2rank() == d / "prank"
    monkeypatch.setenv("VSBENCH_P2RANK", str(d / "prank"))
    assert tools.find_p2rank() == d / "prank"


def test_p2rank_and_java_are_optional_not_required():
    from vsbench import tools
    assert {"prank", "java"} <= set(tools.OPTIONAL_TOOLS)
    assert not ({"prank", "java"} & set(tools.REQUIRED_TOOLS)), \
        "a missing JRE must not stop a docking run"


# -------------------------------------------------- paired box-strategy tests
def _grid_rows(pattern):
    """pattern: {template: {strategy: success}} for one engine."""
    return [{"template": t, "strategy": s, "ok": True,
             "gn_success": v, "ad_success": v}
            for t, d in pattern.items() for s, v in d.items()]


def test_mcnemar_matches_the_hand_computed_exact_value():
    """12 discordant one way against 1 the other is p = 0.0034; this is the published value."""
    from vsbench import gridstrat as G
    rows = _grid_rows({f"a{i}": {"ligand": True, "pocket": False} for i in range(12)}
                      | {"b0": {"ligand": False, "pocket": True}})
    t = G.paired_strategies(rows, engines=("gn",))
    r = t.iloc[0]
    assert r.a_only == 12 and r.b_only == 1 and r.discordant == 13
    assert round(r.p_mcnemar_exact, 4) == 0.0034
    assert r.favours == "ligand"


def test_no_discordant_pairs_is_p_one_and_unresolved():
    from vsbench import gridstrat as G
    rows = _grid_rows({f"a{i}": {"ligand": True, "pocket": True} for i in range(10)})
    r = G.paired_strategies(rows, engines=("gn",)).iloc[0]
    assert r.discordant == 0 and r.p_mcnemar_exact == 1.0 and r.favours == "not resolved"


def test_pairwise_scope_keeps_a_comparison_when_a_strategy_is_added():
    """The reason pairwise is the default: adding a strategy must not silently restate an
    already-published comparison on fewer templates."""
    from vsbench import gridstrat as G
    base = _grid_rows({f"t{i}": {"ligand": True, "pocket": False} for i in range(6)})
    before = G.paired_strategies(base, engines=("gn",)).iloc[0]
    # a third strategy that ran on only one template
    extended = base + _grid_rows({"t0": {"blind": True}})
    after = next(r for _, r in G.paired_strategies(extended, engines=("gn",)).iterrows()
                 if {r.strategy_a, r.strategy_b} == {"ligand", "pocket"})
    assert (after.n_paired, after.p_mcnemar_exact) == (before.n_paired, before.p_mcnemar_exact)

    listwise = G.paired_strategies(extended, engines=("gn",), scope="listwise")
    lw = next(r for _, r in listwise.iterrows()
              if {r.strategy_a, r.strategy_b} == {"ligand", "pocket"})
    assert lw.n_paired == 1, "listwise collapses onto the templates all three share"


def test_templates_that_failed_to_prepare_are_excluded_not_scored_as_failures():
    from vsbench import gridstrat as G
    rows = _grid_rows({f"t{i}": {"ligand": True, "pocket": False} for i in range(4)})
    rows.append({"template": "bad", "strategy": "pocket", "ok": False,
                 "gn_success": False, "ad_success": False})
    rows.append({"template": "bad", "strategy": "ligand", "ok": True,
                 "gn_success": True, "ad_success": True})
    r = G.paired_strategies(rows, engines=("gn",)).iloc[0]
    assert r.n_paired == 4, "a strategy that could not be prepared is not a pose failure"


def test_p2rank_is_a_registered_box_strategy():
    from vsbench import gridstrat as G
    assert "p2rank" in G.STRATEGIES
    assert G.STRATEGIES.index("pocket") < G.STRATEGIES.index("p2rank") < G.STRATEGIES.index("blind")


def test_p2rank_pocket_centre_is_parsed_for_box_seeding():
    """The parser used to drop center_x/y/z, which silently made a P2Rank box impossible."""
    from vsbench import gridstrat as G, sites as S
    csv = ("name,rank,score,probability,sas_points,surf_atoms,center_x,center_y,center_z,"
           "residue_ids,surf_atom_ids\n"
           "pocket1,1,5.58,0.267,59,28,10.5,20.25,-30.0,A_158 A_160,1 2\n")
    pk = S.parse_p2rank_predictions(csv)
    assert pk[0].extra["center"] == (10.5, 20.25, -30.0)
    assert G.top_pocket_center(pk) == (10.5, 20.25, -30.0)
    assert G.pocket_box(G.top_pocket_center(pk)).center == (10.5, 20.25, -30.0)


# ------------------------------- retrospective vs deployable similarity baseline
def test_deployable_baseline_uses_no_evaluation_label(tmp_path):
    """The reviewer's objection in code: the oracle references the evaluation actives, the
    deployable one must reference only compounds held out of the evaluation set."""
    import pytest as _pytest

    from vsbench import similarity
    evaluation = ["CCO", "CCCO", "c1ccccc1"]
    actives = ["CCO", "CCCO"]                 # the reference IS the flagged subset, in order
    held_out = ["CCCCO"]
    is_active = [True, True, False]

    oracle = similarity.reference_similarity(evaluation, actives, query_is_reference=is_active)
    deployable = similarity.reference_similarity(evaluation, held_out)
    assert len(oracle) == len(deployable) == 3
    # leave-one-out: an active is never scored 1.0 against itself
    assert all(v < 1.0 for v in oracle[:2]), oracle
    # the deployable score cannot depend on evaluation labels, because it never sees them
    assert list(deployable) == list(similarity.reference_similarity(evaluation, held_out))

    # a reference set that is not the flagged subset must raise, not silently return self-matches
    with _pytest.raises(ValueError, match="leave-one-out"):
        similarity.reference_similarity(evaluation, evaluation, query_is_reference=is_active)


def test_held_out_actives_exclude_everything_in_the_evaluation_set(tmp_path, monkeypatch):
    import pandas as pd
    from vsbench import data, leakage
    cfg = _cfg(tmp_path)
    fake = data.Target(name="TGT", root=Path(cfg.data_root) / "TGT",
                       actives=pd.DataFrame({"smiles": ["C", "CC", "CCC"],
                                             "id": ["a1", "a2", "a3"]}),
                       inactives=pd.DataFrame({"smiles": [], "id": []}), templates=())
    monkeypatch.setattr(data, "load_target", lambda *a, **k: fake)
    held = leakage.held_out_actives(cfg, "TGT", ["a1"])
    assert list(held["id"]) == ["a2", "a3"]
    assert leakage.held_out_actives(cfg, "TGT", ["a1", "a2", "a3"]).empty


def test_leakage_registered_scores_are_the_three_primary_ones():
    from vsbench import analysis, leakage
    assert set(leakage.REGISTERED_SCORES) <= set(analysis.SCORE_COLUMNS)
    assert "consensus" in leakage.REGISTERED_SCORES
