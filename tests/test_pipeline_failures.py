"""A failing target must not stop the others (no GPU needed: stages are stubbed)."""
import json
from pathlib import Path

import pandas as pd
import pytest

from vsbench import cli, pipeline
from vsbench.redock import RedockResult


@pytest.fixture
def cfg_path(tmp_path):
    cfg = {"data_root": str(tmp_path), "results_dir": str(tmp_path / "res"),
           "targets": ["BAD", "GOOD"], "n_boot": 20}
    p = tmp_path / "cfg.json"
    p.write_text(json.dumps(cfg))
    return p


def test_a_failing_target_is_recorded_and_the_rest_still_run(cfg_path, monkeypatch):
    calls = []

    def fake_templates(cfg, target, log=print):
        calls.append(("templates", target))
        if target == "BAD":
            raise ValueError("no template could be redocked")
        return [], RedockResult(template="1abc", ok=True)

    monkeypatch.setattr(pipeline, "run_templates", fake_templates)
    monkeypatch.setattr(pipeline, "run_screen", lambda cfg, t, tmpl, log=print: calls.append(("screen", t, tmpl)))
    monkeypatch.setattr(pipeline, "run_analysis", lambda cfg, t, df=None, log=print: calls.append(("analyze", t)))
    monkeypatch.setattr(pipeline, "aggregate", lambda cfg: calls.append(("aggregate", tuple(cfg.targets))))

    cfg = pipeline.Config.from_json(cfg_path)
    failed = pipeline.run_all(cfg)
    assert failed == ["BAD"]
    assert ("screen", "GOOD", "1abc") in calls and ("analyze", "GOOD") in calls
    assert ("screen", "BAD", "1abc") not in calls
    assert ("aggregate", ("GOOD",)) in calls                  # only successful targets are aggregated
    err = (cfg.target_dir("BAD") / "error.txt").read_text()
    assert "no template could be redocked" in err
    assert "FAILED" in (pipeline.Path(cfg.results_dir) / "run.log").read_text()


def test_cli_exit_code_reflects_failures(cfg_path, monkeypatch):
    monkeypatch.setattr(pipeline, "run_templates", lambda *a, **k: (_ for _ in ()).throw(RuntimeError("boom")))
    assert cli.main(["run", "--config", str(cfg_path), "--stage", "templates"]) == 1
    monkeypatch.setattr(pipeline, "run_templates", lambda *a, **k: ([], RedockResult("t", ok=True)))
    assert cli.main(["run", "--config", str(cfg_path), "--stage", "templates"]) == 0


def test_interaction_recovery_is_written_in_batches_so_an_interruption_keeps_progress(tmp_path, monkeypatch):
    """A long PLIP pass that only writes at the end loses everything when it is killed. Each batch must
    land in the cache file, and a later call must pick up only what is still missing."""
    import json

    import numpy as np
    from vsbench import interactions, pipeline, prioritize, redock as R, sites

    tdir = tmp_path / "T"
    (tdir / "redock" / "1abc").mkdir(parents=True)
    (tdir / "selected.json").write_text(json.dumps({"template": "1abc"}))

    cfg = pipeline.Config(data_root=str(tmp_path), results_dir=str(tmp_path), targets=["T"],
                          lead_selection={"top_k": 2, "n_jobs": 1, "interaction_batch": 2})
    monkeypatch.setattr(pipeline.Config, "target_dir", lambda self, t: tdir)
    monkeypatch.setattr(pipeline, "data", type("D", (), {
        "load_target": staticmethod(lambda root, name: type("T", (), {
            "protein_mol2": lambda self, s: tmp_path / "p.mol2",
            "ligand_mol2": lambda self, s: tmp_path / "l.mol2"})())})) 
    monkeypatch.setattr(sites, "write_pdb_from_mol2", lambda src, dst: Path(dst).write_text("ATOM\n"))
    monkeypatch.setattr(R, "crystal_ligand", lambda lig, wd: (object(), ""))
    monkeypatch.setattr(interactions, "reference_set", lambda prot, xtal: frozenset({("hb", "A", "1")}))
    monkeypatch.setattr(pipeline, "gnina_top_poses", lambda wd, todo: {n: object() for n in sorted(todo)})
    monkeypatch.setattr(prioritize, "select_pool", lambda sub, k: sub.head(0))

    seen, writes = [], []

    def fake_score(prot, ref, poses, n_jobs=8):
        seen.append(sorted(poses))
        writes.append(len(pd.read_csv(tdir / "interactions.csv")) if (tdir / "interactions.csv").exists() else 0)
        return pd.DataFrame({"name": sorted(poses), "ok": True, "recall": 0.5})

    monkeypatch.setattr(interactions, "score_poses", fake_score)
    df = pd.DataFrame({"name": ["a", "b", "c", "d", "e"], "neg_set": ["active"] * 5, "label": [1] * 5})
    out = pipeline.run_interactions(cfg, "T", df, {}, log=lambda *_: None, all_compounds=True)

    assert seen == [["a", "b"], ["c", "d"], ["e"]]           # batched, not one call for everything
    assert writes == [0, 2, 4]                               # each batch was on disk before the next ran
    assert sorted(out["name"]) == ["a", "b", "c", "d", "e"]

    seen.clear()
    pipeline.run_interactions(cfg, "T", df, {}, log=lambda *_: None, all_compounds=True)
    assert seen == []                                        # nothing recomputed on a resume
