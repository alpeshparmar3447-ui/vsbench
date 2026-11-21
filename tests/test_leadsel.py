"""Pipeline helpers for lead selection, with numeric compound ids (LIT-PCBA ids are integers)."""
import json

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from vsbench import pipeline


def _sdf(path, names):
    w = Chem.SDWriter(str(path))
    for n in names:
        m = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        AllChem.EmbedMolecule(m, randomSeed=1)
        m.SetProp("_Name", str(n))
        w.write(m)
        m2 = Chem.Mol(m)                                   # a second, lower-ranked pose of the same ligand
        m2.SetProp("_Name", str(n))
        m2.SetProp("rank", "2")
        w.write(m2)
    w.close()


def test_top_poses_match_numeric_ids_and_keep_the_first_pose(tmp_path):
    _sdf(tmp_path / "gnina_00000.out.sdf", [845212, 848082])
    _sdf(tmp_path / "gnina_00001.out.sdf", [14724824])
    got = pipeline.gnina_top_poses(tmp_path, {845212, "14724824"})          # ints and strings both work
    assert set(got) == {"845212", "14724824"}
    assert not got["845212"].HasProp("rank")                                 # the top pose, not the second


def _cfg(tmp_path):
    return pipeline.Config(data_root=str(tmp_path), results_dir=str(tmp_path / "r"), targets=["T"])


def test_admet_is_cached_and_numeric_ids_survive_the_round_trip(tmp_path, monkeypatch):
    from vsbench import admet
    calls = []

    def fake(smiles):
        calls.append(len(smiles))
        return pd.DataFrame({e: np.full(len(smiles), 0.1) for e in admet.ENDPOINTS}, index=smiles)
    monkeypatch.setattr(admet, "predict", fake)
    cfg = _cfg(tmp_path)
    cfg.target_dir("T").mkdir(parents=True)
    df = pd.DataFrame({"name": [845212, 848082], "smiles": ["CC(=O)Oc1ccccc1C(=O)O", "CCO"]})
    a = pipeline.run_admet(cfg, "T", df, log=lambda *_: None)
    assert list(a["name"]) == ["845212", "848082"] and (a["tier"] == "PASS").all()
    b = pipeline.run_admet(cfg, "T", df, log=lambda *_: None)                  # second call reads the cache
    assert calls == [2] and list(b["name"]) == ["845212", "848082"]
    c = pipeline.run_admet(cfg, "T", pd.concat([df, pd.DataFrame({"name": [1], "smiles": ["CCN"]})]),
                           log=lambda *_: None)                                # a new compound forces a refresh
    assert calls == [2, 3] and len(c) == 3


def test_config_round_trips_the_lead_selection_block(tmp_path):
    cfg = _cfg(tmp_path)
    assert cfg.lead_selection["top_k"] == 100
    p = tmp_path / "c.json"
    p.write_text(json.dumps({"data_root": "d", "results_dir": "r", "targets": ["T"], "lead_selection": None}))
    assert pipeline.Config.from_json(p).lead_selection is None
