"""The parts that make the harness usable by someone other than its author:
configuration portability, pre-flight data validation and the environment report."""
import json
from pathlib import Path

import pytest

from vsbench import cli
from vsbench import data as D
from vsbench import tools as T
from vsbench.pipeline import Config

MINIMAL = {"data_root": "data", "results_dir": "out", "targets": ["T1"]}


def _mol2(n_res=6):
    """A minimal protein mol2 with readable SUBSTRUCTURE records (chain A, residues 301+)."""
    lines = ["@<TRIPOS>MOLECULE", "x", "@<TRIPOS>ATOM"]
    serial, subs = 1, []
    for i in range(1, n_res + 1):
        x = 4.0 * i
        for name, dx, typ in (("N", 0, "N.am"), ("CA", 1.4, "C.3"), ("C", 2.9, "C.2"), ("O", 2.9, "O.2")):
            lines.append(f"{serial:7d} {name:<4s} {x + dx:9.4f} 0.0000 0.0000 {typ:5s} {i:3d} ALA{300 + i} 0.0")
            serial += 1
        subs.append(f"{i:6d} ALA{300 + i}      1 RESIDUE           4 A     ALA     1 ROOT ALA A {300 + i}")
    return "\n".join(lines + ["@<TRIPOS>BOND", "@<TRIPOS>SUBSTRUCTURE"] + subs) + "\n"


def _write_config(dirpath: Path, **overrides) -> Path:
    dirpath.mkdir(parents=True, exist_ok=True)
    cfg = dict(MINIMAL, **overrides)
    path = dirpath / "cfg.json"
    path.write_text(json.dumps(cfg))
    return path


# ------------------------------------------------------------------ configuration portability
def test_a_relative_path_resolves_against_the_config_file_not_the_working_directory(tmp_path, monkeypatch):
    """The shipped configs carry relative paths so that a checkout works on any machine. Resolving
    them against the current directory instead would silently write results wherever the user
    happened to be standing."""
    path = _write_config(tmp_path / "configs", data_root="../data/raw", results_dir="../results/q1")
    monkeypatch.chdir(tmp_path.parent)                   # deliberately not the repository root
    cfg = Config.from_json(path)
    assert Path(cfg.data_root) == (tmp_path / "data/raw").resolve()
    assert Path(cfg.results_dir) == (tmp_path / "results/q1").resolve()


def test_an_absolute_path_in_a_config_is_left_alone(tmp_path):
    path = _write_config(tmp_path, data_root="/srv/litpcba", results_dir="/srv/out")
    cfg = Config.from_json(path)
    assert (cfg.data_root, cfg.results_dir) == ("/srv/litpcba", "/srv/out")


def test_an_unknown_configuration_key_names_itself(tmp_path):
    """A typo in a hand-edited config otherwise raises a bare TypeError from the dataclass."""
    path = _write_config(tmp_path, n_active=40)           # the real key is n_actives
    with pytest.raises(ValueError, match="n_active"):
        Config.from_json(path)


def test_a_missing_required_configuration_key_names_itself(tmp_path):
    path = tmp_path / "cfg.json"
    path.write_text(json.dumps({"data_root": "d", "results_dir": "o"}))     # no targets
    with pytest.raises(ValueError, match="targets"):
        Config.from_json(path)


def test_a_config_that_is_not_json_says_so(tmp_path):
    path = tmp_path / "cfg.json"
    path.write_text("{not json")
    with pytest.raises(ValueError, match="not valid JSON"):
        Config.from_json(path)


# ------------------------------------------------------------------ data validation
def _target(root: Path, name="T1", n_act=20, n_inact=800, templates=("1abc",), protein=None):
    d = root / name
    d.mkdir(parents=True, exist_ok=True)
    (d / "actives.smi").write_text("".join(f"c1ccccc1{'C' * (i % 5)} ACT{i}\n" for i in range(n_act)))
    (d / "inactives.smi").write_text("".join(f"c1ccncc1{'C' * (i % 7)} INA{i}\n" for i in range(n_inact)))
    for t in templates:
        (d / f"{t}_protein.mol2").write_text(protein if protein is not None else _mol2())
        (d / f"{t}_ligand.mol2").write_text(_mol2(n_res=1))
    return d


def test_a_well_formed_target_directory_validates(tmp_path):
    _target(tmp_path)
    c = D.validate_target_dir(tmp_path / "T1", n_actives=20, per_active_matched=8, n_random=300)
    assert c.ok, c.problems
    assert (c.n_actives, c.n_inactives, c.templates) == (20, 800, ("1abc",))


def test_a_missing_smiles_file_is_reported(tmp_path):
    d = _target(tmp_path)
    (d / "inactives.smi").unlink()
    c = D.validate_target_dir(d)
    assert not c.ok and any("inactives.smi" in p for p in c.problems)


def test_a_protein_without_its_ligand_is_reported(tmp_path):
    """Template selection redocks the co-crystal ligand, so a receptor with no deposited ligand
    cannot be validated and must not be silently skipped."""
    d = _target(tmp_path)
    (d / "1abc_ligand.mol2").unlink()
    c = D.validate_target_dir(d, n_actives=20)
    assert not c.ok
    assert any("1abc_ligand.mol2" in p for p in c.problems)
    assert any("no <pdb>_protein.mol2" in p for p in c.problems)


def test_a_receptor_whose_residues_all_read_as_unk_is_caught_before_any_docking(tmp_path):
    """Correction 1's failure mode. Such a receptor docks and scores normally, but pocket
    detection reports a shattered surface and the interaction pass finds nothing, both without
    raising. Catching it here costs seconds; catching it afterwards cost a production run."""
    headless = _mol2().split("@<TRIPOS>SUBSTRUCTURE")[0] + "@<TRIPOS>SUBSTRUCTURE\n"
    d = _target(tmp_path, protein=headless)
    c = D.validate_target_dir(d, n_actives=20)
    assert not c.ok
    assert any("UNK" in p for p in c.problems)


def test_an_inactive_pool_too_small_for_the_requested_matching_is_reported(tmp_path):
    d = _target(tmp_path, n_act=20, n_inact=200)
    c = D.validate_target_dir(d, n_actives=20, per_active_matched=8, n_random=300)
    assert not c.ok
    assert any("460 are needed" in p for p in c.problems)      # 300 random + 20 x 8 matched


def test_an_unparsable_active_is_counted(tmp_path):
    d = _target(tmp_path, n_act=12)
    with open(d / "actives.smi", "a") as fh:
        fh.write("not_a_smiles BAD1\n")
    c = D.validate_target_dir(d, n_actives=13)
    assert c.n_unparsable == 1 and any("cannot be built by RDKit" in p for p in c.problems)


def test_too_few_actives_for_a_bootstrap_is_reported(tmp_path):
    d = _target(tmp_path, n_act=5)
    c = D.validate_target_dir(d, n_actives=5)
    assert not c.ok and any("too few" in p for p in c.problems)


def test_validate_data_root_walks_every_subdirectory_when_no_target_is_named(tmp_path):
    _target(tmp_path, "T1")
    _target(tmp_path, "T2")
    checks = D.validate_data_root(tmp_path, None, n_actives=20)
    assert [c.name for c in checks] == ["T1", "T2"]
    assert all(c.ok for c in checks)


def test_a_missing_data_root_is_refused(tmp_path):
    with pytest.raises(FileNotFoundError):
        D.validate_data_root(tmp_path / "nope")


# ------------------------------------------------------------------ environment report
def test_check_tools_marks_the_required_ones(monkeypatch):
    monkeypatch.setattr(T, "find_tool", lambda name: None)
    statuses = T.check_tools()
    assert {s.name for s in statuses if s.required} == set(T.REQUIRED_TOOLS)
    assert all(not s.ok for s in statuses)
    assert all(s.purpose for s in statuses)


def test_every_python_dependency_reports_a_purpose():
    rows = T.check_python_packages()
    assert all(purpose for _, _, _, purpose in rows)
    assert dict((n, req) for n, _, req, _ in rows)["numpy"] is True
    assert dict((n, req) for n, _, req, _ in rows)["admet_ai"] is False


def test_doctor_fails_when_a_required_executable_is_missing(monkeypatch, capsys):
    monkeypatch.setattr(T, "find_tool", lambda name: None if name == "gnina" else Path("/usr/bin") / name)
    monkeypatch.setattr(T, "gpu_report", lambda: "test GPU")
    assert cli.main(["doctor"]) == 1
    out = capsys.readouterr().out
    assert "gnina" in out and "MISSING" in out and "VSBENCH_BIN" in out


def test_doctor_succeeds_when_every_required_executable_resolves(monkeypatch, capsys):
    monkeypatch.setattr(T, "find_tool", lambda name: Path("/usr/bin") / name)
    monkeypatch.setattr(T, "gpu_report", lambda: "test GPU")
    assert cli.main(["doctor"]) == 0
    assert "All required components present." in capsys.readouterr().out


def test_a_missing_optional_executable_does_not_fail_the_doctor(monkeypatch, capsys):
    monkeypatch.setattr(T, "find_tool",
                        lambda name: None if name in T.OPTIONAL_TOOLS else Path("/usr/bin") / name)
    monkeypatch.setattr(T, "gpu_report", lambda: "test GPU")
    assert cli.main(["doctor"]) == 0


def test_gpu_report_says_so_when_nvidia_smi_is_absent(monkeypatch):
    monkeypatch.setattr(T.shutil, "which", lambda name: None)
    assert "nvidia-smi not found" in T.gpu_report()


# ------------------------------------------------------------------ the check-data command
def test_check_data_exits_nonzero_when_any_target_is_unusable(tmp_path, capsys):
    _target(tmp_path, "Good", n_act=20)
    bad = _target(tmp_path, "Bad", n_act=20)
    (bad / "actives.smi").unlink()
    code = cli.main(["check-data", "--data-root", str(tmp_path), "--n-actives", "20"])
    out = capsys.readouterr().out
    assert code == 1
    assert "FAIL  Bad" in out and "ok    Good" in out and "1 of 2 target(s) usable" in out


def test_check_data_warns_that_the_decoy_arm_needs_a_second_target(tmp_path, capsys):
    _target(tmp_path, "Only", n_act=20)
    assert cli.main(["check-data", "--data-root", str(tmp_path), "--n-actives", "20"]) == 0
    assert "at least two target directories" in capsys.readouterr().out


def test_check_data_takes_its_targets_and_ratios_from_a_config(tmp_path, capsys):
    _target(tmp_path / "data", "T1", n_act=20, n_inact=200)
    _target(tmp_path / "data", "T2", n_act=20)
    path = _write_config(tmp_path / "configs", data_root="../data", results_dir="../out",
                         targets=["T1"], n_actives=20, per_active_matched=8, n_random=300)
    assert cli.main(["check-data", "--config", str(path)]) == 1
    out = capsys.readouterr().out
    assert "T1" in out and "T2" not in out                 # only the configured target is checked
    assert "460 are needed" in out                         # the config's own ratios, not the defaults


def test_check_data_needs_a_data_root_or_a_config(capsys):
    with pytest.raises(SystemExit):
        cli.main(["check-data"])
    assert "--data-root or --config" in capsys.readouterr().err


def test_the_cli_reports_its_version(capsys):
    with pytest.raises(SystemExit) as e:
        cli.main(["--version"])
    assert e.value.code == 0
    assert cli.__version__ in capsys.readouterr().out
