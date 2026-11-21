import json
import threading
import time

import pytest

from vsbench import timecorr


def _info(tmp_path, walls):
    p = tmp_path / "screen_info.json"
    p.write_text(json.dumps({"autodock": {"wall_s": 5.0}, "autodock_seed_wall_s": walls, "n": 3}))
    return p


def test_correction_subtracts_records_and_leaves_other_fields_alone(tmp_path):
    p = _info(tmp_path, {"42": 1769.0, "43": 1772.0, "44": 36807.0})
    ledger = tmp_path / "ledger.json"
    rec = timecorr.apply_correction(p, 44, 35028.0, "host asleep", ledger)
    d = json.loads(p.read_text())
    assert d["autodock_seed_wall_s"] == {"42": 1769.0, "43": 1772.0, "44": 1779.0}
    assert d["autodock"] == {"wall_s": 5.0} and d["n"] == 3
    assert rec["raw_s"] == 36807.0 and rec["corrected_s"] == 1779.0 and rec["reason"] == "host asleep"
    assert d["timing_corrections"] == [rec]
    assert json.loads(ledger.read_text())[0]["subtracted_s"] == 35028.0


def test_applying_twice_does_not_subtract_twice(tmp_path):
    p = _info(tmp_path, {"44": 36807.0})
    timecorr.apply_correction(p, 44, 35028.0, "r")
    timecorr.apply_correction(p, 44, 35028.0, "r")
    d = json.loads(p.read_text())
    assert d["autodock_seed_wall_s"]["44"] == 1779.0 and len(d["timing_corrections"]) == 1


def test_a_correction_larger_than_the_recorded_time_is_refused(tmp_path):
    p = _info(tmp_path, {"44": 900.0})
    with pytest.raises(ValueError):
        timecorr.apply_correction(p, 44, 35028.0, "r")
    assert json.loads(p.read_text())["autodock_seed_wall_s"]["44"] == 900.0


def test_unknown_seed_is_an_error(tmp_path):
    p = _info(tmp_path, {"42": 1.0})
    with pytest.raises(KeyError):
        timecorr.apply_correction(p, 44, 1.0, "r")


def test_waiting_picks_up_a_file_written_later(tmp_path):
    p = tmp_path / "screen_info.json"

    def writer():
        time.sleep(0.4)
        p.write_text("{partial")                                          # a half-written file must be tolerated
        time.sleep(0.4)
        p.write_text(json.dumps({"autodock_seed_wall_s": {"44": 100.0}}))

    t = threading.Thread(target=writer)
    t.start()
    rec = timecorr.wait_and_apply(p, 44, 40.0, "r", None, timeout_s=20, poll_s=0.1)
    t.join()
    assert rec["corrected_s"] == 60.0 and json.loads(p.read_text())["autodock_seed_wall_s"]["44"] == 60.0


def test_waiting_times_out(tmp_path):
    with pytest.raises(TimeoutError):
        timecorr.wait_and_apply(tmp_path / "never.json", 44, 1.0, "r", None, timeout_s=0.3, poll_s=0.1)
