"""Remove a known host-suspend interval from a recorded wall time, and keep a ledger of what was removed.

When the host sleeps, a wall-clock timer keeps counting while no compute happens. The engine timers now use a
monotonic clock, but runs started before that change carry the gap in their ``screen_info.json``. This module
subtracts a documented interval from one AutoDock-GPU seed timing and records the correction beside the raw
value, so nothing is edited silently and a second application is a no-op.
"""
from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path


def apply_correction(info_path: Path, seed: int, subtract_s: float, reason: str, ledger: Path | None = None) -> dict:
    """Subtract ``subtract_s`` from ``autodock_seed_wall_s[seed]``. Idempotent. Returns the correction record."""
    info_path = Path(info_path)
    info = json.loads(info_path.read_text())
    key = str(seed)
    walls = info.get("autodock_seed_wall_s", {})
    if key not in walls:
        raise KeyError(f"no autodock_seed_wall_s entry for seed {seed} in {info_path}")
    done = info.setdefault("timing_corrections", [])
    for rec in done:
        if rec["field"] == f"autodock_seed_wall_s.{key}" and abs(rec["subtracted_s"] - subtract_s) < 1e-6:
            return rec                                                   # already applied
    raw = float(walls[key])
    corrected = raw - subtract_s
    if corrected < 0:
        raise ValueError(f"correction of {subtract_s:.0f} s exceeds the recorded {raw:.0f} s; refusing")
    rec = {"field": f"autodock_seed_wall_s.{key}", "raw_s": raw, "subtracted_s": float(subtract_s),
           "corrected_s": corrected, "reason": reason}
    walls[key] = corrected
    done.append(rec)
    tmp = info_path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(info))
    os.replace(tmp, info_path)                                           # atomic: readers never see a partial file
    if ledger is not None:
        ledger = Path(ledger)
        entries = json.loads(ledger.read_text()) if ledger.exists() else []
        entries.append({"file": str(info_path), **rec})
        ledger.write_text(json.dumps(entries, indent=2))
    return rec


def wait_and_apply(info_path: Path, seed: int, subtract_s: float, reason: str, ledger: Path | None,
                   timeout_s: float, poll_s: float = 30.0) -> dict:
    """Wait for a running screen to write its info file, then apply the correction."""
    deadline = time.monotonic() + timeout_s
    while True:
        try:
            if str(seed) in json.loads(Path(info_path).read_text()).get("autodock_seed_wall_s", {}):
                return apply_correction(info_path, seed, subtract_s, reason, ledger)
        except (FileNotFoundError, json.JSONDecodeError):
            pass                                                         # not written yet, or mid-write
        if time.monotonic() > deadline:
            raise TimeoutError(f"{info_path} did not appear within {timeout_s:.0f} s")
        time.sleep(poll_s)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--info", required=True, type=Path)
    ap.add_argument("--seed", required=True, type=int)
    ap.add_argument("--subtract-s", required=True, type=float)
    ap.add_argument("--reason", required=True)
    ap.add_argument("--ledger", type=Path)
    ap.add_argument("--wait-s", type=float, default=0.0, help="wait up to this long for the info file")
    a = ap.parse_args(argv)
    rec = (wait_and_apply(a.info, a.seed, a.subtract_s, a.reason, a.ledger, a.wait_s) if a.wait_s
           else apply_correction(a.info, a.seed, a.subtract_s, a.reason, a.ledger))
    print(json.dumps(rec))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
