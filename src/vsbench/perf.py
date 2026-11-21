"""Controlled performance measurements (plan section 4.13).

Every case docks the same ligands into the same receptor and box, in a fresh directory, after one warm-up
run per engine (CUDA context and CNN weights load once). Wall time is the engine call alone: ligand
preparation and file I/O outside the engine are excluded and reported elsewhere. GPU memory and
utilisation are sampled from nvidia-smi every 0.5 s while a case runs.

Cases
  autodock-gpu           GPU
  gnina (CNN, GPU)       GPU
  gnina (CNN, CPU)       gnina with --no_gpu
  gnina (Vina, CPU)      gnina with --cnn_scoring none
  vina (CPU)             AutoDock Vina, one process per ligand, as many in parallel as threads: the most
                         favourable arrangement for a CPU baseline. Used only as a throughput reference.
"""
from __future__ import annotations

import shutil
import statistics
import subprocess
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from . import engines
from .tools import require_tool, run_tool, tool_env

CASES = ("autodock-gpu", "gnina (CNN, GPU)", "gnina (CNN, CPU)", "gnina (Vina, CPU)", "vina (CPU)")


def parse_smi(line: str) -> tuple[float, float]:
    """'1234 MiB, 56 %' style nvidia-smi csv line -> (memory MiB, utilisation %)."""
    mem, util = [p.strip().split()[0] for p in line.strip().split(",")[:2]]
    return float(mem), float(util)


class GpuSampler:
    """Context manager sampling GPU memory and utilisation. Peak memory is reported above the idle baseline."""

    def __init__(self, interval: float = 0.5, index: int = 0):
        self.interval, self.index = interval, index
        self.mem, self.util = [], []
        self._stop = threading.Event()
        self._t = None
        self.available = shutil.which("nvidia-smi") is not None

    def _read(self):
        out = subprocess.run(["nvidia-smi", "-i", str(self.index), "--query-gpu=memory.used,utilization.gpu",
                              "--format=csv,noheader"], capture_output=True, text=True, timeout=10)
        return parse_smi(out.stdout.splitlines()[0])

    def _loop(self):
        while not self._stop.is_set():
            try:
                m, u = self._read()
                self.mem.append(m)
                self.util.append(u)
            except Exception:
                pass
            self._stop.wait(self.interval)

    def __enter__(self):
        if self.available:
            try:
                self.baseline = self._read()[0]
            except Exception:
                self.available = False
        if self.available:
            self._t = threading.Thread(target=self._loop, daemon=True)
            self._t.start()
        return self

    def __exit__(self, *exc):
        if self._t:
            self._stop.set()
            self._t.join(timeout=5)

    def summary(self) -> dict:
        if not self.mem:
            return {"peak_vram_mib": np.nan, "mean_gpu_util": np.nan}
        return {"peak_vram_mib": max(0.0, max(self.mem) - self.baseline), "mean_gpu_util": float(np.mean(self.util))}


def _vina_one(args):
    rec, lig, box, out, seed, exhaustiveness = args
    p = subprocess.run([str(require_tool("vina")), "--receptor", str(rec), "--ligand", str(lig),
                        "--center_x", str(box.center[0]), "--center_y", str(box.center[1]),
                        "--center_z", str(box.center[2]), "--size_x", str(box.size[0]), "--size_y", str(box.size[1]),
                        "--size_z", str(box.size[2]), "--cpu", "1", "--exhaustiveness", str(exhaustiveness),
                        "--seed", str(seed), "--out", str(out)], capture_output=True, text=True)
    return p.returncode == 0 and Path(out).exists()


def run_vina(rec: engines.Receptor, ligands, workdir: Path, threads: int, exhaustiveness: int = 8, seed: int = 42) -> int:
    workdir.mkdir(parents=True, exist_ok=True)
    jobs = []
    for l in ligands:
        lp = workdir / f"{l.name}.pdbqt"
        lp.write_text(l.pdbqt)
        jobs.append((rec.pdbqt, lp, rec.box, workdir / f"{l.name}.out.pdbqt", seed, exhaustiveness))
    with ThreadPoolExecutor(max_workers=threads) as ex:
        return sum(ex.map(_vina_one, jobs))


def run_case(case: str, rec: engines.Receptor, ligands, workdir: Path, threads: int, seed: int = 42,
             nrun: int = 20, exhaustiveness: int = 8) -> dict:
    """Time one case once. Returns wall seconds, the number of ligands that produced output, GPU stats."""
    shutil.rmtree(workdir, ignore_errors=True)
    workdir.mkdir(parents=True)
    with GpuSampler() as g:
        t0 = time.perf_counter()
        if case == "autodock-gpu":
            df, _ = engines.run_autodock_gpu(rec, ligands, workdir, nrun=nrun, seed=seed)
            n_ok = int(df["ad_ok"].sum())
        elif case.startswith("gnina"):
            kw = {"gnina (CNN, GPU)": {}, "gnina (CNN, CPU)": {"no_gpu": True},
                  "gnina (Vina, CPU)": {"cnn_scoring": "none"}}[case]
            engines.run_gnina(rec, ligands, workdir, exhaustiveness=exhaustiveness, num_modes=3, seed=seed,
                              cpu=threads, **kw)
            n_ok = len({k for f in workdir.glob("gnina_*.out.sdf") for k in engines.read_gnina_sdf(f)})
        elif case == "vina (CPU)":
            n_ok = run_vina(rec, ligands, workdir, threads, exhaustiveness, seed)
        else:
            raise ValueError(case)
        wall = time.perf_counter() - t0
    return {"case": case, "n": len(ligands), "n_ok": int(n_ok), "threads": threads, "wall_s": wall, **g.summary()}


def run_suite(rec, ligands, workdir: Path, plan: list[tuple[str, int, int]], repeats: int = 3, log=print,
              **kw) -> pd.DataFrame:
    """``plan`` is a list of (case, n_ligands, threads). A warm-up run precedes each GPU case's first timing."""
    rows, warmed = [], set()
    for case, n, threads in plan:
        subset = list(ligands)[:n]
        if case in ("autodock-gpu", "gnina (CNN, GPU)") and case not in warmed:
            run_case(case, rec, subset[:2], Path(workdir) / "warmup", threads, **kw)
            warmed.add(case)
        for rep in range(repeats):
            r = run_case(case, rec, subset, Path(workdir) / f"{case.split()[0]}_{n}_{threads}_{rep}", threads, **kw)
            r["repeat"] = rep
            rows.append(r)
            log(f"{case}: n={n} threads={threads} rep={rep} {r['wall_s']:.1f}s ({r['n_ok']} ok)")
    return pd.DataFrame(rows)


def summarize(runs: pd.DataFrame, reference: str = "vina (CPU)") -> pd.DataFrame:
    """Mean and SD of wall time over repeats, throughput, and speed-up over the reference case at the same n."""
    g = runs.groupby(["case", "n", "threads"], sort=False)
    t = g.agg(wall_mean_s=("wall_s", "mean"), wall_sd_s=("wall_s", lambda v: statistics.stdev(v) if len(v) > 1 else np.nan),
              n_ok=("n_ok", "min"), peak_vram_mib=("peak_vram_mib", "max"), mean_gpu_util=("mean_gpu_util", "mean"),
              repeats=("repeat", "count")).reset_index()
    t["s_per_compound"] = t["wall_mean_s"] / t["n"]
    t["compounds_per_hour"] = 3600.0 / t["s_per_compound"]
    ref = t[t["case"] == reference]
    if len(ref):
        best = ref.groupby("n")["s_per_compound"].min()                 # the reference at its best thread count
        t["speedup_vs_ref"] = [best.get(n, np.nan) / s for n, s in zip(t["n"], t["s_per_compound"])]
    return t


def projected_hours(s_per_compound: float, n_library: int) -> float:
    return s_per_compound * n_library / 3600.0
