"""Locate external docking executables and run them with a correct environment.

Executables are looked up in ``$VSBENCH_BIN``, then the active conda/mamba environment,
then a few conventional environment locations under the user's home, then ``$PATH``.
Nothing here is specific to one machine: set ``VSBENCH_BIN`` to the directory holding the
binaries and every lookup follows.

gnina is a static CUDA binary that needs its own environment's ``lib/`` on
``LD_LIBRARY_PATH`` (cuDNN, cuBLAS), so every subprocess is launched with that directory
prepended.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

# Conventional environment prefixes, relative to the user's home. Tried in order, after
# $VSBENCH_BIN and $CONDA_PREFIX, so an explicit setting always wins.
ENV_CANDIDATES = (
    "micromamba/envs/dock",
    "mambaforge/envs/dock",
    "miniforge3/envs/dock",
    "miniconda3/envs/dock",
    "anaconda3/envs/dock",
)

#: Executables without which a screening run cannot proceed, and what each is used for.
REQUIRED_TOOLS = {
    "autodock_gpu": "AutoDock-GPU: GPU docking with the AutoDock4 scoring function",
    "autogrid4": "AutoGrid4: precomputes the AutoDock4 affinity maps",
    "gnina": "gnina: Vina-type docking with CNN rescoring (the second scoring function)",
    "mk_prepare_receptor.py": "Meeko: receptor preparation (PDBQT and box definition)",
    "obabel": "Open Babel: file-format conversion and 3D embedding",
}

#: Executables only some stages need. Their absence disables that stage, not the run.
OPTIONAL_TOOLS = {
    "fpocket": "fpocket: binding-site detection (site, grid-strategy and druggability analyses)",
    "vina": "AutoDock Vina: the CPU baseline in the performance study",
    "prank": "P2Rank: local machine-learned pocket ranking, the offline alternative to DoGSiteScorer",
    "java": "Java runtime: required by P2Rank and by nothing else",
}

#: Conda/mamba environments that conventionally hold a JRE, tried after $JAVA_HOME.
JAVA_ENV_CANDIDATES = ("java17", "java", "openjdk", "jdk")


def find_java() -> Path | None:
    """A Java runtime for P2Rank. ``$JAVA_HOME`` wins, then a conventional conda env."""
    if os.environ.get("JAVA_HOME"):
        p = Path(os.environ["JAVA_HOME"]) / "bin" / "java"
        if p.is_file() and os.access(p, os.X_OK):
            return p
    roots = [Path(os.environ["MAMBA_ROOT_PREFIX"])] if os.environ.get("MAMBA_ROOT_PREFIX") else []
    home = Path.home()
    roots += [home / r for r in ("micromamba", "mambaforge", "miniforge3", "miniconda3", "anaconda3")]
    for root in roots:
        for env in JAVA_ENV_CANDIDATES:
            p = root / "envs" / env / "bin" / "java"
            if p.is_file() and os.access(p, os.X_OK):
                return p
    found = shutil.which("java")
    return Path(found) if found else None


def find_p2rank() -> Path | None:
    """The P2Rank launcher. ``$VSBENCH_P2RANK`` wins, then an unpacked copy under ``tools/``.

    P2Rank ships as a tarball holding a ``prank`` shell script beside its models, so it is
    normally not on ``PATH`` and cannot be found by name alone.
    """
    if os.environ.get("VSBENCH_P2RANK"):
        p = Path(os.environ["VSBENCH_P2RANK"])
        p = p / "prank" if p.is_dir() else p
        if p.is_file() and os.access(p, os.X_OK):
            return p
    here = Path(__file__).resolve()
    for base in (*here.parents[:4], Path.cwd()):
        for cand in sorted((base / "tools").glob("p2rank*/prank"), reverse=True):
            if cand.is_file() and os.access(cand, os.X_OK):
                return cand
    found = shutil.which("prank")
    return Path(found) if found else None


def _search_dirs() -> list[Path]:
    dirs = []
    if os.environ.get("VSBENCH_BIN"):
        dirs.append(Path(os.environ["VSBENCH_BIN"]))
    if os.environ.get("CONDA_PREFIX"):
        dirs.append(Path(os.environ["CONDA_PREFIX"]) / "bin")
    if os.environ.get("MAMBA_ROOT_PREFIX"):
        dirs.append(Path(os.environ["MAMBA_ROOT_PREFIX"]) / "envs" / "dock" / "bin")
    home = Path.home()
    dirs.extend(home / c / "bin" for c in ENV_CANDIDATES)
    return dirs


def find_tool(name: str) -> Path | None:
    if name == "java":                      # lives in its own environment, not beside the docking binaries
        return find_java()
    if name == "prank":                     # a script beside its models, normally not on PATH
        return find_p2rank()
    for d in _search_dirs():
        p = d / name
        if p.is_file() and os.access(p, os.X_OK):
            return p
    found = shutil.which(name)
    return Path(found) if found else None


def require_tool(name: str) -> Path:
    p = find_tool(name)
    if p is None:
        raise FileNotFoundError(
            f"required executable {name!r} not found. Set VSBENCH_BIN to the directory "
            f"containing it, or put it on PATH; run 'vsbench doctor' to see what is missing")
    return p


def tool_env(tool: Path) -> dict[str, str]:
    env = dict(os.environ)
    lib = tool.resolve().parent.parent / "lib"
    if lib.is_dir():
        env["LD_LIBRARY_PATH"] = f"{lib}:{env.get('LD_LIBRARY_PATH', '')}".rstrip(":")
    env["PATH"] = f"{tool.parent}:{env.get('PATH', '')}"
    return env


def run_tool(name: str, args: list[str], cwd: Path | None = None, timeout: float | None = None,
             check: bool = False) -> subprocess.CompletedProcess:
    exe = require_tool(name)
    proc = subprocess.run([str(exe), *map(str, args)], cwd=cwd, env=tool_env(exe),
                          capture_output=True, text=True, timeout=timeout)
    if check and proc.returncode != 0:
        raise RuntimeError(f"{name} failed ({proc.returncode}): {(proc.stderr or proc.stdout)[-600:]}")
    return proc


# --------------------------------------------------------------------- environment report
@dataclass(frozen=True)
class ToolStatus:
    name: str
    purpose: str
    required: bool
    path: Path | None

    @property
    def ok(self) -> bool:
        return self.path is not None


def check_tools() -> list[ToolStatus]:
    """Where each external executable resolves, required ones first.

    Reported rather than raised, so that a new user can see every missing piece at once
    instead of discovering them one failed stage at a time.
    """
    out = [ToolStatus(n, p, True, find_tool(n)) for n, p in REQUIRED_TOOLS.items()]
    out += [ToolStatus(n, p, False, find_tool(n)) for n, p in OPTIONAL_TOOLS.items()]
    return out


def check_python_packages() -> list[tuple[str, str | None, bool, str]]:
    """``(import name, version or None, required, purpose)`` for each Python dependency."""
    import importlib.metadata as md

    wanted = [("numpy", True, "arrays and metrics"),
              ("scipy", True, "nearest-neighbour property matching, statistics"),
              ("pandas", True, "every result table"),
              ("rdkit", True, "descriptors, fingerprints, 3D embedding"),
              ("matplotlib", True, "figures"),
              ("joblib", True, "parallel ligand preparation"),
              ("meeko", True, "receptor and ligand PDBQT preparation"),
              ("plip", False, "interaction fingerprints (extra 'interactions')"),
              ("admet_ai", False, "ADMET profiling (extra 'admet')"),
              ("sklearn", False, "used by the test suite only (extra 'dev')")]
    dist = {"rdkit": "rdkit", "admet_ai": "admet-ai", "sklearn": "scikit-learn"}
    rows = []
    for name, required, purpose in wanted:
        try:
            version = md.version(dist.get(name, name))
        except Exception:                       # not installed, or installed without metadata
            version = None
        rows.append((name, version, required, purpose))
    return rows


def gpu_report() -> str:
    """One line describing the visible NVIDIA GPU, or why none was found."""
    exe = shutil.which("nvidia-smi")
    if exe is None:
        return "nvidia-smi not found: AutoDock-GPU and gnina need an NVIDIA GPU"
    try:
        proc = subprocess.run([exe, "--query-gpu=name,memory.total,driver_version",
                               "--format=csv,noheader"], capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError) as exc:
        return f"nvidia-smi failed: {exc}"
    line = proc.stdout.strip().splitlines()
    return line[0].strip() if proc.returncode == 0 and line else (proc.stderr.strip() or "no GPU reported")
