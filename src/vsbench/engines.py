"""Receptor/grid preparation and the two docking engines (AutoDock-GPU, gnina).

Both engines dock into the same prepared receptor and box. They are independent
scoring functions: AutoDock4 empirical free energy (AutoDock-GPU) and Vina-type
minimisation with CNN rescoring (gnina).
"""
from __future__ import annotations

import hashlib
import json
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
from rdkit import Chem

from .box import Box
from .ligprep import Ligand
from .metals import ALL_IONS, Ion, ion_to_pdb_line
from .poses import split_models
from .tools import run_tool

_SAFE_NAME = re.compile(r"^[A-Za-z0-9_.-]+$")
_HIE_TIE_RE = re.compile(r"for residue_key='([^']+)', \d+ have passed.*tied for fewest")
_CLUSTER_RE = re.compile(
    r"^\s*(\d+)\s*\|\s*([+-]?\d+\.\d+)\s*\|\s*(\d+)\s*\|\s*([+-]?\d+\.\d+)\s*\|\s*(\d+)\s*\|", re.M)
_RUN_RE = re.compile(r"USER\s+Run\s*=\s*(\d+)")
_MAX_TIE_RETRIES = 25


def check_names(names) -> None:
    bad = [n for n in names if not _SAFE_NAME.match(n)]
    if bad:
        raise ValueError(f"unsafe ligand names for file output: {bad[:3]}")


# ------------------------------------------------------------------ receptor
_VENDOR = Path(__file__).parent / "vendor" / "ad4zn"
_WATERS = {"HOH", "WAT", "DOD"}
# AutoDock4Zn pairwise potentials for tetrahedral zinc coordination (Santos-Martins
# et al. 2014, Table 2), appended to the grid parameter file as in the upstream script.
_AD4ZN_NBP_LINES = (
    "nbp_r_eps 0.25 23.2135 12 6 NA TZ\n"
    "nbp_r_eps 2.1   3.8453 12 6 OA Zn\n"
    "nbp_r_eps 2.25  7.5914 12 6 SA Zn\n"
    "nbp_r_eps 1.0   0.0    12 6 HD Zn\n"
    "nbp_r_eps 2.0   0.2966 12 6  N Zn\n"
    "nbp_r_eps 2.0   0.0060 12 6 NA Zn\n"
)


@dataclass(frozen=True)
class PrepReport:
    """What receptor preparation kept and what it removed, so nothing is dropped silently."""
    n_protein_atoms: int
    waters_removed: int
    hetero_removed: dict
    ions_kept: tuple
    autogrid_attempts: int = 1       # >1 records that AutoGrid4 failed once and succeeded on an identical retry


@dataclass(frozen=True)
class Receptor:
    workdir: Path
    pdbqt: Path
    fld: Path
    box: Box
    ions: tuple = ()
    zn_mode: bool = False
    report: PrepReport | None = None


def source_to_pdb(src: Path, out_pdb: Path, keep_ions: Sequence[Ion] = ()) -> PrepReport:
    """Protein (mol2 or PDB) to a clean PDB: protein atoms, waters and other hetero groups
    removed, and exactly ``keep_ions`` added as HETATM records. Returns a PrepReport."""
    src = Path(src)
    if src.suffix.lower() == ".mol2":
        tmp = out_pdb.with_suffix(".all.pdb")
        run_tool("obabel", [str(src), "-O", str(tmp)], check=True)
        lines = tmp.read_text().splitlines()
    else:
        lines = src.read_text().splitlines()
    atoms, waters, hetero, last_serial = [], 0, {}, 0
    for line in lines:
        if line.startswith("ATOM"):
            if len(line) > 16 and line[16] not in (" ", "A"):
                continue                      # alternate locations other than A
            atoms.append(line)
            last_serial = max(last_serial, int(line[6:11]))
        elif line.startswith("TER"):
            atoms.append(line)
        elif line.startswith("HETATM"):
            res = line[17:20].strip()
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            is_kept_ion = any(np.allclose(xyz, i.xyz, atol=0.01) for i in keep_ions)
            if res in _WATERS:
                waters += 1
            elif not is_kept_ion:                 # a retained ion is re-written below, not removed
                hetero[res] = hetero.get(res, 0) + 1
    out = list(atoms)
    kept = []
    for n, ion in enumerate(keep_ions, start=1):
        out.append(ion_to_pdb_line(ion, last_serial + n))
        kept.append(ion.key())
    out.append("END")
    out_pdb.write_text("\n".join(out) + "\n")
    # Ions that were in the source but not requested are reported as removed hetero groups.
    return PrepReport(sum(1 for l in atoms if l.startswith("ATOM")), waters, hetero, tuple(kept))


def _restore_ion_type_case(original_pdbqt: Path, corrupted_pdbqt: Path) -> None:
    """zinc_pseudo.py upper-cases every AutoDock type (Zn -> ZN); restore the original
    casing by atom serial, which the script preserves for pre-existing atoms."""
    orig = {int(l[6:11]): l[77:79] for l in original_pdbqt.read_text().splitlines()
            if l.startswith(("ATOM", "HETATM")) and l[17:20].strip().upper() in ALL_IONS | {"FE2", "CU1"}}
    fixed = []
    for l in corrupted_pdbqt.read_text().splitlines(keepends=True):
        if l.startswith(("ATOM", "HETATM")):
            t = orig.get(int(l[6:11]))
            if t is not None and l[77:79] != t:
                l = l[:77] + t + l[79:]
        fixed.append(l)
    corrupted_pdbqt.write_text("".join(fixed))


def _add_zinc_pseudo_atoms(receptor_pdbqt: Path, out_pdbqt: Path, log: Path) -> None:
    proc = subprocess.run([sys.executable, str(_VENDOR / "zinc_pseudo.py"), "-r", str(receptor_pdbqt),
                           "-o", str(out_pdbqt)], capture_output=True, text=True)
    log.write_text(proc.stdout + proc.stderr)
    if proc.returncode != 0 or not out_pdbqt.exists():
        raise RuntimeError(f"zinc_pseudo.py failed (see {log}); is there a Zn in {receptor_pdbqt}?")
    _restore_ion_type_case(receptor_pdbqt, out_pdbqt)


def _zn_augment_gpf(gpf: Path, receptor_pdbqt: Path, workdir: Path) -> None:
    """Rewrite the grid parameter file for AutoDock4Zn: TZ receptor, AD4Zn.dat, coordination terms."""
    tz = workdir / "receptor_tz.pdbqt"
    _add_zinc_pseudo_atoms(receptor_pdbqt, tz, workdir / "zinc_pseudo.log")
    dat = workdir / "AD4Zn.dat"
    shutil.copyfile(_VENDOR / "AD4Zn.dat", dat)
    out, saw = [], False
    for line in gpf.read_text().splitlines(keepends=True):
        if line.startswith("parameter_file "):
            out.append(f"parameter_file {dat.name}\n")
            saw = True
        elif line.startswith("receptor_types "):
            out.append(line.rstrip("\n") + " TZ\n")
        elif line.startswith("receptor "):
            out.append(f"receptor {tz.name}\n")
        else:
            out.append(line)
    if not saw:
        out.insert(0, f"parameter_file {dat.name}\n")
    out.append(_AD4ZN_NBP_LINES)
    gpf.write_text("".join(out))


def _verify_ions_in_pdbqt(pdbqt: Path, ions: Sequence[Ion]) -> None:
    """Raise if a requested ion did not survive preparation (no silent apo receptors)."""
    present = {}
    for l in pdbqt.read_text().splitlines():
        if l.startswith(("ATOM", "HETATM")):
            res = l[17:20].strip().upper()
            res = {"FE2": "FE", "CU1": "CU"}.get(res, res)
            present[res] = present.get(res, 0) + 1
    want = {}
    for ion in ions:
        want[ion.element] = want.get(ion.element, 0) + 1
    missing = {el: n for el, n in want.items() if present.get(el, 0) < n}
    if missing:
        raise RuntimeError(f"ions requested but missing from the prepared receptor: {missing}")


def prepare_receptor(src: Path, workdir: Path, box: Box, keep_ions: Sequence[Ion] = (),
                     zn_mode: bool | None = None) -> Receptor:
    """Protein file -> receptor PDBQT + AutoDock4 grid maps, retaining ``keep_ions``.

    With a zinc ion present (``zn_mode`` None means auto) AutoDock-GPU maps are built with
    AutoDock4Zn; gnina uses the plain receptor PDBQT, which carries the same ions.
    Cached on the source bytes, box, ions and mode.
    """
    keep_ions = tuple(keep_ions)
    zn = any(i.element == "ZN" for i in keep_ions) if zn_mode is None else bool(zn_mode)
    if zn and not any(i.element == "ZN" for i in keep_ions):
        raise ValueError("zn_mode requested but no zinc ion was supplied")
    workdir.mkdir(parents=True, exist_ok=True)
    stamp = workdir / "stamp.json"
    fld = workdir / "receptor.maps.fld"
    key = json.dumps({"src": hashlib.sha256(Path(src).read_bytes()).hexdigest(), "center": box.center,
                      "size": box.size, "spacing": box.spacing, "zn": zn,
                      "ions": [(i.element, i.xyz) for i in keep_ions]})
    report_path = workdir / "prep_report.json"
    if (stamp.exists() and stamp.read_text() == key and fld.exists() and (workdir / "receptor.pdbqt").exists()
            and report_path.exists()):
        rep = PrepReport(**{**json.loads(report_path.read_text()),
                            "ions_kept": tuple(json.loads(report_path.read_text())["ions_kept"])})
        return Receptor(workdir, workdir / "receptor.pdbqt", fld, box, tuple(i.key() for i in keep_ions), zn, rep)

    pdb = workdir / "receptor.pdb"
    report = source_to_pdb(Path(src), pdb, keep_ions)
    forced: dict[str, str] = {}
    for _ in range(_MAX_TIE_RETRIES):
        args = ["--read_pdb", str(pdb), "-o", "receptor", "-p", "-a", "--default_altloc", "A",
                "--box_center", *map(str, box.center), "--box_size", *map(str, box.size),
                "-g", "receptor.gpf"]
        if forced:
            args += ["--set_template", ",".join(f"{k}={v}" for k, v in forced.items())]
        proc = run_tool("mk_prepare_receptor.py", args, cwd=workdir)
        (workdir / "prep.log").write_text(proc.stdout + proc.stderr)
        if proc.returncode == 0:
            break
        m = _HIE_TIE_RE.search(proc.stderr)
        if m and m.group(1) not in forced:   # histidine tautomer tie: resolve as HIE, as the module does
            forced[m.group(1)] = "HIE"
            continue
        raise RuntimeError(f"mk_prepare_receptor.py failed: {proc.stderr[-500:]}")
    else:
        raise RuntimeError("mk_prepare_receptor.py kept failing on tautomer ties")
    _verify_ions_in_pdbqt(workdir / "receptor.pdbqt", keep_ions)
    if zn:
        _zn_augment_gpf(workdir / "receptor.gpf", workdir / "receptor.pdbqt", workdir)
    # AutoGrid4 once failed on a valid PDBQT ("ATOM/HETATM line bad x,y,z ... line 1") that it then read
    # correctly on 64 of 64 identical re-runs, so one identical retry is allowed. A second failure is
    # deterministic and raises; the number of attempts is kept in the preparation report.
    errors = []
    for attempt in (1, 2):
        proc = run_tool("autogrid4", ["-p", "receptor.gpf", "-l", "receptor.glg"], cwd=workdir)
        if proc.returncode == 0 and fld.exists():
            break
        errors.append((proc.stdout + proc.stderr)[-300:])
    else:
        raise RuntimeError("autogrid4 failed twice: " + " | ".join(errors))
    report = replace(report, autogrid_attempts=len(errors) + 1)
    if zn:
        # TZ is a receptor-only pseudo-type: it has no map file of its own and acts through the
        # other maps, so verify the pseudo-atoms and that autogrid actually loaded AD4Zn.
        tz = [l for l in (workdir / "receptor_tz.pdbqt").read_text().splitlines()
              if l.startswith(("ATOM", "HETATM")) and l[77:79].strip() == "TZ"]
        glg = (workdir / "receptor.glg").read_text()
        if not tz or "AD4Zn.dat" not in glg or "rec_index(TZ)" not in glg:
            raise RuntimeError("AutoDock4Zn requested but TZ pseudo-atoms or the AD4Zn parameters were not applied")
    report_path.write_text(json.dumps({**report.__dict__, "ions_kept": list(report.ions_kept)}))
    stamp.write_text(key)
    return Receptor(workdir, workdir / "receptor.pdbqt", fld, box, tuple(i.key() for i in keep_ions), zn, report)


# ----------------------------------------------------------------- AutoDock-GPU
def parse_clusters(dlg_text: str) -> list[dict]:
    """Rows of the DLG clustering histogram, best cluster first."""
    return [{"rank": int(r), "energy": float(e), "run": int(run), "mean_energy": float(me),
             "size": int(n)} for r, e, run, me, n in _CLUSTER_RE.findall(dlg_text)]


def summarize_dlg(dlg_text: str) -> dict:
    """Lowest-energy pose and largest-cluster pose of one docking."""
    clusters = parse_clusters(dlg_text)
    if not clusters:
        raise ValueError("no clustering histogram in DLG")
    best = min(clusters, key=lambda c: c["energy"])
    largest = max(clusters, key=lambda c: (c["size"], -c["energy"]))
    return {"energy": best["energy"], "run": best["run"], "n_clusters": len(clusters),
            "largest_size": largest["size"], "largest_energy": largest["energy"],
            "largest_run": largest["run"]}


def dlg_run_blocks(dlg_text: str) -> dict[int, list[str]]:
    """DOCKED lines of each run, keyed by run number."""
    out = {}
    for block in split_models(dlg_text):
        for line in block:
            m = _RUN_RE.search(line)
            if m:
                out[int(m.group(1))] = block
                break
    return out


def dlg_complete(path: Path) -> bool:
    if not path.exists():
        return False
    tail = path.read_text(errors="ignore")
    return "CLUSTERING HISTOGRAM" in tail and "RANKING" in tail


def run_autodock_gpu(receptor: Receptor, ligands: Sequence[Ligand], workdir: Path, nrun: int = 20,
                     seed: int = 42, autostop: bool = True, heuristics: bool = True,
                     chunk: int = 200, log=None) -> tuple[pd.DataFrame, dict]:
    """Batch-dock ligands; resumable (finished DLGs are kept and not redone).

    Every ligand PDBQT lives in ``workdir/ligands`` because AutoDock-GPU writes
    each result log beside its input, so a separate directory per settings
    combination is what keeps runs from overwriting each other.
    """
    ok = [l for l in ligands if l.ok]
    check_names(l.name for l in ok)
    lig_dir = workdir / "ligands"
    lig_dir.mkdir(parents=True, exist_ok=True)
    for l in ok:
        p = lig_dir / f"{l.name}.pdbqt"
        if not p.exists():
            p.write_text(l.pdbqt)
    todo = [l.name for l in ok if not dlg_complete(lig_dir / f"{l.name}.dlg")]
    t0 = time.monotonic()
    for i in range(0, len(todo), chunk):
        part = todo[i:i + chunk]
        listing = workdir / f"filelist_{i // chunk:05d}.txt"
        listing.write_text(str(receptor.fld.resolve()) + "\n"
                           + "\n".join(str((lig_dir / f"{n}.pdbqt").resolve()) for n in part) + "\n")
        proc = run_tool("autodock_gpu", [
            "--filelist", str(listing), "--nrun", str(nrun), "--autostop", str(int(autostop)),
            "--heuristics", str(int(heuristics)), "--seed", f"{seed},{seed},{seed}",
            "--xmloutput", "0"])
        (workdir / f"autodock_gpu_{i // chunk:05d}.log").write_text(proc.stdout + proc.stderr)
        if log:
            log(f"autodock-gpu chunk {i // chunk}: {len(part)} ligands, exit {proc.returncode}")
    wall = time.monotonic() - t0
    rows = []
    for l in ligands:
        row = {"name": l.name, "ad_ok": False, "ad_energy": np.nan, "ad_run": np.nan,
               "ad_n_clusters": np.nan, "ad_largest_size": np.nan, "ad_largest_energy": np.nan,
               "ad_largest_run": np.nan}
        dlg = lig_dir / f"{l.name}.dlg"
        if l.ok and dlg_complete(dlg):
            s = summarize_dlg(dlg.read_text())
            row.update(ad_ok=True, ad_energy=s["energy"], ad_run=s["run"],
                       ad_n_clusters=s["n_clusters"], ad_largest_size=s["largest_size"],
                       ad_largest_energy=s["largest_energy"], ad_largest_run=s["largest_run"])
        rows.append(row)
    info = {"engine": "autodock-gpu", "n_requested": len(ligands), "n_docked_now": len(todo),
            "wall_s": wall, "nrun": nrun, "autostop": autostop, "heuristics": heuristics,
            "seed": seed}
    return pd.DataFrame(rows), info


# ----------------------------------------------------------------------- gnina
_GNINA_FIELDS = {"minimizedAffinity": "gn_vina", "CNNscore": "gn_cnn_score",
                 "CNNaffinity": "gn_cnn_affinity", "CNN_VS": "gn_cnn_vs"}


def read_gnina_sdf(path: Path) -> dict[str, list[Chem.Mol]]:
    """Poses per ligand name, in gnina's output order (top-ranked first)."""
    out: dict[str, list[Chem.Mol]] = {}
    for mol in Chem.SDMolSupplier(str(path), removeHs=False, sanitize=False):
        if mol is None:
            continue
        out.setdefault(mol.GetProp("_Name"), []).append(mol)
    return out


def gnina_top_pose_row(name: str, poses: list[Chem.Mol] | None) -> dict:
    row = {"name": name, "gn_ok": False, **{v: np.nan for v in _GNINA_FIELDS.values()}}
    if poses:
        top = poses[0]
        if all(k in top.GetPropsAsDict() for k in _GNINA_FIELDS):
            row.update(gn_ok=True, **{v: float(top.GetProp(k)) for k, v in _GNINA_FIELDS.items()})
    return row


def run_gnina(receptor: Receptor, ligands: Sequence[Ligand], workdir: Path, exhaustiveness: int = 8,
              num_modes: int = 3, seed: int = 42, cpu: int = 16, device: int = 0,
              chunk: int = 250, log=None, cnn_scoring: str = "rescore",
              no_gpu: bool = False) -> tuple[pd.DataFrame, dict]:
    """Batch-dock via multi-molecule SDF; resumable per chunk. ``cnn_scoring='none'`` and ``no_gpu`` exist
    for the performance study (CPU-only baselines); the benchmark itself always uses rescoring on the GPU."""
    ok = [l for l in ligands if l.ok]
    check_names(l.name for l in ok)
    workdir.mkdir(parents=True, exist_ok=True)
    box = receptor.box
    t0 = time.monotonic()
    done_now = 0
    top: dict[str, list[Chem.Mol]] = {}
    for i in range(0, len(ok), chunk):
        part = ok[i:i + chunk]
        stem = workdir / f"gnina_{i // chunk:05d}"
        sdf_out, marker = stem.with_suffix(".out.sdf"), stem.with_suffix(".done")
        if not marker.exists():
            sdf_in = stem.with_suffix(".in.sdf")
            w = Chem.SDWriter(str(sdf_in))
            for l in part:
                m = Chem.MolFromMolBlock(l.molblock, removeHs=False)
                m.SetProp("_Name", l.name)
                w.write(m)
            w.close()
            proc = run_tool("gnina", [
                "-r", str(receptor.pdbqt), "-l", str(sdf_in),
                "--center_x", box.center[0], "--center_y", box.center[1], "--center_z", box.center[2],
                "--size_x", box.size[0], "--size_y", box.size[1], "--size_z", box.size[2],
                "--exhaustiveness", exhaustiveness, "--num_modes", num_modes, "--seed", seed,
                "--cpu", cpu, "--device", device, "--cnn_scoring", cnn_scoring,
                *(["--no_gpu"] if no_gpu else []), "-o", sdf_out])
            stem.with_suffix(".log").write_text(proc.stdout + proc.stderr)
            if proc.returncode == 0 and sdf_out.exists():
                marker.write_text("ok")
                done_now += len(part)
            if log:
                log(f"gnina chunk {i // chunk}: {len(part)} ligands, exit {proc.returncode}")
        if sdf_out.exists():
            cached = read_gnina_sdf(sdf_out)
            # The .done marker is keyed by chunk index alone, never by the ligands in the chunk, so a
            # directory reused for a different ligand list serves that other list's poses: the names
            # asked for come back absent and are recorded as docking failures rather than as a stale
            # cache. Re-docking the seed-variance subset after its composition changed is precisely
            # that case. A name present that was never requested proves it; a merely missing name can
            # be an honest per-ligand failure, so only the former is reported.
            foreign = set(cached) - {l.name for l in part}
            if foreign and log:
                log(f"  WARNING: {stem.name}.out.sdf holds {len(foreign)} compound(s) absent from this "
                    f"chunk, so the cached poses belong to a different ligand list. Delete {workdir} "
                    f"and re-dock instead of trusting them.")
            for name, poses in cached.items():
                top[name] = poses
    wall = time.monotonic() - t0
    rows = [gnina_top_pose_row(l.name, top.get(l.name) if l.ok else None) for l in ligands]
    info = {"engine": "gnina", "n_requested": len(ligands), "n_docked_now": done_now, "wall_s": wall,
            "exhaustiveness": exhaustiveness, "num_modes": num_modes, "seed": seed, "cpu": cpu}
    return pd.DataFrame(rows), info


def clean_workdir(path: Path) -> None:
    shutil.rmtree(path, ignore_errors=True)
