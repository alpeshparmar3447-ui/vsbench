"""Binding-site identification and druggability assessment (plan section 4.4).

Three independent pocket detectors are run on each template and judged by one question that
needs no assumption about druggability: does the detector recover the co-crystal site?
* fpocket    - geometry (alpha spheres), local
* P2Rank     - machine-learned ligandability, local (needs Java)
* DoGSiteScorer - difference-of-Gaussians pockets with a druggability score, public web API

Every detector's pockets are reduced to a set of (chain, residue number) and compared with
the residues that line the co-crystal ligand, so the methods are directly comparable. Control
structures (a validated druggable pocket and a polar, difficult one) calibrate how far each
score can be trusted, following Michel et al. (Front Chem 2020).
"""
from __future__ import annotations

import csv
import io
import json
import re
import subprocess
import time
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .structure_qc import read_mol2_residues, substructure_residue

SITE_CUTOFF = 4.5          # A: protein residues with an atom this close to the ligand line the site
HIT_RECOVERY = 0.5         # a pocket "recovers" the site if it holds at least half of the lining residues
DOGSITE_URL = "https://proteins.plus/api/dogsite_rest"


@dataclass
class Pocket:
    rank: int
    score: float                      # the detector's own ranking score
    druggability: float               # fpocket druggability score / P2Rank probability / DoGSite drugScore
    volume: float
    residues: frozenset = frozenset()
    extra: dict = field(default_factory=dict)


# ---------------------------------------------------------------------- native site
def native_site(residues, ligand_xyz, cutoff: float = SITE_CUTOFF) -> frozenset:
    """(chain, residue number) of every residue with a heavy atom within ``cutoff`` of the ligand."""
    lig = np.asarray(ligand_xyz)
    out = set()
    for r in residues:
        pts = np.array([v for k, v in r.atoms.items()])
        if len(pts) and np.min(np.linalg.norm(pts[:, None, :] - lig[None, :, :], axis=2)) <= cutoff:
            out.add((r.chain, r.number))
    return frozenset(out)


def write_pdb_from_mol2(mol2: Path, out_pdb: Path) -> int:
    """Heavy-atom PDB of a mol2 protein with the ORIGINAL chain and residue numbers, so every
    detector (including DoGSite, which reads the deposition) reports the same residue identity."""
    section, subs, atoms = None, {}, []
    for line in Path(mol2).read_text(errors="ignore").splitlines():
        if line.startswith("@<TRIPOS>"):
            section = line.strip()
            continue
        t = line.split()
        if section == "@<TRIPOS>ATOM" and len(t) >= 8 and not t[5].startswith("H"):
            atoms.append((t[1], t[5].split(".")[0], float(t[2]), float(t[3]), float(t[4]), int(t[6])))
        elif section == "@<TRIPOS>SUBSTRUCTURE":
            res = substructure_residue(t)
            if res is not None:
                subs[int(t[0])] = res
    if atoms and not subs:
        # Every residue would be written as UNK. A pocket detector then reports a shattered surface
        # instead of the real site, and PLIP reports no interactions at all, with no error either way.
        raise ValueError(f"{Path(mol2).name}: no residue could be read from its SUBSTRUCTURE records, "
                         f"so chain and residue identity would be lost")
    lines = []
    for n, (name, el, x, y, z, sid) in enumerate(atoms, 1):
        chain, num, res = subs.get(sid, ("A", sid, "UNK"))
        an = f" {name:<3s}" if len(name) < 4 else name[:4]
        lines.append(f"ATOM  {n:5d} {an} {res[:3]:>3s} {(chain or 'A')[:1]}{num % 10000:4d}    "
                     f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el.upper():>2s}")
    Path(out_pdb).write_text("\n".join(lines) + "\nEND\n")
    return len(atoms)


def residues_from_pdb_text(text: str) -> frozenset:
    out = set()
    for l in text.splitlines():
        if l.startswith(("ATOM", "HETATM")) and len(l) >= 26:
            try:
                out.add((l[21].strip() or "A", int(l[22:26])))
            except ValueError:
                pass
    return frozenset(out)


def recovery(native: frozenset, pocket: frozenset) -> float:
    """Fraction of the native-site residues that a pocket contains."""
    return len(native & pocket) / len(native) if native else float("nan")


def site_recovery(pockets: list[Pocket], native: frozenset, hit: float = HIT_RECOVERY) -> dict:
    """How well a detector's ranked pockets recover the co-crystal site."""
    ordered = sorted(pockets, key=lambda p: p.rank)
    recs = [recovery(native, p.residues) for p in ordered]
    first_hit = next((p.rank for p, r in zip(ordered, recs) if r >= hit), None)
    best = int(np.argmax(recs)) if recs else None
    out = {"n_pockets": len(ordered), "best_recovery": float(max(recs)) if recs else float("nan"),
           "rank_of_site": first_hit if first_hit is not None else float("nan"),
           "top1_hit": bool(recs and recs[0] >= hit), "top3_hit": bool(any(r >= hit for r in recs[:3])),
           "site_pocket_druggability": float("nan"), "site_pocket_volume": float("nan")}
    if best is not None and recs[best] >= hit:
        out["site_pocket_druggability"] = ordered[best].druggability
        out["site_pocket_volume"] = ordered[best].volume
    return out


# ---------------------------------------------------------------------- fpocket
def parse_fpocket_info(text: str) -> list[dict]:
    pockets = []
    parts = re.split(r"^Pocket\s+(\d+)\s*:\s*$", text, flags=re.M)
    for i in range(1, len(parts), 2):
        d = {"pocket": int(parts[i])}
        for line in parts[i + 1].splitlines():
            m = re.match(r"\s*(.+?)\s*:\s*([-0-9.eE]+)\s*$", line)
            if m:
                d[m.group(1).strip()] = float(m.group(2))
        pockets.append(d)
    return pockets


def alpha_sphere_center(pqr_text: str):
    """Centroid of fpocket's alpha-sphere centres (the pocket's own geometric centre), or None."""
    pts = [(float(l[30:38]), float(l[38:46]), float(l[46:54])) for l in pqr_text.splitlines()
           if l.startswith(("ATOM", "HETATM")) and len(l) >= 54]
    return tuple(float(v) for v in np.mean(pts, axis=0)) if pts else None


def run_fpocket(pdb: Path, workdir: Path, binary: str = "fpocket") -> list[Pocket]:
    from .tools import run_tool
    workdir.mkdir(parents=True, exist_ok=True)
    local = workdir / Path(pdb).name
    local.write_text(Path(pdb).read_text())
    proc = run_tool(binary, ["-f", str(local)], cwd=workdir)
    out = workdir / f"{local.stem}_out"
    info = out / f"{local.stem}_info.txt"
    if proc.returncode != 0 or not info.exists():
        raise RuntimeError(f"fpocket failed: {(proc.stdout + proc.stderr)[-300:]}")
    pockets = []
    for d in parse_fpocket_info(info.read_text()):
        atm = out / "pockets" / f"pocket{d['pocket']}_atm.pdb"
        res = residues_from_pdb_text(atm.read_text()) if atm.exists() else frozenset()
        extra = {k: v for k, v in d.items() if k != "pocket"}
        vert = out / "pockets" / f"pocket{d['pocket']}_vert.pqr"
        if vert.exists():
            extra["center"] = alpha_sphere_center(vert.read_text())
        pockets.append(Pocket(d["pocket"], d.get("Score", float("nan")), d.get("Druggability Score", float("nan")),
                              d.get("Volume", float("nan")), res, extra))
    return pockets


# ---------------------------------------------------------------------- P2Rank
def parse_p2rank_predictions(text: str) -> list[Pocket]:
    rows = list(csv.DictReader(io.StringIO(text), skipinitialspace=True))
    out = []
    for r in rows:
        r = {k.strip(): (v or "").strip() for k, v in r.items()}
        res = frozenset((tok.split("_")[0], int(re.match(r"-?\d+", tok.split("_")[1]).group(0)))
                        for tok in r["residue_ids"].split() if "_" in tok)
        extra = {"sas_points": float(r["sas_points"]), "surf_atoms": float(r["surf_atoms"])}
        # P2Rank reports the pocket centre directly, where fpocket needs it derived from the
        # alpha spheres. Both land in extra["center"], which is what box seeding reads.
        try:
            extra["center"] = (float(r["center_x"]), float(r["center_y"]), float(r["center_z"]))
        except (KeyError, TypeError, ValueError):
            pass
        out.append(Pocket(int(r["rank"]), float(r["score"]), float(r["probability"]),
                          float("nan"), res, extra))
    return out


def run_p2rank(pdb: Path, workdir: Path, prank: Path, java_bin: Path, threads: int = 4) -> list[Pocket]:
    import os
    workdir.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env["PATH"] = f"{Path(java_bin).parent}:{env['PATH']}"
    env["JAVA_HOME"] = str(Path(java_bin).parent.parent)
    proc = subprocess.run([str(prank), "predict", "-f", str(pdb), "-o", str(workdir), "-threads", str(threads)],
                          capture_output=True, text=True, env=env)
    csvs = list(workdir.glob("*_predictions.csv"))
    if proc.returncode != 0 or not csvs:
        raise RuntimeError(f"P2Rank failed: {(proc.stdout + proc.stderr)[-400:]}")
    return parse_p2rank_predictions(csvs[0].read_text())


# ---------------------------------------------------------------------- DoGSiteScorer
def _http(method: str, url: str, body: dict | None = None, timeout: float = 60.0):
    data = json.dumps(body).encode() if body is not None else None
    req = urllib.request.Request(url, data=data, method=method,
                                 headers={"Accept": "application/json", "Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.status, r.read().decode()
    except urllib.error.HTTPError as e:
        return e.code, e.read().decode()


def parse_dogsite_table(text: str) -> list[dict]:
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    return [r for r in rows if re.fullmatch(r"P_\d+", r["name"])]


def run_dogsite(pdb_code: str, cache: Path, http=_http, poll_s: float = 5.0, max_wait_s: float = 600.0,
                sleep=time.sleep, chain: str = "") -> list[Pocket]:
    """DoGSiteScorer pockets for a public PDB entry (whole-protein pocket detection, no ligand given).
    Results are cached as JSON so a repeated run never calls the server again."""
    cache.mkdir(parents=True, exist_ok=True)
    f = cache / f"{pdb_code.upper()}.json"
    if f.exists():
        return [Pocket(p["rank"], p["score"], p["druggability"], p["volume"], frozenset(map(tuple, p["residues"])), p["extra"])
                for p in json.loads(f.read_text())]
    body = {"dogsite": {"pdbCode": pdb_code.lower(), "analysisDetail": "1", "bindingSitePredictionGranularity": "1",
                        "ligand": "", "chain": chain}}
    waited = 0.0
    while True:
        code, txt = http("POST", DOGSITE_URL, body)
        if code in (200, 202):
            loc = json.loads(txt)["location"]
            break
        if code == 429 and waited < max_wait_s:
            sleep(30)
            waited += 30
            continue
        raise RuntimeError(f"DoGSite submission failed: HTTP {code} {txt[:200]}")
    while waited <= max_wait_s:
        code, txt = http("GET", loc)
        d = json.loads(txt)
        if d.get("status_code") == 200:
            break
        sleep(poll_s)
        waited += poll_s
    else:
        raise TimeoutError(f"DoGSite job for {pdb_code} did not finish in {max_wait_s:.0f}s")
    _, table = http("GET", d["result_table"])
    rows = parse_dogsite_table(table)
    res_urls = {Path(u).name[:-8]: u for u in d["residues"] if u.endswith("_res.pdb")}
    pockets = []
    for r in rows:
        name = r["name"]
        url = next((u for k, u in res_urls.items() if k.endswith(f"_{name}")), None)
        residues = residues_from_pdb_text(http("GET", url)[1]) if url else frozenset()
        pockets.append(Pocket(0, float(r.get("simpleScore", "nan")), float(r["drugScore"]), float(r["volume"]), residues,
                              {k: r[k] for k in ("enclosure", "depth", "surface") if k in r} | {"name": name}))
    pockets.sort(key=lambda p: (-p.druggability, -p.volume))
    for i, p in enumerate(pockets, 1):
        p.rank = i
    f.write_text(json.dumps([{"rank": p.rank, "score": p.score, "druggability": p.druggability, "volume": p.volume,
                              "residues": sorted(p.residues), "extra": p.extra} for p in pockets]))
    return pockets


# ---------------------------------------------------------------------- ligand burial
def ligand_buried_fraction(protein_pdb: Path, ligand_xyz) -> float:
    """1 - SASA(ligand in the complex) / SASA(ligand alone): how enclosed the co-crystal ligand is."""
    from Bio.PDB import PDBParser
    from Bio.PDB.SASA import ShrakeRupley
    lig = np.asarray(ligand_xyz)
    lines = Path(protein_pdb).read_text().splitlines()
    het = [f"HETATM{9000 + i:5d}  C   LIG Z   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C" for i, (x, y, z) in enumerate(lig)]
    body = [l for l in lines if l.startswith("ATOM")]

    def sasa(pdb_lines):
        p = Path(protein_pdb).with_suffix(".tmp_sasa.pdb")
        p.write_text("\n".join(pdb_lines) + "\nEND\n")
        s = PDBParser(QUIET=True).get_structure("x", str(p))
        ShrakeRupley(probe_radius=1.4, n_points=200).compute(s, level="A")
        v = sum(a.sasa for a in s.get_atoms() if a.get_parent().get_resname() == "LIG")
        p.unlink()
        return v
    alone, bound = sasa(het), sasa(body + het)
    return float(1.0 - bound / alone) if alone > 0 else float("nan")


# ---------------------------------------------------------------------- control structures
# Michel et al. 2020 calibrate pocket tools on a validated druggable target and a difficult one.
# Each control carries the residue name of a bound organic ligand that defines its native site.
CONTROLS = {
    "1HXB": ("positive", "ROC", "HIV-1 protease with saquinavir"),
    "1M17": ("positive", "AQ4", "EGFR kinase with erlotinib"),
    "1CIL": ("positive", "ETS", "carbonic anhydrase II with a sulfonamide (zinc enzyme)"),
    "1BZC": ("negative", "TPI", "PTP1B with an inhibitor: a polar, historically difficult pocket"),
}


def control_native_site(pdb_text: str, ligand_resname: str, cutoff: float = SITE_CUTOFF):
    """(site residues, ligand coordinates) for a control deposition."""
    lig, prot = [], []
    for l in pdb_text.splitlines():
        if l.startswith("HETATM") and l[17:20].strip() == ligand_resname and l[16] in (" ", "A"):
            lig.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
        elif l.startswith("ATOM") and l[16] in (" ", "A") and l[76:78].strip() != "H":
            prot.append(((l[21].strip() or "A", int(l[22:26])), (float(l[30:38]), float(l[38:46]), float(l[46:54]))))
    if not lig:
        raise ValueError(f"ligand {ligand_resname} not found")
    lig = np.array(lig)
    site = {k for k, xyz in prot if np.min(np.linalg.norm(lig - np.array(xyz), axis=1)) <= cutoff}
    return frozenset(site), lig


def heavy_protein_pdb(pdb_text: str, out: Path) -> Path:
    keep = [l for l in pdb_text.splitlines() if l.startswith("ATOM") and l[16] in (" ", "A") and l[76:78].strip() != "H"]
    out.write_text("\n".join(keep) + "\nEND\n")
    return out
