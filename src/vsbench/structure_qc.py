"""Structure assessment of a receptor template (plan section 4.3).

Reads the benchmark protein file (mol2) and, when available, the original deposition and
reports what a careful reviewer of the receptor would ask: chains, chain breaks and
numbering gaps, missing backbone atoms, a Ramachandran flag, resolution, engineered
residues, alternate locations and any cofactor or other hetero group next to the ligand.

The Ramachandran regions are coarse ellipses around the three populated areas of the
general-residue plot. They are a screening flag, not an outlier call: compare templates
with each other, never read one fraction in isolation.
"""
from __future__ import annotations

import math
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

from . import metals as K

PEPTIDE_BOND_MAX = 2.0           # A: C(i)-N(i+1) beyond this is a chain break
BACKBONE = ("N", "CA", "C", "O")
# (name, centre phi/psi in degrees, core half-widths); 'allowed' doubles the half-widths.
RAMA_REGIONS = [
    ("alpha_R", (-65.0, -40.0), (35.0, 30.0)),
    ("beta", (-110.0, 130.0), (45.0, 35.0)),
    ("alpha_L", (55.0, 30.0), (27.5, 27.5)),
]
WATERS = {"HOH", "WAT", "DOD"}


@dataclass(frozen=True)
class Residue:
    chain: str
    number: int
    name: str
    atoms: dict


_SUBST_NAME = re.compile(r"^([A-Za-z]{1,3})(-?\d+)[A-Za-z]?$")


def substructure_residue(t: list[str]) -> tuple[str, int, str] | None:
    """Chain, residue number and residue name from one mol2 SUBSTRUCTURE line.

    Two layouts occur in this deposition. The common one ends ``... PRO A 206``, with the chain and the
    residue number in separate fields, so the number is the twelfth. A structure whose residue numbers
    reach four digits writes them joined to the chain instead, ``... GLU B1385``, which is one field
    shorter; in this panel that is the mTOR kinase domain, residues 1385 to 2549. Reading the number
    from the trailing field, and falling back to the subst_name, covers both layouts without depending
    on the field count.
    """
    if len(t) < 7:
        return None
    chain = t[5] if t[5] and not t[5].startswith("*") else "A"
    name = t[6] if not t[6].startswith("*") else (t[9] if len(t) > 9 else "UNK")   # '****': unusual residue
    if len(t) >= 12:
        try:
            return chain, int(t[11]), name
        except ValueError:
            pass
    for field_ in (t[-1], t[1]):                                        # 'B1385', then 'GLU1385'
        m = _SUBST_NAME.match(field_)
        if m:
            return chain, int(m.group(2)), name if name != "UNK" else m.group(1).upper()
    return None


def read_mol2_residues(path: Path) -> list[Residue]:
    """Residues in file order, with chain and original numbering from the SUBSTRUCTURE block."""
    text = Path(path).read_text(errors="ignore")
    atoms_by_sub: dict[int, dict] = {}
    section = None
    subs: dict[int, tuple[str, int, str]] = {}
    for line in text.splitlines():
        if line.startswith("@<TRIPOS>"):
            section = line.strip()
            continue
        t = line.split()
        if section == "@<TRIPOS>ATOM" and len(t) >= 8:
            sid = int(t[6])
            if not t[5].startswith("H"):
                atoms_by_sub.setdefault(sid, {})[t[1]] = np.array([float(t[2]), float(t[3]), float(t[4])])
        elif section == "@<TRIPOS>SUBSTRUCTURE":
            res = substructure_residue(t)
            if res is not None:
                subs[int(t[0])] = res                       # chain, original number, residue name
    out = []
    for sid in sorted(atoms_by_sub):
        chain, num, name = subs.get(sid, ("?", sid, "UNK"))
        out.append(Residue(chain, num, name, atoms_by_sub[sid]))
    return out


def dihedral(p0, p1, p2, p3) -> float:
    """Dihedral angle in degrees (IUPAC sign convention)."""
    b0, b1, b2 = (np.asarray(p0) - p1), (np.asarray(p2) - p1), (np.asarray(p3) - p2)   # b0 points from p1 back to p0
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return float(math.degrees(math.atan2(np.dot(np.cross(b1n, v), w), np.dot(v, w))))


def classify_rama(phi: float, psi: float) -> str:
    """'core', 'allowed' or 'outside' under the coarse regions above."""
    def inside(scale):
        for _, (c_phi, c_psi), (a, b) in RAMA_REGIONS:
            d_phi = ((phi - c_phi + 180.0) % 360.0) - 180.0
            d_psi = ((psi - c_psi + 180.0) % 360.0) - 180.0
            if (d_phi / (a * scale)) ** 2 + (d_psi / (b * scale)) ** 2 <= 1.0:
                return True
        return False
    return "core" if inside(1.0) else "allowed" if inside(2.0) else "outside"


NON_POLYMER = WATERS | K.ALL_IONS | {"FE2", "CU1"}


def backbone_report(residues: list[Residue]) -> dict:
    """Backbone geometry of the polypeptide only.

    Waters and free ions are dropped first. They are not part of the chain, and a crystallographic
    water numbered in the thousands otherwise reads as a numbering gap at the end of the protein, or,
    where its residue name was lost, as a whole extra chain.
    """
    residues = [r for r in residues if r.name.upper() not in NON_POLYMER]
    n_breaks = n_gaps = n_unexplained = missing = 0
    classes = {"core": 0, "allowed": 0, "outside": 0}
    chains = sorted({r.chain for r in residues})
    for i, r in enumerate(residues):
        terminal = i == len(residues) - 1 or residues[i + 1].chain != r.chain
        # Chain termini are excluded: their carboxylate oxygens are named differently (OXT) in many files.
        if not terminal and not all(a in r.atoms for a in BACKBONE):
            missing += 1
    for prev, cur in zip(residues, residues[1:]):
        if prev.chain != cur.chain:
            continue
        broken = "C" in prev.atoms and "N" in cur.atoms and np.linalg.norm(prev.atoms["C"] - cur.atoms["N"]) > PEPTIDE_BOND_MAX
        jump = cur.number - prev.number > 1
        n_breaks += bool(broken)
        n_gaps += bool(jump)
        n_unexplained += bool(broken and not jump)
    for i in range(1, len(residues) - 1):
        a, r, b = residues[i - 1], residues[i], residues[i + 1]
        if r.name in ("GLY", "PRO") or not (a.chain == r.chain == b.chain):
            continue
        try:
            joined = (np.linalg.norm(a.atoms["C"] - r.atoms["N"]) <= PEPTIDE_BOND_MAX
                      and np.linalg.norm(r.atoms["C"] - b.atoms["N"]) <= PEPTIDE_BOND_MAX)
            if not joined:
                continue
            phi = dihedral(a.atoms["C"], r.atoms["N"], r.atoms["CA"], r.atoms["C"])
            psi = dihedral(r.atoms["N"], r.atoms["CA"], r.atoms["C"], b.atoms["N"])
        except KeyError:
            continue
        classes[classify_rama(phi, psi)] += 1
    n_eval = sum(classes.values())
    frac = lambda k: (classes[k] / n_eval) if n_eval else float("nan")
    return {"n_chains": len(chains), "n_residues": len(residues),
            "n_chain_breaks": n_breaks, "n_numbering_gaps": n_gaps, "n_unexplained_breaks": n_unexplained,
            "n_missing_backbone": missing, "n_rama_evaluated": n_eval,
            "rama_core": frac("core"), "rama_allowed": frac("allowed"), "rama_outside": frac("outside")}


# ------------------------------------------------------------------ deposition
def parse_resolution(text: str, is_cif: bool) -> float:
    if is_cif:
        m = re.search(r"^_refine\.ls_d_res_high\s+([0-9.]+)", text, re.M) or \
            re.search(r"^_reflns\.d_resolution_high\s+([0-9.]+)", text, re.M)
    else:
        m = re.search(r"REMARK\s+2\s+RESOLUTION\.\s+([0-9.]+)", text)
    return float(m.group(1)) if m else float("nan")


def count_engineered(text: str, is_cif: bool) -> dict:
    """Engineered residues recorded by the depositors: expression tags, mutations, conflicts."""
    tag = mut = other = 0
    if is_cif:
        for m in re.finditer(r"^\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(.*)$", text, re.M):
            pass  # mmCIF engineered-residue tables are free-form; counted from SEQADV-equivalent text below
        up = text.upper()
        tag, mut = up.count("EXPRESSION TAG"), up.count("ENGINEERED MUTATION")
        return {"expression_tag": tag, "engineered_mutation": mut, "other_seqadv": 0}
    for line in text.splitlines():
        if line.startswith("SEQADV"):
            comment = line[49:].upper()
            if "EXPRESSION TAG" in comment:
                tag += 1
            elif "ENGINEERED" in comment or "MUTATION" in comment:
                mut += 1
            else:
                other += 1
    return {"expression_tag": tag, "engineered_mutation": mut, "other_seqadv": other}


def count_altlocs(text: str, is_cif: bool) -> int:
    if is_cif:
        return 0
    return sum(1 for l in text.splitlines() if l.startswith(("ATOM", "HETATM")) and len(l) > 16 and l[16] not in (" ",))


def hetero_near_ligand(deposition: Path, ligand_xyz: np.ndarray, cutoff: float = 6.0) -> list[str]:
    """Residue names of non-water, non-ion hetero groups within ``cutoff`` of the ligand, excluding
    the ligand itself (its atoms coincide with the template ligand)."""
    text = Path(deposition).read_text(errors="ignore")
    groups: dict[tuple, list] = {}
    if Path(deposition).suffix.lower() == ".pdb":
        for l in text.splitlines():
            if l.startswith("HETATM") and l[17:20].strip() not in WATERS | K.ALL_IONS | {"FE2", "CU1"}:
                key = (l[17:20].strip(), l[21], int(l[22:26]))
                groups.setdefault(key, []).append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    else:
        ix, rows = K._read_cif_atoms(text)
        for r in rows:
            if r[ix["group_PDB"]] == "HETATM" and r[ix["pdbx_PDB_model_num"]] == "1":
                res = r[ix["label_comp_id"]]
                if res not in WATERS | K.ALL_IONS | {"FE2", "CU1"}:
                    key = (res, r[ix["auth_asym_id"]], int(r[ix["auth_seq_id"]]))
                    groups.setdefault(key, []).append(tuple(float(r[ix[c]]) for c in ("Cartn_x", "Cartn_y", "Cartn_z")))
    found = []
    lig = np.asarray(ligand_xyz)
    for (res, _, _), pts in groups.items():
        pts = np.asarray(pts)
        d = np.linalg.norm(pts[:, None, :] - lig[None, :, :], axis=2)
        if np.mean(d.min(axis=1) < 0.4) > 0.8:      # this group is the template ligand itself
            continue
        if d.min() <= cutoff:
            found.append(res)
    return sorted(set(found))


@dataclass
class TemplateQC:
    template: str
    n_chains: int = 0
    n_residues: int = 0
    n_chain_breaks: int = 0
    n_numbering_gaps: int = 0
    n_unexplained_breaks: int = 0
    n_missing_backbone: int = 0
    n_rama_evaluated: int = 0
    rama_core: float = float("nan")
    rama_allowed: float = float("nan")
    rama_outside: float = float("nan")
    resolution: float = float("nan")
    expression_tag: int = 0
    engineered_mutation: int = 0
    other_seqadv: int = 0
    altloc_atoms: int = 0
    hetero_near_ligand: str = ""
    flags: str = ""
    extra: dict = field(default_factory=dict)

    def as_dict(self) -> dict:
        d = asdict(self)
        d.pop("extra")
        return d


def qc_template(template: str, protein_mol2: Path, ligand_xyz=None, deposition: Path | None = None) -> TemplateQC:
    qc = TemplateQC(template=template)
    for k, v in backbone_report(read_mol2_residues(protein_mol2)).items():
        setattr(qc, k, v)
    if deposition is not None and Path(deposition).exists():
        text = Path(deposition).read_text(errors="ignore")
        is_cif = Path(deposition).suffix.lower() == ".cif"
        qc.resolution = parse_resolution(text, is_cif)
        for k, v in count_engineered(text, is_cif).items():
            setattr(qc, k, v)
        qc.altloc_atoms = count_altlocs(text, is_cif)
        if ligand_xyz is not None:
            qc.hetero_near_ligand = ";".join(hetero_near_ligand(deposition, np.asarray(ligand_xyz)))
    flags = []
    if qc.n_unexplained_breaks:
        flags.append("unexplained_chain_break")
    if qc.n_missing_backbone:
        flags.append("missing_backbone_atoms")
    if qc.n_chains > 1:
        flags.append("multiple_chains")
    if qc.hetero_near_ligand:
        flags.append("hetero_group_near_ligand")
    if qc.engineered_mutation:
        flags.append("engineered_mutation")
    if qc.resolution == qc.resolution and qc.resolution > 3.0:
        flags.append("resolution_gt_3A")
    qc.flags = ";".join(flags)
    return qc


def flag_rama_outliers(outside_fractions: dict[str, float], k: float = 3.0) -> list[str]:
    """Templates whose outside-region fraction exceeds median + k robust SDs of the set."""
    vals = np.array([v for v in outside_fractions.values() if v == v])
    if len(vals) < 4:
        return []
    med = float(np.median(vals))
    mad = float(np.median(np.abs(vals - med))) * 1.4826
    cut = med + k * max(mad, 0.005)
    return sorted(t for t, v in outside_fractions.items() if v == v and v > cut)
