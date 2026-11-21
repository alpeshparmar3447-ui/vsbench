"""Metal-ion awareness: find ions, judge whether they matter for a pocket, check coordination.

Why this exists: receptor preparation that keeps only protein atoms silently
turns a metalloenzyme into its apo form, and the run still exits cleanly. Every
template is therefore audited, ions are retained or restored on purpose, and
chelating ligands are recognised. Policy follows the project rules: zinc gets
AutoDock4Zn (tetrahedral pseudo-atoms) in AutoDock-GPU; other ions are kept as
plain AutoDock4/Vina non-bonded atoms, which carry no coordination-specific
parameters, and results for such pockets are flagged.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from rdkit import Chem
from scipy.spatial import cKDTree

# Divalent or transition-metal ions that are catalytic or structural.
RETAIN_ELEMENTS = frozenset({"ZN", "MG", "MN", "FE", "CA", "CU", "NI", "CO"})
# Monovalent ions are usually crystallisation buffer, but can be functional (K+ in
# pyruvate kinase). They are reported always and kept only when close to the ligand.
MONOVALENT = frozenset({"NA", "K"})
_ALIASES = {"FE2": "FE", "CU1": "CU", "ZN2": "ZN", "MG2": "MG"}
ALL_IONS = RETAIN_ELEMENTS | MONOVALENT
DONOR_ELEMENTS = frozenset({"N", "O", "S"})

DEFAULT_POCKET_CUTOFF = 8.0       # A from any ligand heavy atom: the ion can shape binding
DEFAULT_MONOVALENT_CUTOFF = 6.0


@dataclass(frozen=True)
class Ion:
    element: str
    resname: str
    chain: str
    resnum: int
    xyz: tuple[float, float, float]
    source: str

    def key(self) -> str:
        return f"{self.element}{self.resnum}{self.chain}"


def _norm(sym: str) -> str:
    s = sym.strip().upper()
    return _ALIASES.get(s, s)


# ------------------------------------------------------------------- readers
def read_pdb_ions(text: str) -> list[Ion]:
    """Metal ions from HETATM records of a PDB file (first alternate location only)."""
    out = []
    for line in text.splitlines():
        if not line.startswith("HETATM") or len(line) < 54:
            continue
        if line[16] not in (" ", "A"):
            continue
        res = line[17:20].strip()
        elem = line[76:78].strip() if len(line) >= 78 else ""
        sym = _norm(elem or res)
        if sym not in ALL_IONS or _norm(res) not in ALL_IONS:
            continue
        out.append(Ion(sym, res, line[21].strip(), int(line[22:26]),
                       (float(line[30:38]), float(line[38:46]), float(line[46:54])), "pdb"))
    return out


def read_mol2_ions(text: str) -> list[Ion]:
    out = []
    in_atoms = False
    for line in text.splitlines():
        if line.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            continue
        if line.startswith("@<TRIPOS>") and in_atoms:
            break
        if in_atoms and line.strip():
            t = line.split()
            sym = _norm(t[5].split(".")[0])
            if sym in ALL_IONS and len(t) >= 6 and not t[5].upper().startswith(("C.", "N.", "O.", "S.")):
                out.append(Ion(sym, t[7] if len(t) > 7 else sym, "", 0,
                               (float(t[2]), float(t[3]), float(t[4])), "mol2"))
    return out


def _read_cif_atoms(text: str):
    cols, rows, inloop = [], [], False
    for line in text.splitlines():
        if line.startswith("_atom_site."):
            cols.append(line.strip().split(".")[1])
            inloop = True
            continue
        if inloop:
            if line.startswith(("#", "loop_", "_")):
                if rows:
                    break
                continue
            if line.strip():
                rows.append(line.split())
    ix = {c: i for i, c in enumerate(cols)}
    return ix, rows


def read_cif_ions(text: str) -> list[Ion]:
    ix, rows = _read_cif_atoms(text)
    out = []
    for r in rows:
        if r[ix["group_PDB"]] != "HETATM" or r[ix["pdbx_PDB_model_num"]] != "1":
            continue
        if r[ix["label_alt_id"]] not in (".", "A"):
            continue
        sym = _norm(r[ix["type_symbol"]])
        if sym in ALL_IONS and _norm(r[ix["label_comp_id"]]) in ALL_IONS:
            out.append(Ion(sym, r[ix["label_comp_id"]], r[ix["auth_asym_id"]],
                           int(r[ix["auth_seq_id"]]),
                           (float(r[ix["Cartn_x"]]), float(r[ix["Cartn_y"]]), float(r[ix["Cartn_z"]])), "cif"))
    return out


def read_structure_ions(path: Path) -> list[Ion]:
    text = Path(path).read_text(errors="ignore")
    suffix = Path(path).suffix.lower()
    if suffix == ".mol2":
        return read_mol2_ions(text)
    if suffix in (".cif", ".mmcif"):
        return read_cif_ions(text)
    return read_pdb_ions(text)


def read_ca_coordinates(path: Path) -> np.ndarray:
    """C-alpha positions of a PDB, mmCIF or mol2 protein (used to test the frame)."""
    text = Path(path).read_text(errors="ignore")
    suffix = Path(path).suffix.lower()
    xyz = []
    if suffix == ".mol2":
        in_atoms = False
        for line in text.splitlines():
            if line.startswith("@<TRIPOS>ATOM"):
                in_atoms = True
                continue
            if line.startswith("@<TRIPOS>") and in_atoms:
                break
            if in_atoms and line.split()[1:2] == ["CA"] and line.split()[5].startswith("C."):
                t = line.split()
                xyz.append((float(t[2]), float(t[3]), float(t[4])))
    elif suffix in (".cif", ".mmcif"):
        ix, rows = _read_cif_atoms(text)
        for r in rows:
            if (r[ix["group_PDB"]] == "ATOM" and r[ix["label_atom_id"]] == "CA"
                    and r[ix["pdbx_PDB_model_num"]] == "1" and r[ix["label_alt_id"]] in (".", "A")):
                xyz.append((float(r[ix["Cartn_x"]]), float(r[ix["Cartn_y"]]), float(r[ix["Cartn_z"]])))
    else:
        for line in text.splitlines():
            if line.startswith("ATOM") and line[12:16].strip() == "CA" and line[16] in (" ", "A"):
                xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return np.asarray(xyz, dtype=float).reshape(-1, 3)


# ------------------------------------------------------------------ geometry
def frame_overlap(ca_a: np.ndarray, ca_b: np.ndarray, tol: float = 0.5) -> float:
    """Fraction of points in ``ca_a`` that have a point of ``ca_b`` within ``tol`` angstrom.

    Near 1.0 means the two files share a coordinate frame, so ions read from the
    original deposition can be transplanted into the benchmark template.
    """
    if len(ca_a) == 0 or len(ca_b) == 0:
        return 0.0
    d, _ = cKDTree(ca_b).query(ca_a)
    return float(np.mean(d <= tol))


def min_distance(ion: Ion, ligand_xyz: np.ndarray) -> float:
    return float(np.min(np.linalg.norm(np.asarray(ligand_xyz) - np.asarray(ion.xyz), axis=1)))


def ions_near(ions: list[Ion], ligand_xyz: np.ndarray, cutoff: float) -> list[tuple[Ion, float]]:
    """Ions within ``cutoff`` of any ligand heavy atom, nearest first."""
    hits = [(i, min_distance(i, ligand_xyz)) for i in ions]
    return sorted([h for h in hits if h[1] <= cutoff], key=lambda h: h[1])


def in_box(xyz, center, size, margin: float = 2.0) -> bool:
    """Is a point inside the search box, enlarged by ``margin`` on every side?"""
    return bool(all(abs(p - c) <= s / 2.0 + margin for p, c, s in zip(xyz, center, size)))


def relevant_ions(ions: list[Ion], ligand_xyz: np.ndarray, pocket_cutoff: float = DEFAULT_POCKET_CUTOFF,
                  monovalent_cutoff: float = DEFAULT_MONOVALENT_CUTOFF, box=None,
                  box_margin: float = 2.0) -> list[tuple[Ion, float]]:
    """Ions to keep in the receptor.

    Functional ions are kept when within ``pocket_cutoff`` of the ligand or, if a search
    ``box`` is given, anywhere inside it (they then shape the grid maps). Monovalent
    ions are kept only when very close to the ligand. Ions elsewhere in the protein
    are reported by the audit and not kept. Result is nearest-first.
    """
    out = []
    for ion in ions:
        d = min_distance(ion, ligand_xyz)
        inside = box is not None and in_box(ion.xyz, box.center, box.size, box_margin)
        if ion.element in RETAIN_ELEMENTS and (d <= pocket_cutoff or inside):
            out.append((ion, d))
        elif ion.element in MONOVALENT and d <= monovalent_cutoff:
            out.append((ion, d))
    return sorted(out, key=lambda h: h[1])


def fetch_deposition(pdb_id: str, cache_dir: Path, timeout: float = 60.0,
                     offline: bool = False) -> Path:
    """Original deposited structure (PDB format, else mmCIF), cached on disk.

    The cache is checked first, so a populated cache makes this function offline in practice.
    With ``offline=True`` a cache miss raises instead of reaching the network, and the message
    names the file to drop in, so a run can be guaranteed to make no outbound request.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("pdb", "cif"):
        dst = cache_dir / f"{pdb_id.upper()}.{ext}"
        if dst.exists() and dst.stat().st_size > 1000:
            return dst
    if offline:
        raise FileNotFoundError(
            f"{pdb_id}: no cached deposition and offline mode is on. "
            f"Place {pdb_id.upper()}.pdb (or .cif) in {cache_dir} to use it.")

    import urllib.error
    import urllib.request
    for ext in ("pdb", "cif"):
        dst = cache_dir / f"{pdb_id.upper()}.{ext}"
        try:
            with urllib.request.urlopen(
                    f"https://files.rcsb.org/download/{pdb_id.upper()}.{ext}", timeout=timeout) as r:
                dst.write_bytes(r.read())
            return dst
        except urllib.error.HTTPError:
            continue
    raise FileNotFoundError(f"no deposition found for {pdb_id}")


def coordination(ion_xyz, atoms, cutoff: float = 2.8) -> list[tuple[str, str, int, float]]:
    """Donor atoms (N, O, S) within ``cutoff`` of an ion.

    ``atoms`` is an iterable of (atom_name, resname, resnum, element, (x, y, z)).
    Returns (atom_name, resname, resnum, distance) sorted by distance.
    """
    out = []
    p = np.asarray(ion_xyz)
    for name, res, num, el, xyz in atoms:
        if el.upper() in DONOR_ELEMENTS:
            d = float(np.linalg.norm(np.asarray(xyz) - p))
            if d <= cutoff:
                out.append((name, res, num, d))
    return sorted(out, key=lambda t: t[3])


def read_pdb_protein_atoms(text: str, include_hetero: bool = False):
    """(name, resname, resnum, element, xyz) for protein ATOM records.

    Hetero atoms (waters, ligands) are excluded unless ``include_hetero`` is set, so a
    coordination count means protein donors by default; ions are never returned.
    """
    for line in text.splitlines():
        if not line.startswith(("ATOM", "HETATM")) or len(line) < 54 or line[16] not in (" ", "A"):
            continue
        res = line[17:20].strip()
        if line.startswith("HETATM") and (not include_hetero or _norm(res) in ALL_IONS):
            continue
        el = line[76:78].strip() if len(line) >= 78 else line[12:16].strip()[:1]
        yield (line[12:16].strip(), res, int(line[22:26]), el or line[12:16].strip()[:1],
               (float(line[30:38]), float(line[38:46]), float(line[46:54])))


# ------------------------------------------------------------------ PDB output
def ion_to_pdb_line(ion: Ion, serial: int) -> str:
    """A PDB HETATM record for an ion, in the layout meeko and AutoDock tools expect."""
    name = ion.element.title() if len(ion.element) > 1 else ion.element
    atom_name = f"{ion.element:>2}" if len(ion.element) == 2 else f" {ion.element}"
    return (f"HETATM{serial:5d} {atom_name:<4s} {ion.element:>3s} {(ion.chain or 'Z')[:1]}{ion.resnum % 10000:4d}    "
            f"{ion.xyz[0]:8.3f}{ion.xyz[1]:8.3f}{ion.xyz[2]:8.3f}  1.00  0.00          {ion.element:>2s}  ")


# --------------------------------------------------------------------- ligands
# Metal-chelating groups that must be deprotonated to bind a zinc ion in AutoDock4Zn.
# SMARTS and target indices follow the project's tested docking module
# (_METAL_CHELATOR_PATTERNS in insilico_dock/stages/prep_ligands.py).
CHELATOR_PATTERNS: dict[str, tuple[str, int]] = {
    "hydroxamic_acid": ("[CX3](=O)[#7X3][OX2H1]", -1),
    "catechol": ("[c;R]([OX2H1])[c;R][OX2H1]", 1),
    "hydroxypyridinone_or_tropolone": ("[OX2H1][c;R][c;R]=[OX1]", 0),
}
# Chelating groups that are recognised but NOT auto-ionised: they need a hand-checked
# input state, so a ligand carrying one is flagged instead of silently mis-handled.
MANUAL_IONISATION_PATTERNS: dict[str, str] = {
    "thiol": "[SX2H1]",
    "phosphonate_or_phosphate_OH": "[PX4](=O)[OX2H1]",
    "8_hydroxyquinoline": "[OX2H1]c1cccc2cccnc12",
    "primary_sulfonamide": "[SX4](=O)(=O)[NX3H2]",
}


def chelator_classes(mol: Chem.Mol) -> dict[str, list[str]]:
    """Which chelator groups a molecule carries: {'auto': [...], 'manual': [...]}."""
    auto = [n for n, (s, _) in CHELATOR_PATTERNS.items() if mol.HasSubstructMatch(Chem.MolFromSmarts(s))]
    manual = [n for n, s in MANUAL_IONISATION_PATTERNS.items() if mol.HasSubstructMatch(Chem.MolFromSmarts(s))]
    return {"auto": auto, "manual": manual}


def deprotonate_chelators(mol: Chem.Mol) -> tuple[Chem.Mol, list[str]]:
    """Deprotonate the coordinating oxygen of each recognised chelator group."""
    idx: dict[int, str] = {}
    for name, (smarts, pos) in CHELATOR_PATTERNS.items():
        for match in mol.GetSubstructMatches(Chem.MolFromSmarts(smarts)):
            idx.setdefault(match[pos], name)
    if not idx:
        return mol, []
    rw = Chem.RWMol(mol)
    for i in idx:
        atom = rw.GetAtomWithIdx(i)
        atom.SetFormalCharge(atom.GetFormalCharge() - 1)
        atom.SetNumExplicitHs(max(atom.GetTotalNumHs() - 1, 0))
        atom.SetNoImplicit(True)
    out = rw.GetMol()
    Chem.SanitizeMol(out)
    return out, sorted(set(idx.values()))


_ION_LINE = re.compile(r"^HETATM")


# ------------------------------------------------------------------ planning
@dataclass(frozen=True)
class MetalPlan:
    """Which ions a template's receptor must carry, and why."""
    ions: tuple[Ion, ...]
    source: str                 # 'none' | 'benchmark_file' | 'restored' | 'mixed'
    audited: bool               # the original deposition was checked and shares the file's frame
    zn: bool                    # a zinc ion is retained: AutoDock4Zn and chelator ionisation apply
    n_deposition_ions: int
    remote: tuple[str, ...]     # functional ions seen but not retained, as 'ZN378A@29.0'
    note: str

    @property
    def label(self) -> str:
        return ";".join(i.key() for i in self.ions)


def plan_metals(protein_file: Path, ligand_xyz, box, pdb_id: str | None = None, cache_dir: Path | None = None,
                mode: str = "auto", pocket_cutoff: float = DEFAULT_POCKET_CUTOFF,
                monovalent_cutoff: float = DEFAULT_MONOVALENT_CUTOFF, restore_from_pdb: bool = True,
                min_frame_overlap: float = 0.95, fetch=fetch_deposition) -> MetalPlan:
    """Decide the ions to retain for one template.

    Ions already in ``protein_file`` are used as they are. If ``restore_from_pdb`` and a PDB
    code is given, the original deposition is fetched, its C-alpha atoms are compared with the
    file's, and only when the frames agree are its ions added (a benchmark file that stripped
    them would otherwise yield an apo receptor). Functional ions within ``pocket_cutoff`` of
    the ligand or inside the search box are kept; others are listed as remote, not dropped
    silently. Every outcome carries a human-readable note.
    """
    if mode == "off":
        return MetalPlan((), "none", False, False, 0, (), "metal handling disabled by configuration")
    ligand_xyz = np.asarray(ligand_xyz, dtype=float)
    file_ions = list(read_structure_ions(protein_file))
    candidates = list(file_ions)
    notes, audited, n_dep = [], False, 0
    if restore_from_pdb and pdb_id and cache_dir is not None:
        try:
            dep = fetch(pdb_id, Path(cache_dir))
            overlap = frame_overlap(read_ca_coordinates(protein_file), read_ca_coordinates(dep))
            if overlap >= min_frame_overlap:
                audited = True
                dep_ions = read_structure_ions(dep)
                n_dep = len(dep_ions)
                for ion in dep_ions:
                    if not any(np.linalg.norm(np.array(ion.xyz) - np.array(c.xyz)) < 0.5 for c in candidates):
                        candidates.append(ion)
            else:
                notes.append(f"deposition {pdb_id} not used: frame overlap {overlap:.2f} < {min_frame_overlap}")
        except Exception as exc:
            notes.append(f"deposition {pdb_id} unavailable ({type(exc).__name__}); ions audited from the file only")
    kept = relevant_ions(candidates, ligand_xyz, pocket_cutoff, monovalent_cutoff, box)
    kept_keys = {(i.element, i.xyz) for i, _ in kept}
    remote = []
    for ion in candidates:
        if (ion.element, ion.xyz) not in kept_keys and ion.element in RETAIN_ELEMENTS:
            d = min_distance(ion, ligand_xyz)
            if d <= 30.0:
                remote.append(f"{ion.key()}@{d:.1f}")
    ions = tuple(i for i, _ in kept)
    n_restored = sum(1 for i in ions if not any(
        np.linalg.norm(np.array(i.xyz) - np.array(f.xyz)) < 0.5 for f in file_ions))
    if not ions:
        source = "none"
    elif n_restored == 0:
        source = "benchmark_file"
    elif n_restored == len(ions):
        source = "restored"
    else:
        source = "mixed"
    if ions:
        notes.append("retained: " + ", ".join(f"{i.key()} ({d:.1f} A from ligand)" for i, d in kept))
    if remote:
        notes.append("functional ions kept out (outside pocket and search box): " + ", ".join(remote))
    if not ions and not remote:
        notes.append("no functional ions near the ligand" + (" (deposition checked)" if audited else ""))
    zn = any(i.element == "ZN" for i in ions)
    return MetalPlan(ions, source, audited, zn, n_dep, tuple(remote), "; ".join(notes))
