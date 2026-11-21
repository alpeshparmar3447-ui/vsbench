"""Metal-aware preparation and docking on two real metalloenzymes.

HDAC8 with the hydroxamate SAHA (PDB 1T69): a catalytic zinc, so AutoDock4Zn and chelator
ionisation apply. FEN1 template 5fv7 from LIT-PCBA: two catalytic Mg2+ ions that the benchmark
file lost, restored from the original deposition. The original deposition is downloaded once
and cached; tests skip if it is unavailable and the network is not.
"""
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from vsbench import data, engines, ligprep, metals as K, redock
from vsbench.box import box_from_coords
from vsbench.poses import pose_to_mol
from vsbench.rmsd import symmetry_rmsd
from vsbench.tools import find_tool

ROOT = Path(__file__).parents[1]
CACHE = ROOT / "data" / "pdb_cache"
SAHA = "ONC(=O)CCCCCCC(=O)Nc1ccccc1"

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(
        not all(find_tool(t) for t in ("autodock_gpu", "autogrid4", "gnina", "obabel", "obrms", "mk_prepare_receptor.py"))
        or shutil.which("nvidia-smi") is None, reason="docking executables or GPU not available"),
]


def _deposition(pdb_id):
    try:
        return K.fetch_deposition(pdb_id, CACHE)
    except Exception:
        pytest.skip(f"deposition {pdb_id} not cached and not downloadable")


@pytest.fixture(scope="module")
def hdac8(tmp_path_factory):
    dep = _deposition("1T69")
    text = dep.read_text()
    block = "\n".join(l for l in text.splitlines() if l.startswith("HETATM") and l[17:20].strip() == "SHH")
    xtal = AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles(SAHA), Chem.MolFromPDBBlock(block, proximityBonding=True))
    d = tmp_path_factory.mktemp("hdac8")
    sdf = d / "saha.sdf"
    Chem.MolToMolFile(xtal, str(sdf))
    heavy = Chem.RemoveHs(xtal)
    box = box_from_coords(heavy.GetConformer().GetPositions())
    ions = K.read_pdb_ions(text)
    return dict(pdb=dep, sdf=sdf, heavy=heavy, box=box, ions=ions, dir=d)


def test_zinc_is_found_in_the_catalytic_site_and_coordinated_by_the_expected_residues(hdac8):
    plan = redock.plan_for(hdac8["pdb"], hdac8["heavy"].GetConformer().GetPositions(), hdac8["box"], "1T69", None)
    assert plan.zn and [i.element for i in plan.ions] == ["ZN"] and plan.audited
    atoms = list(K.read_pdb_protein_atoms(hdac8["pdb"].read_text()))
    shell = {(r, n) for _, r, n, _ in K.coordination(plan.ions[0].xyz, atoms, cutoff=2.8)}
    assert {("HIS", 180), ("ASP", 178), ("ASP", 267)} <= shell        # the HDAC8 catalytic triad


def test_zinc_receptor_uses_autodock4zn_and_keeps_the_ion(hdac8):
    d = hdac8["dir"] / "rec_zn"
    rec = engines.prepare_receptor(hdac8["pdb"], d, hdac8["box"], keep_ions=hdac8["ions"])
    assert rec.zn_mode and rec.ions == ("ZN378A",) and rec.report.ions_kept == ("ZN378A",)
    assert "ZN" not in rec.report.hetero_removed and rec.report.hetero_removed == {"SHH": 19}
    tz = [l for l in (d / "receptor_tz.pdbqt").read_text().splitlines() if l[77:79].strip() == "TZ"]
    zn = [l for l in (d / "receptor_tz.pdbqt").read_text().splitlines() if l[77:79].strip() == "Zn"]
    assert len(tz) >= 1 and len(zn) == 1                                # pseudo-atom added, Zn casing restored
    dist = np.linalg.norm(np.array([float(tz[0][30:38]), float(tz[0][38:46]), float(tz[0][46:54])])
                          - np.array([float(zn[0][30:38]), float(zn[0][38:46]), float(zn[0][46:54])]))
    assert dist == pytest.approx(2.0, abs=0.15)                         # TZ sits 2 A from the zinc
    glg = (d / "receptor.glg").read_text()
    assert "AD4Zn.dat" in glg and "rec_index(TZ)" in glg
    assert "Zn" in (d / "receptor.pdbqt").read_text()                   # gnina's receptor carries the ion too


def test_an_apo_receptor_is_a_different_receptor(hdac8):
    zn_dir, apo_dir = hdac8["dir"] / "rec_zn", hdac8["dir"] / "rec_apo"
    rec = engines.prepare_receptor(hdac8["pdb"], zn_dir, hdac8["box"], keep_ions=hdac8["ions"])
    apo = engines.prepare_receptor(hdac8["pdb"], apo_dir, hdac8["box"], keep_ions=[])
    assert apo.ions == () and not apo.zn_mode
    assert "Zn" not in (apo_dir / "receptor.pdbqt").read_text()
    assert (zn_dir / "receptor.OA.map").read_bytes() != (apo_dir / "receptor.OA.map").read_bytes()


def test_chelator_is_deprotonated_only_for_a_zinc_pocket(hdac8):
    zn = ligprep.prepare_ligand(SAHA, "saha", 42, chelate_zn=True)
    other = ligprep.prepare_ligand(SAHA, "saha", 42, chelate_zn=False)
    assert zn.deprotonated and not other.deprotonated and zn.chelators_auto == other.chelators_auto == "hydroxamic_acid"
    charge = lambda L: sum(float(l[70:76]) for l in L.pdbqt.splitlines() if l.startswith(("ATOM", "HETATM")))
    assert charge(zn) == pytest.approx(-1.0, abs=0.02) and charge(other) == pytest.approx(0.0, abs=0.02)


def test_template_redock_records_the_zinc_state_and_the_failed_pose_is_not_called_validated(hdac8):
    res = redock.redock_template(hdac8["pdb"], hdac8["sdf"], "HDAC8", hdac8["dir"] / "redock", seed=42, nrun=20,
                                 cpu=8, pdb_id="1T69")
    assert res.ok, res.error
    assert res.ions == "ZN378A" and res.zn_mode and res.ion_source == "benchmark_file" and res.metal_audited
    # Regression guard. AutoDock4Zn leaves SAHA head-to-tail reversed in this pentacoordinate site
    # (about 7.8 A in place, consistent across seeds and sampling depth). A superposition-based metric
    # reports this very pose as 1.7 A ("validated"); in-place RMSD must keep rejecting it.
    assert res.ad_lowest_rmsd > 5.0
    assert not redock.is_success(res)


def test_in_place_rmsd_agrees_with_open_babel_and_is_not_the_superposed_value(hdac8):
    d = hdac8["dir"] / "obrms"
    rec = engines.prepare_receptor(hdac8["pdb"], hdac8["dir"] / "rec_zn", hdac8["box"], keep_ions=hdac8["ions"])
    lig = ligprep.prepare_ligand(SAHA, "saha", 42, chelate_zn=True)
    engines.run_autodock_gpu(rec, [lig], d, nrun=10, seed=42)
    dlg = (d / "ligands" / "saha.dlg").read_text()
    s = engines.summarize_dlg(dlg)
    pose = pose_to_mol(dlg, engines.dlg_run_blocks(dlg)[s["run"]])
    xtal = hdac8["heavy"]
    Chem.MolToMolFile(xtal, str(d / "xtal.sdf"))
    Chem.MolToMolFile(pose, str(d / "pose.sdf"))

    def obrms(*flags):
        out = subprocess.run([str(find_tool("obrms")), *flags, str(d / "xtal.sdf"), str(d / "pose.sdf")],
                             capture_output=True, text=True).stdout
        import re
        return float(re.search(r"RMSD\s*(?:=:=|:)\s*([0-9.]+)", out).group(1))   # format differs by input type

    mine, in_place, superposed = symmetry_rmsd(xtal, pose), obrms(), obrms("-m")
    assert mine == pytest.approx(in_place, abs=0.05)         # independent implementation agrees
    assert superposed <= in_place + 1e-6                     # aligning can only lower the RMSD


def test_fen1_loses_its_catalytic_magnesium_in_the_benchmark_file_and_the_plan_restores_it(tmp_path):
    root = ROOT / "data" / "raw" / "full_data"
    if not (root / "FEN1").is_dir():
        pytest.skip("LIT-PCBA data not downloaded")
    t = data.load_target(root, "FEN1")
    tm = t.templates[0]
    assert K.read_structure_ions(t.protein_mol2(tm)) == []                     # the benchmark file has no ions
    xtal, _ = redock.crystal_ligand(t.ligand_mol2(tm), tmp_path)
    heavy = Chem.RemoveHs(xtal).GetConformer().GetPositions()
    box = box_from_coords(heavy)
    plan = redock.plan_for(t.protein_mol2(tm), heavy, box, tm, None)
    assert plan.source == "restored" and plan.audited and not plan.zn
    assert sorted(i.element for i in plan.ions) == ["MG", "MG"]
    assert all(K.min_distance(i, heavy) < 3.0 for i in plan.ions)             # the ions bind the ligand
    rec = engines.prepare_receptor(t.protein_mol2(tm), tmp_path / "rec", box, keep_ions=plan.ions)
    assert not rec.zn_mode and rec.report.ions_kept == tuple(i.key() for i in plan.ions)
    assert (tmp_path / "rec" / "receptor.pdbqt").read_text().count(" Mg") == 2
    blind = redock.plan_for(t.protein_mol2(tm), heavy, box, tm, {"mode": "off"})
    assert blind.ions == () and "disabled" in blind.note
