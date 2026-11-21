"""Template validation by redocking the co-crystal ligand, and label-free template choice.

Every receptor/box combination is validated against its own experimental pose
before any screening score is trusted. The ligand is rebuilt from SMILES with a
random conformer (never started from the crystal geometry), docked with both
engines using the production settings, and compared with the crystal pose by
in-place, symmetry-aware heavy-atom RMSD.

The template used for screening is chosen from these redocking results alone;
no active or inactive compound is consulted, so the choice cannot leak
benchmark labels.
"""
from __future__ import annotations

import traceback
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
from rdkit import Chem

from . import engines, metals
from .box import Box, box_from_coords
from .ligprep import prepare_ligand
from .poses import pose_to_mol
from .rmsd import symmetry_rmsd
from .tools import run_tool

SUCCESS_RMSD = 2.0  # conventional redocking success threshold, in angstrom

DEFAULT_METALS = {
    "mode": "auto",                       # 'auto' audits and retains ions; 'off' disables metal handling
    "pocket_cutoff": 8.0,                 # A from the ligand: functional ions kept
    "monovalent_cutoff": 6.0,             # A from the ligand: Na/K kept
    "restore_from_pdb": True,             # re-add ions the benchmark file lost, from the deposition
    "min_frame_overlap": 0.95,
    "pdb_cache": str(Path(__file__).resolve().parents[2] / "data" / "pdb_cache"),
}


def plan_for(protein_file: Path, ligand_xyz, box, pdb_id: str | None, metals_cfg: dict | None):
    """Metal plan for one template under the (partially overridden) default metal settings."""
    c = {**DEFAULT_METALS, **(metals_cfg or {})}
    return metals.plan_metals(protein_file, ligand_xyz, box, pdb_id, Path(c["pdb_cache"]), mode=c["mode"],
                              pocket_cutoff=c["pocket_cutoff"], monovalent_cutoff=c["monovalent_cutoff"],
                              restore_from_pdb=c["restore_from_pdb"], min_frame_overlap=c["min_frame_overlap"])


@dataclass
class RedockResult:
    template: str
    ok: bool
    error: str = ""
    n_heavy: int = 0
    ad_energy: float = np.nan
    ad_lowest_rmsd: float = np.nan      # lowest-energy pose
    ad_largest_rmsd: float = np.nan     # lowest-energy pose of the most populated cluster
    gn_vina: float = np.nan
    gn_top_rmsd: float = np.nan         # gnina top-ranked pose
    gn_best_rmsd: float = np.nan        # best of the returned poses (oracle, for reference)
    box_edge: float = np.nan
    ions: str = ""                      # ions retained in the receptor, e.g. 'ZN378A'
    ion_source: str = "none"            # benchmark_file | restored | mixed | none
    zn_mode: bool = False               # AutoDock4Zn used for AutoDock-GPU
    metal_audited: bool = False         # original deposition checked
    metal_note: str = ""
    extra: dict = field(default_factory=dict)

    def as_dict(self) -> dict:
        d = asdict(self)
        d.pop("extra")
        return d


def crystal_ligand(ligand_mol2: Path, workdir: Path) -> tuple[Chem.Mol, str]:
    """Crystal ligand as an RDKit molecule (heavy atoms + hydrogens) and its SMILES."""
    workdir.mkdir(parents=True, exist_ok=True)
    sdf = workdir / "crystal_ligand.sdf"
    run_tool("obabel", [str(ligand_mol2), "-O", str(sdf)], check=True)
    mol = Chem.MolFromMolFile(str(sdf), removeHs=False, sanitize=True)
    if mol is None:
        raise ValueError("RDKit could not sanitize the crystal ligand")
    smiles = Chem.MolToSmiles(Chem.RemoveHs(mol))
    return mol, smiles


def redock_template(protein_mol2: Path, ligand_mol2: Path, template: str, workdir: Path,
                    seed: int = 42, nrun: int = 20, exhaustiveness: int = 8, num_modes: int = 3,
                    cpu: int = 16, box_kwargs: dict | None = None, metals_cfg: dict | None = None,
                    pdb_id: str | None = None, box_override: Box | None = None) -> RedockResult:
    """Redock the co-crystal ligand. The box is normally built from the ligand itself (the oracle choice
    the screening campaign uses); ``box_override`` substitutes another definition for the grid-strategy study.
    The RMSD is always measured in place against the deposited pose."""
    res = RedockResult(template=template, ok=False)
    try:
        xtal, smiles = crystal_ligand(ligand_mol2, workdir)
        heavy = Chem.RemoveHs(xtal)
        res.n_heavy = heavy.GetNumAtoms()
        box = box_override or box_from_coords(heavy.GetConformer().GetPositions(), **(box_kwargs or {}))
        res.box_edge = max(box.size)
        plan = plan_for(protein_mol2, heavy.GetConformer().GetPositions(), box, pdb_id or template, metals_cfg)
        res.ions, res.ion_source, res.zn_mode = plan.label, plan.source, plan.zn
        res.metal_audited, res.metal_note = plan.audited, plan.note
        rec = engines.prepare_receptor(protein_mol2, workdir / "receptor", box, keep_ions=plan.ions)
        lig = prepare_ligand(smiles, f"{template}_xtal", seed, chelate_zn=plan.zn)
        if not lig.ok:
            raise ValueError(f"ligand preparation failed: {lig.error}")

        ad, _ = engines.run_autodock_gpu(rec, [lig], workdir / "adgpu", nrun=nrun, seed=seed)
        if not bool(ad.loc[0, "ad_ok"]):
            raise RuntimeError("AutoDock-GPU produced no result")
        dlg = (workdir / "adgpu" / "ligands" / f"{lig.name}.dlg").read_text()
        blocks = engines.dlg_run_blocks(dlg)
        res.ad_energy = float(ad.loc[0, "ad_energy"])
        res.ad_lowest_rmsd = symmetry_rmsd(heavy, pose_to_mol(dlg, blocks[int(ad.loc[0, "ad_run"])]))
        res.ad_largest_rmsd = symmetry_rmsd(heavy, pose_to_mol(dlg, blocks[int(ad.loc[0, "ad_largest_run"])]))

        gn, _ = engines.run_gnina(rec, [lig], workdir / "gnina", exhaustiveness=exhaustiveness,
                                  num_modes=num_modes, seed=seed, cpu=cpu)
        if not bool(gn.loc[0, "gn_ok"]):
            raise RuntimeError("gnina produced no result")
        poses = engines.read_gnina_sdf(workdir / "gnina" / "gnina_00000.out.sdf")[lig.name]
        rmsds = [symmetry_rmsd(heavy, p) for p in poses]
        res.gn_vina = float(gn.loc[0, "gn_vina"])
        res.gn_top_rmsd, res.gn_best_rmsd = float(rmsds[0]), float(min(rmsds))
        res.ok = True
    except Exception as exc:  # a template that cannot be validated is reported, not hidden
        res.error = f"{type(exc).__name__}: {exc}"
        res.extra["traceback"] = traceback.format_exc()[-800:]
    return res


def worst_engine_rmsd(r: RedockResult) -> float:
    return max(r.ad_lowest_rmsd, r.gn_top_rmsd)


def is_success(r: RedockResult, threshold: float = SUCCESS_RMSD) -> bool:
    return r.ok and worst_engine_rmsd(r) < threshold


def select_template(results: list[RedockResult]) -> RedockResult:
    """Template whose worst-engine redocking RMSD is smallest (ties: name order).

    Only validated templates are eligible; if none redock below the success
    threshold the best available is still returned and the caller must report
    that the target has no validated receptor.
    """
    valid = [r for r in results if r.ok and np.isfinite(worst_engine_rmsd(r))]
    if not valid:
        raise ValueError("no template could be redocked")
    return min(valid, key=lambda r: (worst_engine_rmsd(r), r.template))
