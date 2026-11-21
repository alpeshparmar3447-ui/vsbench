"""Ligand preparation: SMILES -> 3-D conformer -> docking-ready PDBQT.

Ionisation states are used exactly as the input SMILES supply them, and recorded;
nothing is re-protonated here, because a benchmark whose inputs were prepared at
physiological pH would otherwise be silently altered. A
failed embedding or PDBQT conversion is returned as a failure record, never
dropped, so the caller can report failure rates by class.
"""
from __future__ import annotations

import zlib
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

from . import metals as K

RDLogger.DisableLog("rdApp.*")


@dataclass(frozen=True)
class Ligand:
    name: str
    smiles: str
    ok: bool
    pdbqt: str = ""
    molblock: str = ""      # 3-D structure with explicit hydrogens (for gnina / RDKit)
    n_heavy: int = 0
    n_torsions: int = 0
    error: str = ""
    chelators_auto: str = ""      # recognised chelator classes (deprotonated for a Zn pocket)
    chelators_manual: str = ""    # chelating groups that need a hand-checked ionisation state
    deprotonated: bool = False


def ligand_seed(name: str, base_seed: int) -> int:
    """Deterministic per-ligand seed that does not depend on processing order."""
    return (zlib.crc32(name.encode()) ^ (base_seed * 2654435761)) & 0x7FFFFFFF


EMBED_TIMEOUT_S = 5


def embed(smiles: str, seed: int, max_attempts: int = 3, timeout: int = EMBED_TIMEOUT_S) -> Chem.Mol:
    """ETKDG embedding with a random-coordinate fallback, then MMFF relaxation.

    Each attempt is bounded by ``timeout`` seconds: RDKit's embedding has no default
    limit, and a single pathological molecule (observed once in 3,200 compounds)
    can otherwise stall a whole batch for minutes.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("RDKit could not parse the SMILES")
    if mol.GetNumHeavyAtoms() < 2:
        raise ValueError("fewer than two heavy atoms")
    mol = Chem.AddHs(mol)
    for attempt in range(max_attempts):
        params = AllChem.ETKDGv3()
        params.randomSeed = seed + attempt
        params.useRandomCoords = attempt > 0
        params.timeout = int(timeout)
        if AllChem.EmbedMolecule(mol, params) == 0:
            try:
                AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            except Exception:  # MMFF lacks parameters for some chemotypes; keep ETKDG geometry
                pass
            return mol
    raise ValueError("3-D embedding failed")


def to_pdbqt(mol: Chem.Mol) -> tuple[str, int]:
    """Meeko PDBQT string and the number of active torsions."""
    setups = MoleculePreparation().prepare(mol)
    if not setups:
        raise ValueError("meeko produced no molecule setup")
    pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise ValueError(f"PDBQT conversion failed: {err}")
    n_tors = sum(1 for line in pdbqt.splitlines() if line.startswith("BRANCH"))
    return pdbqt, n_tors


def prepare_ligand(smiles: str, name: str, base_seed: int = 42, timeout: int = EMBED_TIMEOUT_S,
                   chelate_zn: bool = False) -> Ligand:
    """SMILES -> docking-ready ligand.

    ``chelate_zn`` is set for a pocket that holds a zinc ion: recognised chelator groups
    (hydroxamic acid, catechol, hydroxypyridinone/tropolone) are then deprotonated, as
    AutoDock4Zn requires. Groups that need a hand-checked state (thiol, phosphonate,
    8-hydroxyquinoline, primary sulfonamide) are flagged and left unchanged. Chelator
    classes are recorded for every ligand, whatever the pocket.
    """
    auto = manual = ""
    try:
        parsed = Chem.MolFromSmiles(smiles)
        use, flipped = smiles, False
        if parsed is not None:
            flags = K.chelator_classes(parsed)
            auto, manual = ";".join(flags["auto"]), ";".join(flags["manual"])
            if chelate_zn and flags["auto"]:
                deprot, _ = K.deprotonate_chelators(parsed)
                use, flipped = Chem.MolToSmiles(deprot), True
        mol = embed(use, ligand_seed(name, base_seed), timeout=timeout)
        pdbqt, n_tors = to_pdbqt(mol)
        return Ligand(name=name, smiles=smiles, ok=True, pdbqt=pdbqt,
                      molblock=Chem.MolToMolBlock(mol), n_heavy=mol.GetNumHeavyAtoms(),
                      n_torsions=n_tors, chelators_auto=auto, chelators_manual=manual,
                      deprotonated=flipped)
    except Exception as exc:  # recorded, not raised: failures are data
        return Ligand(name=name, smiles=smiles, ok=False, error=str(exc)[:200],
                      chelators_auto=auto, chelators_manual=manual)


def _prepare_star(args):
    return prepare_ligand(*args)


def prepare_many(items, base_seed: int = 42, n_jobs: int = 1, chelate_zn: bool = False) -> list[Ligand]:
    """Prepare ``(smiles, name)`` pairs; output order matches input order."""
    jobs = [(smi, name, base_seed, EMBED_TIMEOUT_S, chelate_zn) for smi, name in items]
    if n_jobs > 1 and len(jobs) > 8:
        with ProcessPoolExecutor(n_jobs) as ex:
            return list(ex.map(_prepare_star, jobs, chunksize=16))
    return [_prepare_star(j) for j in jobs]
