"""Chemical-similarity and analogue-bias audit (plan section 4.10.4).

Docking has no training set, so classical train/test leakage cannot occur. What can occur is
analogue bias: actives that are close analogues of one another, and negatives that are or are
not close to the actives. These functions measure that, so a benchmark's difficulty can be
read next to its docking results. Fingerprints are Morgan radius 2, 2048 bits (ECFP4-like).
"""
from __future__ import annotations

import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.ML.Cluster import Butina

from . import metrics as M

RDLogger.DisableLog("rdApp.*")


def morgan_fps(smiles, radius: int = 2, nbits: int = 2048) -> list:
    out = []
    for s in smiles:
        m = Chem.MolFromSmiles(s)
        out.append(AllChem.GetMorganFingerprintAsBitVect(m, radius, nBits=nbits) if m is not None else None)
    return out


def max_similarity(query_fps, reference_fps, exclude_self: bool = False) -> np.ndarray:
    """Maximum Tanimoto similarity of each query to the reference set (NaN for unparsable queries).
    With ``exclude_self`` the query set is the reference set and a compound is not compared with itself."""
    refs = [f for f in reference_fps if f is not None]
    out = np.full(len(query_fps), np.nan)
    for i, q in enumerate(query_fps):
        if q is None:
            continue
        sims = np.array(DataStructs.BulkTanimotoSimilarity(q, refs))
        if exclude_self:
            sims = np.delete(sims, i) if len(sims) == len(query_fps) else sims
        out[i] = sims.max() if len(sims) else np.nan
    return out


def butina_cluster_sizes(fps, similarity_cutoff: float = 0.4) -> list[int]:
    """Cluster sizes (largest first) at a Tanimoto similarity cutoff, by Butina clustering."""
    fps = [f for f in fps if f is not None]
    n = len(fps)
    if n == 0:
        return []
    dists = []
    for i in range(1, n):
        dists.extend(1.0 - np.array(DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])))
    clusters = Butina.ClusterData(dists, n, 1.0 - similarity_cutoff, isDistData=True)
    return sorted((len(c) for c in clusters), reverse=True)


def murcko(smiles) -> list[str]:
    out = []
    for s in smiles:
        m = Chem.MolFromSmiles(s)
        try:
            out.append(MurckoScaffold.MurckoScaffoldSmiles(mol=m) if m is not None else "")
        except Exception:
            out.append("")
    return out


def scaffold_sharing(actives, others) -> float:
    """Fraction of ``others`` whose non-empty Bemis-Murcko scaffold also occurs among the actives."""
    act = {s for s in murcko(actives) if s}
    sc = [s for s in murcko(others) if s]
    return float(np.mean([s in act for s in sc])) if sc else float("nan")


def analogue_baseline_auc(active_smiles, negative_smiles, n_boot: int = 500, seed: int = 0):
    """ROC-AUC of 'similarity to the actives' used as a score (actives scored leave-one-out).
    Well above 0.5 means the negatives can be told from the actives by analogue structure alone."""
    act_fps, neg_fps = morgan_fps(active_smiles), morgan_fps(negative_smiles)
    s_act = max_similarity(act_fps, act_fps, exclude_self=True)
    s_neg = max_similarity(neg_fps, act_fps)
    y = np.r_[np.ones(len(s_act)), np.zeros(len(s_neg))].astype(int)
    s = np.r_[s_act, s_neg]
    ok = np.isfinite(s)
    return M.bootstrap_ci(M.roc_auc, y[ok], s[ok], n_boot=n_boot, seed=seed)


def reference_similarity(query_smiles, reference_smiles, query_is_reference=None) -> np.ndarray:
    """Max ECFP4 Tanimoto of each query to a reference set of known actives.

    Two reference sets give two very different baselines, and the distinction decides what the
    comparison means:

    *Retrospective (oracle).* The reference is the evaluation set's own actives, so a compound is
    scored against the very labels the benchmark is testing. Actives must then be scored
    leave-one-out, which ``query_is_reference`` does. This is an upper bound on what ligand
    similarity could extract from this data, not a method anyone could run prospectively.

    *Deployable.* The reference is a disjoint set of known actives, so no evaluation label is
    used and the score is one a prospective screen could actually compute. Pass the held-out
    actives as ``reference_smiles`` and leave ``query_is_reference`` None.
    """
    ref_fps = morgan_fps(reference_smiles)
    q_fps = morgan_fps(query_smiles)
    if query_is_reference is None:
        return max_similarity(q_fps, ref_fps)

    out = max_similarity(q_fps, ref_fps)
    mask = np.asarray(query_is_reference, dtype=bool)
    if not mask.any():
        return out
    # Leave-one-out drops the self-match by position, so the masked queries must be exactly the
    # reference set, in order. A mismatched call would otherwise score every active against a set
    # containing it and return 1.0, which reads as a strong baseline rather than as an error.
    if int(mask.sum()) != len(ref_fps):
        raise ValueError(
            f"leave-one-out needs the flagged queries to be the reference set: "
            f"{int(mask.sum())} flagged against {len(ref_fps)} references")
    out[mask] = max_similarity([f for f, m in zip(q_fps, mask) if m], ref_fps, exclude_self=True)
    return out


def audit(active_smiles, negative_sets: dict) -> dict:
    """One-call audit: actives' internal diversity, and each negative set's similarity to the actives."""
    act_fps = morgan_fps(active_smiles)
    internal = max_similarity(act_fps, act_fps, exclude_self=True)
    sizes = butina_cluster_sizes(act_fps, 0.4)
    out = {"actives": {
        "n": len(active_smiles), "median_nn_similarity": float(np.nanmedian(internal)),
        "frac_nn_ge_0.55": float(np.nanmean(internal >= 0.55)),
        "n_clusters_at_0.4": len(sizes), "largest_cluster": sizes[0] if sizes else 0,
        "scaffold_uniqueness": len({s for s in murcko(active_smiles) if s}) / max(len(active_smiles), 1)}}
    for name, smi in negative_sets.items():
        sims = max_similarity(morgan_fps(smi), act_fps)
        point, lo, hi = analogue_baseline_auc(active_smiles, smi)
        out[name] = {"n": len(smi), "median_max_similarity": float(np.nanmedian(sims)),
                     "frac_ge_0.4": float(np.nanmean(sims >= 0.4)), "frac_ge_0.55": float(np.nanmean(sims >= 0.55)),
                     "scaffold_sharing": scaffold_sharing(active_smiles, smi),
                     "analogue_baseline_auc": point, "analogue_baseline_ci": (lo, hi)}
    return out
