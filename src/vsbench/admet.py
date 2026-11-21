"""ADMET profiling and risk tiers (plan section 4.12).

Predictions come from ADMET-AI (Swanson et al., Bioinformatics 2024), run locally. Drug-likeness is
reported two ways: the Veber criteria (rotatable bonds <= 10 and TPSA <= 140 A^2) rather than
Lipinski's rule of five, and QED (Bickerton et al., Nat Chem 2012) as a continuous companion.
Against measured human oral bioavailability in a 588-compound clinical set, Veber compliance
outranks zero-violation Lipinski (Spearman 0.234 against 0.172) and QED outranks both (0.264),
which is why both are reported and Lipinski is not.

Risk tiers, fixed in advance, in this order:
  EXCLUDE    any PAINS alert, or Ames mutagenicity probability above 0.5
  HIGH_RISK  hERG inhibition probability above 0.5
  MODERATE   Veber violation
  PASS       none of the above
The 0.5 thresholds are the conventional cut-off of the binary classifiers; they are rules for ordering
compounds, not a validated clinical model. QED is reported alongside the tier and does not enter it:
the tier rule was fixed before the production run and is not re-tuned here.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import QED, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

# Safety and disposition endpoints kept from ADMET-AI's full output. The pharmacokinetic block is
# retained because endpoint values, not a single composite, are what identify a specific liability:
# a compound can be rejected for efflux or clearance while its summary score looks ordinary.
ENDPOINTS = ["hERG", "AMES", "DILI", "HIA_Hou", "BBB_Martins", "CYP3A4_Veith", "CYP2D6_Veith",
             "Solubility_AqSolDB", "Caco2_Wang", "PPBR_AZ", "Clearance_Hepatocyte_AZ",
             "Clearance_Microsome_AZ", "Pgp_Broccatelli", "VDss_Lombardo", "Half_Life_Obach"]
TIER_RISK = {"PASS": 0, "MODERATE": 1, "HIGH_RISK": 2, "EXCLUDE": 3}
HERG_CUT = 0.5
AMES_CUT = 0.5
VEBER_RB, VEBER_TPSA = 10, 140.0

_model = None
_catalog = None


def _pains_catalog() -> FilterCatalog:
    global _catalog
    if _catalog is None:
        p = FilterCatalogParams()
        for c in (FilterCatalogParams.FilterCatalogs.PAINS_A, FilterCatalogParams.FilterCatalogs.PAINS_B,
                  FilterCatalogParams.FilterCatalogs.PAINS_C):
            p.AddCatalog(c)
        _catalog = FilterCatalog(p)
    return _catalog


def pains_alerts(smiles: str) -> list[str]:
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return []
    return sorted({e.GetDescription() for e in _pains_catalog().GetMatches(m)})


def veber(smiles: str) -> dict:
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return {"rotb_veber": np.nan, "tpsa": np.nan, "veber_pass": False}
    rb = rdMolDescriptors.CalcNumRotatableBonds(m)
    tpsa = rdMolDescriptors.CalcTPSA(m)
    return {"rotb_veber": rb, "tpsa": tpsa, "veber_pass": bool(rb <= VEBER_RB and tpsa <= VEBER_TPSA)}


def qed(smiles: str) -> float:
    """QED, the continuous drug-likeness estimate of Bickerton et al. NaN for an unreadable SMILES.

    Continuous where Veber is a pass/fail, so two compounds that both pass, or both fail, are still
    ordered. It is computed from the structure alone and needs no prediction model.
    """
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return float("nan")
    try:
        return float(QED.qed(m))
    except Exception:                       # QED rejects a few exotic valences rather than scoring them
        return float("nan")


def assign_tier(row) -> str:
    if row.get("n_pains", 0) > 0 or row.get("AMES", 0.0) > AMES_CUT:
        return "EXCLUDE"
    if row.get("hERG", 0.0) > HERG_CUT:
        return "HIGH_RISK"
    if not row.get("veber_pass", True):
        return "MODERATE"
    return "PASS"


def predict(smiles: list[str]) -> pd.DataFrame:
    """ADMET-AI predictions for unique SMILES, indexed by SMILES (all endpoints it returns)."""
    global _model
    import warnings
    warnings.filterwarnings("ignore")
    from admet_ai import ADMETModel
    if _model is None:
        _model = ADMETModel()
    uniq = list(dict.fromkeys(smiles))
    pred = _model.predict(smiles=uniq)
    pred.index = uniq
    return pred


def profile(smiles: list[str], predictions: pd.DataFrame | None = None) -> pd.DataFrame:
    """One row per input SMILES: selected endpoints, Veber, PAINS and the risk tier."""
    pred = predictions if predictions is not None else predict(smiles)
    rows = []
    for s in smiles:
        row = {"smiles": s}
        row.update({e: float(pred.loc[s, e]) if (s in pred.index and e in pred.columns) else np.nan for e in ENDPOINTS})
        row.update(veber(s))
        row["qed"] = qed(s)
        alerts = pains_alerts(s)
        row["n_pains"], row["pains"] = len(alerts), ";".join(alerts)
        row["tier"] = assign_tier(row)
        row["risk"] = TIER_RISK[row["tier"]]
        rows.append(row)
    return pd.DataFrame(rows)
