import numpy as np
import pandas as pd
import pytest

from vsbench import admet as A

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
RHODANINE = "O=C1C(=Cc2ccccc2)SC(=S)N1"                 # ene-rhodanine, a textbook PAINS chemotype
# 12 rotatable bonds, low polarity
GREASY_CHAIN = "CCCCCCCCCCCCCCCC(=O)OC"
# TPSA far above 140 A^2 from many hydroxyls and amides
POLAR = "OCC(O)C(O)C(O)C(O)C(=O)NC(CO)C(=O)NC(CO)C(=O)NC(CO)C(=O)N"


def test_pains_flags_the_known_chemotype_and_not_aspirin():
    assert A.pains_alerts(ASPIRIN) == []
    hits = A.pains_alerts(RHODANINE)
    assert hits and any("ene_rhod" in h for h in hits)
    assert A.pains_alerts("not a smiles") == []


def test_veber_uses_rotatable_bonds_and_tpsa_not_lipinski():
    a = A.veber(ASPIRIN)
    assert a["veber_pass"] and a["tpsa"] == pytest.approx(63.6, abs=0.1) and a["rotb_veber"] <= 10
    g = A.veber(GREASY_CHAIN)
    assert g["rotb_veber"] > 10 and not g["veber_pass"]
    p = A.veber(POLAR)
    assert p["tpsa"] > 140 and not p["veber_pass"]
    bad = A.veber("not a smiles")
    assert np.isnan(bad["tpsa"]) and bad["veber_pass"] is False


@pytest.mark.parametrize("row,tier", [
    ({"n_pains": 0, "AMES": 0.1, "hERG": 0.1, "veber_pass": True}, "PASS"),
    ({"n_pains": 0, "AMES": 0.1, "hERG": 0.1, "veber_pass": False}, "MODERATE"),
    ({"n_pains": 0, "AMES": 0.1, "hERG": 0.9, "veber_pass": True}, "HIGH_RISK"),
    ({"n_pains": 1, "AMES": 0.1, "hERG": 0.1, "veber_pass": True}, "EXCLUDE"),
    ({"n_pains": 0, "AMES": 0.9, "hERG": 0.1, "veber_pass": True}, "EXCLUDE"),
    ({"n_pains": 0, "AMES": 0.9, "hERG": 0.9, "veber_pass": False}, "EXCLUDE"),        # worst tier wins
    ({"n_pains": 0, "AMES": 0.1, "hERG": 0.9, "veber_pass": False}, "HIGH_RISK"),      # hERG beats Veber
    ({"n_pains": 0, "AMES": 0.5, "hERG": 0.5, "veber_pass": True}, "PASS"),            # thresholds are strict ">"
])
def test_tier_rules(row, tier):
    assert A.assign_tier(row) == tier


def test_tier_risk_ordering():
    assert A.TIER_RISK["PASS"] < A.TIER_RISK["MODERATE"] < A.TIER_RISK["HIGH_RISK"] < A.TIER_RISK["EXCLUDE"]


def test_profile_with_injected_predictions_makes_no_model_call():
    pred = pd.DataFrame({e: [0.1, 0.9] for e in A.ENDPOINTS}, index=[ASPIRIN, RHODANINE])
    pred.loc[ASPIRIN, "hERG"] = 0.2
    out = A.profile([ASPIRIN, RHODANINE, "not a smiles"], predictions=pred)
    assert list(out.smiles) == [ASPIRIN, RHODANINE, "not a smiles"]
    assert out.loc[0, "tier"] == "PASS" and out.loc[1, "tier"] == "EXCLUDE"
    assert out.loc[2, "tier"] in {"PASS", "MODERATE"} or np.isnan(out.loc[2, "AMES"])
    assert out.loc[0, "risk"] == 0 and out.loc[1, "risk"] == 3
    assert out.loc[2, "veber_pass"] == False                                             # noqa: E712
    assert np.isnan(out.loc[2, "AMES"])                                                  # unknown SMILES -> NaN, not 0


@pytest.mark.slow
def test_admet_ai_predicts_probabilities_and_known_toxicity_ordering():
    pred = A.predict([ASPIRIN, "c1ccc2c(c1)ccc1ccccc12", "CCN(CC)CCCC(C)Nc1ccnc2cc(Cl)ccc12"])
    for e in A.ENDPOINTS:
        assert e in pred.columns, e
    sub = pred[["hERG", "AMES", "DILI", "HIA_Hou"]]
    assert ((sub >= 0) & (sub <= 1)).all().all()
    # chloroquine-like amine: far higher hERG probability than aspirin
    assert pred.loc["CCN(CC)CCCC(C)Nc1ccnc2cc(Cl)ccc12", "hERG"] > pred.loc[ASPIRIN, "hERG"]
    again = A.predict([ASPIRIN])
    assert again.loc[ASPIRIN, "hERG"] == pytest.approx(pred.loc[ASPIRIN, "hERG"], abs=1e-6)   # deterministic


def test_qed_is_continuous_and_orders_a_drug_above_a_greasy_chain():
    """Veber passes or fails; QED separates compounds that fall on the same side of it."""
    assert A.veber(ASPIRIN)["veber_pass"] and A.veber(POLAR)["veber_pass"] is False
    q_asp, q_greasy = A.qed(ASPIRIN), A.qed(GREASY_CHAIN)
    assert 0.0 <= q_asp <= 1.0 and 0.0 <= q_greasy <= 1.0
    assert q_asp > q_greasy
    assert np.isnan(A.qed("not a smiles"))


def test_qed_is_reported_but_does_not_change_the_pre_registered_tier():
    pred = pd.DataFrame({e: [0.1] for e in A.ENDPOINTS}, index=[GREASY_CHAIN])
    out = A.profile([GREASY_CHAIN], predictions=pred)
    assert "qed" in out.columns and out.loc[0, "qed"] == pytest.approx(A.qed(GREASY_CHAIN))
    assert out.loc[0, "tier"] == "MODERATE"            # Veber violation, unchanged by a low QED


def test_the_pharmacokinetic_endpoints_are_kept_not_discarded():
    """A liability shows up in one endpoint; keeping only the safety block hides it."""
    for e in ("Caco2_Wang", "PPBR_AZ", "Clearance_Hepatocyte_AZ", "Pgp_Broccatelli"):
        assert e in A.ENDPOINTS
