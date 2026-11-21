import numpy as np
import pytest

from vsbench import gates as G
from vsbench.ligprep import prepare_ligand


def _lig(smiles="CC(C)Cc1ccc(cc1)C(C)C(=O)O"):
    l = prepare_ligand(smiles, "x", 42)
    assert l.ok
    return l


def test_torsion_gate_flags_a_flexible_ligand_with_no_torsions_only():
    flex = _lig()
    assert G.torsion_gate(flex.smiles, flex.n_torsions)["frozen_flexible"] is False
    assert G.torsion_gate(flex.smiles, 0)["frozen_flexible"] is True                 # same ligand, torsions lost
    benz = _lig("c1ccccc1")                                                            # rigid: zero torsions is correct
    assert G.torsion_gate(benz.smiles, benz.n_torsions)["frozen_flexible"] is False
    amide = _lig("CC(=O)Nc1ccccc1")                                                    # an amide bond is legitimately frozen
    assert G.torsion_gate(amide.smiles, amide.n_torsions)["frozen_flexible"] is False


def test_charge_gate_only_fires_on_large_charges():
    p = _lig().pdbqt
    g = G.charge_gate(p)
    assert g["max_abs_charge"] < 1.0 and g["charge_flag"] is False
    stressed = G.scale_charges(p, 4.0)
    assert G.charge_gate(stressed)["charge_flag"] is True and G.charge_gate(stressed)["max_abs_charge"] > 1.0


def test_scaling_multiplies_charges_and_keeps_the_file_valid():
    p = _lig().pdbqt
    q0, q4 = np.array(G.pdbqt_charges(p)), np.array(G.pdbqt_charges(G.scale_charges(p, 4.0)))
    np.testing.assert_allclose(q4, np.clip(q0 * 4, -9.999, 9.999), atol=2e-3)
    assert G.scale_charges(p, 1.0).count("ATOM") == p.count("ATOM")
    assert all(len(l) >= 77 for l in G.scale_charges(p, 4.0).splitlines() if l.startswith("ATOM"))   # columns intact


def test_capping_clamps_only_the_excess_and_never_raises_a_charge():
    p = G.scale_charges(_lig().pdbqt, 4.0)
    capped = G.cap_charges(p, 1.0)
    q, qc = np.array(G.pdbqt_charges(p)), np.array(G.pdbqt_charges(capped))
    assert np.abs(qc).max() <= 1.0 + 1e-9
    below = np.abs(q) <= 1.0
    np.testing.assert_allclose(qc[below], q[below], atol=1e-9)                       # small charges untouched
    assert np.all(np.abs(qc) <= np.abs(q) + 1e-9)


def test_freezing_removes_the_torsion_tree_but_keeps_every_atom_and_coordinate():
    lig = _lig()
    assert lig.n_torsions >= 3
    frozen = G.freeze_torsions(lig.pdbqt)
    assert G.torsion_count(frozen) == 0 and "BRANCH" not in frozen and "ROOT" in frozen
    atoms = lambda t: [l for l in t.splitlines() if l.startswith("ATOM")]
    assert [l[30:54] for l in atoms(frozen)] == [l[30:54] for l in atoms(lig.pdbqt)]
    assert G.torsion_gate(lig.smiles, G.torsion_count(frozen))["frozen_flexible"] is True     # the gate catches it


def test_apply_gates_reports_every_ligand_including_failures():
    good, bad = _lig(), prepare_ligand("not a smiles", "bad", 42)
    rows = G.apply_gates([good, bad])
    assert [r["name"] for r in rows] == ["x", "bad"]
    assert rows[0]["rotb_strict"] >= 3 and not rows[0]["frozen_flexible"] and rows[0]["max_abs_charge"] < 1.0
    assert rows[1]["max_abs_charge"] != rows[1]["max_abs_charge"] and rows[1]["frozen_flexible"] is False   # NaN, not flagged
