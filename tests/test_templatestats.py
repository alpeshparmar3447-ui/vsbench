from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem

from vsbench import templatestats as T

DATA = Path(__file__).parents[1] / "data" / "raw" / "full_data"


def _write_mol2(path: Path, smiles: str):
    """RDKit cannot write mol2, so round-trip through a PDB block, which MolFromMol2File cannot read.
    Instead write a minimal mol2 by hand from the embedded molecule."""
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(m, randomSeed=7)
    conf = m.GetConformer()
    lines = ["@<TRIPOS>MOLECULE", "x", f"{m.GetNumAtoms()} {m.GetNumBonds()} 0 0 0", "SMALL", "USER_CHARGES",
             "", "@<TRIPOS>ATOM"]
    for a in m.GetAtoms():
        p = conf.GetAtomPosition(a.GetIdx())
        sym = a.GetSymbol()
        lines.append(f"{a.GetIdx()+1:>7} {sym}{a.GetIdx()+1:<6} {p.x:>9.4f} {p.y:>9.4f} {p.z:>9.4f} "
                     f"{sym:<5} 1 LIG 0.0000")
    lines.append("@<TRIPOS>BOND")
    order = {Chem.BondType.SINGLE: "1", Chem.BondType.DOUBLE: "2", Chem.BondType.TRIPLE: "3",
             Chem.BondType.AROMATIC: "ar"}
    for i, b in enumerate(m.GetBonds(), 1):
        lines.append(f"{i:>6} {b.GetBeginAtomIdx()+1:>4} {b.GetEndAtomIdx()+1:>4} "
                     f"{order.get(b.GetBondType(), '1')}")
    path.write_text("\n".join(lines) + "\n")


def test_properties_of_a_small_rigid_ligand(tmp_path):
    p = tmp_path / "a_ligand.mol2"
    _write_mol2(p, "c1ccccc1")
    d = T.ligand_properties(p)
    assert d["sanitized"] and d["heavy_atoms"] == 6 and d["n_rings"] == 1
    assert d["max_ring_size"] == 6 and d["macrocycle"] is False


def test_a_macrocycle_is_flagged(tmp_path):
    p = tmp_path / "m_ligand.mol2"
    _write_mol2(p, "C1CCCCCCCCCCCCC1")                           # 14-membered ring
    d = T.ligand_properties(p)
    assert d["max_ring_size"] == 14 and d["macrocycle"] is True


def test_a_missing_file_yields_unknown_properties_not_an_error(tmp_path):
    d = T.ligand_properties(tmp_path / "nope_ligand.mol2")
    assert d["sanitized"] is False and np.isnan(d["heavy_atoms"]) and d["macrocycle"] is None


def test_an_unsanitizable_ligand_keeps_its_atom_count_only(tmp_path):
    p = tmp_path / "bad_ligand.mol2"
    # A neutral nitrogen with five bonds: the same valence error the real 4jt5 ligand raises.
    atoms = ["      1 N1      0.0000   0.0000   0.0000 N.3   1 LIG 0.0000"]
    atoms += [f"{i:>7} C{i}      {i:.4f}   0.0000   0.0000 C.3   1 LIG 0.0000" for i in range(2, 7)]
    bonds = [f"{i-1:>6} {1:>4} {i:>4} 1" for i in range(2, 7)]
    p.write_text("@<TRIPOS>MOLECULE\nx\n6 5 0 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n"
                 + "\n".join(atoms) + "\n@<TRIPOS>BOND\n" + "\n".join(bonds) + "\n")
    d = T.ligand_properties(p)
    assert d["sanitized"] is False and np.isnan(d["mol_weight"])
    assert d["heavy_atoms"] == 6                                 # the count survives, nothing else does


def test_spearman_finds_a_monotone_relation_and_rejects_a_flat_one():
    x = np.arange(20.0)
    rho, p, n = T.spearman(x, x ** 3, n_perm=500, seed=0)
    assert rho == pytest.approx(1.0) and p < 0.01 and n == 20
    rho2, p2, _ = T.spearman(x, np.ones(20), n_perm=100)
    assert np.isnan(rho2)


def test_spearman_p_value_is_never_zero():
    rho, p, _ = T.spearman(np.arange(30.0), np.arange(30.0), n_perm=100, seed=0)
    assert p > 0


def test_spearman_drops_pairs_with_a_missing_value():
    x = np.array([1.0, 2.0, 3.0, 4.0, np.nan])
    y = np.array([1.0, 2.0, 3.0, 4.0, 99.0])
    rho, _, n = T.spearman(x, y, n_perm=100)
    assert n == 4 and rho == pytest.approx(1.0)


def _table():
    return pd.DataFrame([
        {"target": "A", "template": "1", "ran": True, "worst_rmsd": 0.5, "validated": True,
         "sanitized": True, "macrocycle": False, "heavy_atoms": 30.0, "mol_weight": 400.0,
         "rotatable_bonds": 4.0, "n_rings": 3.0, "max_ring_size": 6.0},
        {"target": "A", "template": "2", "ran": True, "worst_rmsd": 1.2, "validated": True,
         "sanitized": True, "macrocycle": False, "heavy_atoms": 32.0, "mol_weight": 430.0,
         "rotatable_bonds": 5.0, "n_rings": 3.0, "max_ring_size": 6.0},
        {"target": "B", "template": "3", "ran": True, "worst_rmsd": 9.0, "validated": False,
         "sanitized": True, "macrocycle": True, "heavy_atoms": 65.0, "mol_weight": 914.0,
         "rotatable_bonds": 6.0, "n_rings": 4.0, "max_ring_size": 29.0},
        {"target": "B", "template": "4", "ran": False, "worst_rmsd": np.nan, "validated": False,
         "sanitized": False, "macrocycle": None, "heavy_atoms": 23.0, "mol_weight": np.nan,
         "rotatable_bonds": np.nan, "n_rings": np.nan, "max_ring_size": np.nan},
    ])


def test_macrocycle_split_counts_every_template_including_the_unknown_one():
    t = T.by_macrocycle(_table())
    s = t.set_index("group")
    assert s.loc["non-macrocyclic", "n_templates"] == 2 and s.loc["non-macrocyclic", "validated"] == 2
    assert s.loc["macrocyclic", "n_templates"] == 1 and s.loc["macrocyclic", "validated"] == 0
    assert s.loc["ligand not sanitizable", "n_templates"] == 1
    assert int(t["n_templates"].sum()) == 4                      # nothing is lost
    assert 0 <= s.loc["non-macrocyclic", "ci_low"] <= 1


def test_correlations_use_only_sanitized_templates_with_an_rmsd():
    c = T.correlations(_table(), n_perm=200, seed=0)
    assert set(c.property) == set(T.PROPERTIES)
    assert (c.n_templates == 3).all()                            # the unsanitizable row is excluded


def test_per_target_flags_a_target_with_no_dockable_structure():
    p = T.per_target(_table()).set_index("target")
    assert p.loc["A", "validated"] == 2 and p.loc["A", "has_non_macrocyclic_template"]
    assert p.loc["B", "validated"] == 0 and p.loc["B", "macrocyclic"] == 1
    assert not p.loc["B", "has_non_macrocyclic_template"]        # only a macrocycle and an unreadable ligand
    assert p.loc["B", "ran"] == 1


@pytest.mark.skipif(not (DATA / "MTORC1").is_dir(), reason="LIT-PCBA data not downloaded")
def test_the_real_mtorc1_templates_are_mostly_macrocyclic():
    """The finding this module exists to record, measured from the deposited ligands themselves."""
    props = [T.ligand_properties(f) for f in sorted((DATA / "MTORC1").glob("*_ligand.mol2"))]
    macro = [p for p in props if p["macrocycle"]]
    assert len(props) == 11 and len(macro) == 9
    assert all(p["max_ring_size"] >= 29 for p in macro)
    assert all(900 <= p["mol_weight"] <= 1000 for p in macro)
