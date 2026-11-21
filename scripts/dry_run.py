"""Data-path dry run at production scale: sampling, matching, ligand preparation and the
docking-free property baselines. No docking is performed."""
import sys, time, json
from pathlib import Path
import pandas as pd
from rdkit import Chem
from vsbench import data, ligprep, metrics as M
from vsbench.pipeline import Config, build_compound_table

cfg = Config.from_json(Path(sys.argv[1]))
rows, base_rows, bal_rows = [], [], []
for t in cfg.targets:
    t0 = time.time()
    target = data.load_target(Path(cfg.data_root), t)
    table, sets = build_compound_table(cfg, target)
    ligs = ligprep.prepare_many(list(zip(table.smiles, table.name)), cfg.seed, cfg.n_jobs)
    table["prep_ok"] = [l.ok for l in ligs]
    table["n_heavy"] = [Chem.MolFromSmiles(s).GetNumHeavyAtoms() for s in table.smiles]
    fails = table[~table.prep_ok]
    rows.append(dict(target=t, templates=len(target.templates), actives=int((table.neg_set=='active').sum()),
                     random=int((table.neg_set=='random').sum()), matched=int((table.neg_set=='matched').sum()),
                     unparsable=sets.n_unparsable, prep_failed=len(fails),
                     failed_by_class=dict(fails.neg_set.value_counts()), seconds=round(time.time()-t0,1)))
    bal = sets.balance().abs(); bal_rows.append(dict(target=t, **{f"max|SMD| {k}": round(bal[k].max(),2) for k in ("random","matched")},
                                                      worst_prop_matched=bal["matched"].idxmax()))
    for neg in ("random", "matched"):
        sub = pd.concat([table[table.neg_set=="active"], table[table.neg_set==neg]])
        y = sub.label.to_numpy()
        for prop in ("n_heavy","mw","clogp","tpsa"):
            p, lo, hi = M.bootstrap_ci(M.roc_auc, y, sub[prop].to_numpy(float), n_boot=200, seed=1)
            base_rows.append(dict(target=t, negative_set=neg, property=prop, auc=round(p,2), lo=round(lo,2), hi=round(hi,2)))
print(pd.DataFrame(rows).to_string(index=False)); print()
print(pd.DataFrame(bal_rows).to_string(index=False)); print()
b = pd.DataFrame(base_rows)
for prop in ("n_heavy","mw"):
    w = b[b.property==prop].pivot(index="target", columns="negative_set", values="auc")
    print(f"property-only ROC-AUC ({prop}):"); print(w.to_string()); print()
