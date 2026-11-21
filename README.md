# vsbench

**How much of a docking campaign's enrichment survives when the inactives are matched to the
actives on simple molecular properties?**

Retrospective virtual-screening benchmarks reward whatever separates actives from inactives, and
for docking that can be a one-dimensional property such as molecular size rather than recognition
of the binding site. `vsbench` measures the size of that confound. For each target it draws three
disjoint negative sets from **experimentally confirmed inactives** — a random draw, a
property-matched draw, and DUD-E-protocol decoys taken from other targets — docks them with two
independent scoring functions, and reports retrieval next to the performance that single molecular
properties achieve on their own and next to a size-controlled (conditional) AUC.

It is a measuring instrument, not a hit-finder. Use it to audit a screening protocol before you
trust its enrichment, or as a harness for a benchmark of your own.

What it measures in one run:

| | |
|---|---|
| **Retrieval** | ROC-AUC, PR-AUC, EF/NEF at 1% and 5%, BEDROC(α=20), each with stratified bootstrap intervals |
| **Bias controls** | standardized mean differences per property, single-property baselines scored as if they were docking scores, size-controlled (conditional) AUC |
| **Template choice** | redocking RMSD of every candidate template's own co-crystal ligand, measured **in place**; the template is chosen without any benchmark label |
| **Receptor quality** | chains, chain breaks against numbering gaps, missing backbone atoms, engineered residues, altlocs, resolution, coarse Ramachandran flags, metal-ion audit against the original deposition |
| **Binding site** | fpocket and DoGSiteScorer pockets, ligand-site recovery, druggability scores with positive and negative controls |
| **Protocol choices** | preparation-gate ablation, docking-box strategy (ligand-centred, pocket-centred, whole-protein blind), seed variance, score fusion variants, CPU-against-GPU throughput |
| **Downstream** | interaction-fingerprint recovery against the co-crystal pose, ADMET profiling, integrated lead prioritization |

Every analysis decision is fixed in [`PROTOCOL.md`](PROTOCOL.md) **before** the run it governs,
amendments included, so the reported numbers are not the best of several analyses.

---

## Install

```bash
python -m pip install vsbench                 # core: metrics, matching, docking, sites
python -m pip install 'vsbench[all]'          # adds interaction fingerprints, ADMET, test deps
```

From a checkout:

```bash
git clone <repository-url> && cd vsbench
python -m pip install -e '.[all]'
```

Python 3.10 or newer. Optional extras: `interactions` (PLIP), `admet` (ADMET-AI, pulls in torch),
`dev` (pytest, scikit-learn).

### External executables

Docking itself is done by programs that are not pip-installable. Put them on `PATH`, or in one
directory and set `VSBENCH_BIN` to it:

| Executable | Needed for | Source |
|---|---|---|
| `autodock_gpu` | AutoDock4 scoring on the GPU | [AutoDock-GPU](https://github.com/ccsb-scripps/AutoDock-GPU) |
| `autogrid4` | AutoDock4 affinity maps | [AutoDock4](https://autodock.scripps.edu/) |
| `gnina` | Vina-type docking with CNN rescoring | [gnina](https://github.com/gnina/gnina) |
| `mk_prepare_receptor.py` | receptor PDBQT preparation | [Meeko](https://github.com/forlilab/Meeko) (pip) |
| `obabel` | format conversion, 3D embedding | [Open Babel](https://openbabel.org/) |
| `fpocket` | binding-site detection (optional) | [fpocket](https://github.com/Discngine/fpocket) |
| `vina` | the CPU baseline in the performance study (optional) | [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina) |
| `prank` | local machine-learned pocket ranking (optional) | [P2Rank](https://github.com/rdk/p2rank) |
| `java` | required by P2Rank and nothing else (optional) | any JRE 11+ |

An NVIDIA GPU is required for AutoDock-GPU and gnina. P2Rank ships as a tarball holding a
`prank` script beside its models rather than a binary on `PATH`, so unpack it under `tools/` or
point `VSBENCH_P2RANK` at it; its Java runtime is found through `$JAVA_HOME` or a conda
environment named `java17`, `java`, `openjdk` or `jdk`.

Conda users: the search path includes `$CONDA_PREFIX/bin` and `~/*/envs/dock/bin`, so activating
an environment that holds the binaries is enough.

### Check the machine first

```bash
vsbench doctor
```

It prints where every executable resolved, which Python packages are present, and what the GPU is.
It exits non-zero if anything required is missing, so it works as a CI or setup gate.

---

## Running your own target

A **target** is a directory. Nothing is specific to LIT-PCBA:

```
my_data/
  MY_TARGET/
    actives.smi            "<SMILES> <ID>", one compound per line
    inactives.smi          same format; experimentally confirmed inactives, not generated decoys
    1abc_protein.mol2      receptor, including the residue SUBSTRUCTURE records
    1abc_ligand.mol2       that entry's co-crystal ligand, used to validate the template
    2def_protein.mol2      further candidate templates, as many as you have
    2def_ligand.mol2
```

The file stem before `_protein` is treated as a PDB identifier and is used to look up the original
deposition for the structure and metal audits. Every `*_protein.mol2` needs its matching
`*_ligand.mol2`: template selection works by redocking the deposited ligand, so a receptor without
one cannot be validated.

**Validate before you spend GPU time.** This opens and parses every file rather than checking that
it exists:

```bash
vsbench check-data --data-root my_data
```

It reports compound counts, templates, unparsable structures, an inactive pool too small for the
matching ratio you asked for, and — the failure that is otherwise silent — a receptor whose residue
identity cannot be read, which docks and scores normally while pocket detection and interaction
analysis quietly return nothing.

Then write a configuration and run it:

```json
{
  "data_root": "../my_data",
  "results_dir": "../results/mine",
  "targets": ["MY_TARGET"],
  "n_actives": 40,
  "n_random": 300,
  "per_active_matched": 8,
  "decoy_per_active": 0,
  "seed": 42
}
```

```bash
vsbench run    --config configs/mine.json
vsbench report --config configs/mine.json
```

Relative `data_root` and `results_dir` resolve against the configuration file's own directory, so
a checkout runs the same from any working directory. Unknown or missing keys are named, not
silently ignored.

Two things to know about the negative sets. The property-matched set is drawn from a candidate pool
**disjoint** from the random set, so the two comparisons are independent. The DUD-E-protocol decoy
arm draws from the *other* targets' inactives, so it needs at least two target directories; set
`decoy_per_active` to 0 with a single target.

---

## Commands

```bash
vsbench doctor                                 # environment check, no config needed
vsbench check-data   --data-root DIR           # validate data before docking
vsbench run          --config cfg.json         # templates -> screen -> analyze
vsbench report       --config cfg.json         # tables.md and figures/
vsbench ablation     --config cfg.json         # preparation-gate ablation
vsbench grid         --config cfg.json         # ligand / pocket / blind box strategies
vsbench perf         --config cfg.json         # CPU against GPU throughput
vsbench rescore      --config cfg.json         # size-controlled AUC, score fusion variants
vsbench interactions --config cfg.json         # interaction recovery for every compound
vsbench qc           --config cfg.json         # structure quality of every candidate template
vsbench sites        --config cfg.json         # pocket detection and druggability scoring
vsbench metals       --config cfg.json         # metal-ion audit against the depositions
vsbench boxdiag      --config cfg.json --target T --template 1abc --edges 16 20 24 30
```

`run` has three stages, selectable with `--stage`: `templates` redocks every candidate and selects
one using redocking RMSD alone, `screen` samples, prepares and docks, `analyze` builds the metric
tables. **Re-running resumes**: completed dockings are not repeated and timings accumulate, so an
interrupted run costs minutes rather than the whole campaign. `--targets` restricts any command to
a subset.

Outputs land under `results_dir`: `metrics_all_targets.csv` and `metrics_macro.csv` (retrieval with
intervals), `size_control_all_targets.csv`, `score_variants_*.csv`, `redock.csv` and `selected.json`
per target, and `figures/` as PNG, PDF and SVG. Tables are written publication-ready three ways:
`tables.md` carries all of them, `tables_main.md` the main-text set and `tables_supplementary.md`
the additional file. Numbering follows position, so main tables run from 1 and supplementary from
S1 and an absent optional table closes the gap rather than leaving one.

Two comparisons are worth reading before any enrichment figure. `baselines.csv` scores seven
physical descriptors and the nearest-neighbour similarity to the actives as if each were a docking
score, and `baseline_comparison.csv` puts every docking score against the strongest of them by
paired bootstrap, reporting both a pointwise and a multiplicity-adjusted interval. A docking result
that does not clear its own property and similarity baselines has not been shown to use the
structure, whatever its absolute ROC-AUC.

### Assessing the receptors

`qc`, `sites` and `metals` describe the *templates* rather than the screen, so they cover every
candidate template of a target rather than the one redocking selected, and their tables are
written beside the run directory and shared by every tier run against the same data. Add
`--all-targets` to cover every target under the data root, including ones this config does not
screen.

```jsonc
"receptors": {
  "fpocket": true,      // local pocket detection; on by default
  "p2rank": "auto",     // local ML pocket ranking; "auto" enables it when P2Rank and Java are found
  "dogsite": false,     // DoGSiteScorer: public web service, public PDB codes only, minutes each
  "controls": true,     // the four fixed druggability control structures
  "restore_metals": true,
  "offline": false,     // true forbids every outbound request
  "n_jobs": 2           // keep low when a web service is involved
}
```

### Running with no network at all

Everything that docks, scores or analyses is local. Two optional steps can reach out, and both
have a local answer:

| Reaches out | For | Instead |
|---|---|---|
| `proteins.plus` | DoGSiteScorer pocket ranking | **P2Rank**, which does the same job locally |
| `files.rcsb.org` | deposited structures for the metal audit and the QC fields taken from the header | a populated `pdb_cache`; `qc` never downloads at all |

`vsbench doctor` prints both, and says whether the local ranker is available. Setting
`"offline": true` forbids them: DoGSiteScorer is forced off even if the config asked for it, and
a deposition that is not already cached raises an error naming the file to drop in rather than
fetching it. Cached depositions keep working, so once the cache is filled a run makes no
outbound request at all.

`sites` writes `sites.csv`, one row per template and detector giving whether the detector's top
pocket, top three, or any pocket recovered the site the co-crystal ligand occupies, with that
pocket's druggability score and volume. With `controls` on it also writes `site_controls.csv`;
`receptors.druggability_verdict()` reads it and applies the rule fixed in advance, that a
druggability score counts as informative only where the detector scores the negative control
below every positive. On the published run neither detector passed that test, which is why no
druggability claim is made from them.

Two caveats worth knowing. DoGSiteScorer resolves a public accession, so it is only attempted
for templates whose identifier is a PDB code; a template named anything else is simply skipped
by that detector rather than failing. And fpocket estimates pocket volume by Monte Carlo, so
`site_pocket_volume` moves by a few parts in a thousand between runs on identical input, while
pocket count, ranking, recovery and druggability score are stable.

---

## Reproducing the paper

```bash
vsbench check-data --config configs/q1.json    # 5 targets, ~1 min
vsbench run        --config configs/q1.json    # the production screen
vsbench ablation   --config configs/q1.json
vsbench grid       --config configs/q1.json
vsbench perf       --config configs/q1.json
vsbench rescore    --config configs/q1.json
vsbench report     --config configs/q1.json
```

The paper's panel is the five LIT-PCBA targets with at least 25 actives and at most 70 000
inactives: ESR1_ant, MAPK1, MTORC1, PPARG, TP53. Download `full_data.tgz` from
[LIT-PCBA](https://drugdesign.unistra.fr/LIT-PCBA/) and extract it to `data/raw/full_data/`.
Smaller configurations are provided for trying the pipeline out: `configs/pilot.json` (2 targets),
`configs/reduced.json`, `configs/standard.json`, and `configs/metal_arm.json` (the FEN1 Mg²⁺ arm).

Reference cost on one RTX 3050 6 GB laptop GPU with 20 CPU threads: about 13 h of compute for
`configs/q1.json`, roughly 9 to 10 s per compound with both engines and three AutoDock-GPU seeds.
Wall-clock was 23 h because the host slept for 10 h of it, which is worth knowing if you time a
run: compare elapsed process time against the clock rather than trusting the clock.
`requirements-lock.txt` pins the exact Python environment the published numbers came from.

### Tests

```bash
pytest -q -m "not slow"     # no docking executables needed, about 40 s
pytest -q                   # whole suite, runs the real engines on the GPU, about 20 min
```

The real-engine tests build receptor grids and write docking output into pytest's temporary
directory, which costs roughly 500 MB per run and is kept for the last few runs. On a system
whose `/tmp` is a tmpfs, as it is by default under WSL2, that comes out of RAM and a few
consecutive runs will fill it. Point the scratch at a disk-backed directory first:

```bash
TMPDIR=~/.cache/vsbench-pytest pytest -q
```

---

## Design decisions worth knowing

* Higher score always means "more likely active". Docking energies are negated once, at the edge.
* A failed preparation or docking gets the engine's worst score and is tabulated, never dropped —
  dropping failures inflates every metric.
* Template choice uses redocking RMSD only. No benchmark label can influence it.
* Redocking RMSD is measured **in place**. Superposing the pose onto the reference first, which
  `obrms -m` does by default, called a pose 7.9 Å from the crystal "1.7 Å".
* The default box is generous, 6 Å buffer and 24 Å minimum edge: a tight box penalised large
  ligands with clash energies above +60 kcal/mol in a pilot.
* Consensus is the mean of per-engine z-scores. Treat a narrow gap between ranks accordingly.
* Metal ions are audited against the original deposition. Ions the benchmark file stripped are
  restored when the coordinate frames agree, retained ions must survive preparation or the stage
  fails loudly, and zinc pockets use AutoDock4Zn.
* Index groups, pocket identities and receptor groups are referenced by name, never by position.

## Known limitations

* **Statistical power.** With a few dozen actives per target, a bootstrap AUC interval is roughly
  ±0.10 wide. A result consistent with chance is also consistent with modest enrichment; read the
  interval, not the point estimate.
* **Bootstrap replicates resample compounds, not scaffolds.** Actives for one target are often
  analogues, so the effective sample size is below the nominal one and intervals are if anything
  too narrow. That bias favours finding a difference, so it weakens a positive result rather than
  a null one.
* Comparing many score variants against one reference inflates the chance of a spurious winner.
  The variant comparison reports a multiplicity-adjusted interval beside the pointwise one and
  promotes only on the adjusted one, which needs deeper bootstrap tails — see `--variant-n-boot`.
* Property matching covers seven one-dimensional descriptors. It does not match on shape,
  scaffold or synthetic accessibility, so a confound in those dimensions would survive it.
* One receptor template per target by default. Ensemble docking is not implemented.
* The DUD-E-protocol decoys are confirmed inactive against *other* targets, not against this one.
  That is the DUD-E assumption, stated rather than removed.
* Interaction fingerprints depend on PLIP recognising the residue identity in the receptor file.
  Where it cannot, the pass returns nothing instead of failing; `vsbench check-data` is the guard.
* The druggability scores were not calibrated: four controls cannot calibrate a score, and in this
  panel they did not separate the positive from the negative control.
* DoGSiteScorer needs network access and a public PDB identifier, so it cannot score a private or
  modelled structure.

## Citing

See [`CITATION.cff`](CITATION.cff). Please cite the accompanying paper as well as the software.

## License

MIT — see [`LICENSE`](LICENSE).
