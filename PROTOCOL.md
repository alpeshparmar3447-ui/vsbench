# Pre-specified protocol

Written before any benchmark compound was docked. Nothing in this file may be
changed after the production run starts without recording the change and the
reason in the manuscript. (The pilot dockings used only to measure timing and to
test the code, 20 ESR1_ant compounds, were not scored for retrieval.)

## Question
How much of the retrieval performance that two independent docking scoring
functions achieve on experimentally confirmed inactives survives when the
inactives are matched to the actives on simple one-dimensional properties, and
how does it compare with what those properties achieve on their own?

## Data
LIT-PCBA, full data set (Tran-Nguyen, Jacquemard and Rognan, J Chem Inf Model
2020, 60, 4263). Actives and inactives are confirmed in PubChem dose-response
assays. Ligand ionisation states are used as supplied.

## Target panel (rule, not a list chosen by results)
Every LIT-PCBA target with at least 25 actives and at most 70,000 inactives:
ESR1_ant, MAPK1, MTORC1, PPARG, TP53. All five are reported whatever their
results. Targets with fewer than 25 actives (ADRB2, ESR1_ago, OPRK1) are
excluded because a ROC-AUC interval from fewer than 25 actives is uninformative.
Targets with more than 70,000 inactives are excluded only for compute cost.

## Receptor template (label-free)
Each template's co-crystal ligand is redocked with both engines under the
production settings, starting from a random conformer built from SMILES. The
template with the smallest worst-engine RMSD (in-place, symmetry-aware,
heavy atoms) is used; ties break by template name. No active or inactive
compound is consulted. A target whose best template does not reach 2.0 A on both
engines is reported as having no validated receptor and its screening results
are flagged.

## Box
Cubic, centred on the heavy-atom centroid of the selected template's ligand;
edge = largest ligand span + 12 A, at least 24 A; grid spacing 0.375 A.

## Compounds per target
All actives if fewer than 40, otherwise 40 drawn at random (seed 42). Two
separately drawn sets of inactives, disjoint from each other:
* random: 300 inactives drawn uniformly from all inactives;
* matched: 8 inactives per active, nearest unused neighbours in standardised
  molecular weight, cLogP, HBD, HBA, rotatable bonds, TPSA and net charge,
  drawn from a disjoint candidate pool of up to 20,000 inactives.

## Engines and primary scores (fixed in advance)
* AutoDock-GPU v1.6 (AutoDock4 scoring): 20 LGA runs, heuristics on, autostop
  on, seed 42. Primary score: lowest-energy pose. Secondary: lowest energy of
  the most populated cluster.
* gnina v1.3.3: exhaustiveness 8, 3 modes, CNN rescoring, seed 42. Primary
  score: CNN_VS (CNN pose score x CNN affinity) of the top-ranked pose.
  Secondary: the Vina-type minimised affinity.
* Consensus: mean of the two primary z-scores.
A failed preparation or docking receives that engine's worst finite score and
is counted in the failure table.

## Metrics
ROC-AUC (primary), average precision, BEDROC (alpha = 20), EF1%, EF5%, NEF1%,
each with a 1,000-replicate stratified percentile bootstrap (95%). The primary
comparison is ROC-AUC on the matched set against ROC-AUC on the random set.

## Controls
* Single-property baselines (heavy atoms, MW, cLogP, TPSA, HBA, HBD, rotatable
  bonds) scored as if they were docking scores, on both negative sets.
* Paired bootstrap of the pose-selection effect (largest cluster minus lowest
  energy) and of the engine difference (gnina minus AutoDock-GPU).

## Interpretation rules
* A docking AUC is described as retrieval above chance only if its interval
  excludes 0.5 on the matched set.
* A docking AUC is described as exceeding a property baseline only if the
  paired comparison supports it; otherwise the baseline is reported as
  matching it.
* Macro-averages over the five targets are descriptive; no pooled significance
  test across targets is claimed.
* Results are reported for all five targets, and for every metric listed above.

## Amendment 1 (2026-11-10, before the production run): metal ions
Made at the request of the project lead after the harness was found to have no metal
logic. No benchmark compound had been docked. The panel, sample sizes, template rule
(apart from the receptor now carrying ions where the rules below require them), engines,
primary scores, metrics and interpretation rules are unchanged.

### Template audit and ion retention
* Each template's original PDB deposition is downloaded. Ions are considered only if at
  least 95% of the benchmark file's C-alpha atoms lie within 0.5 A of a deposited C-alpha
  (shared coordinate frame). Otherwise nothing is restored and the template is flagged.
* Retained: Zn, Mg, Mn, Fe, Ca, Cu, Ni, Co within 8 A of any heavy atom of the template
  ligand or inside the search box enlarged by 2 A; Na and K only within 6 A of the ligand.
  All other ions are reported in the audit and not retained.
* Retained ions are written into the receptor as explicit atoms. Preparation fails if a
  requested ion is absent from the receptor PDBQT. An unavailable deposition is recorded in
  the output, never silently ignored.

### Treatment
* Zinc pocket: AutoDock-GPU grids use AutoDock4Zn (tetrahedral pseudo-atom, zinc charge
  removed, six coordination terms). gnina uses the plain receptor with the ion.
  Hydroxamic acids, catechols and hydroxypyridinones/tropolones are deprotonated; thiols,
  phosphonate hydroxyls, 8-hydroxyquinolines and primary sulfonamides are flagged and left
  unchanged.
* Any other retained ion: generic non-bonded atom in both engines (no dedicated parameters
  exist). The result is flagged as such.
* A template that cannot be validated by redocking is reported as having no validated
  receptor, whatever its ion state. Metal-aware preparation is not a claim of metal-aware
  accuracy.

### Redocking metric
In-place, graph-matched, heavy-atom RMSD with no superposition. A superposing metric (for
example `obrms -m`) measures conformation, not placement: it reports an AutoDock-GPU pose
of SAHA in HDAC8 that is 7.9 A from the crystal pose, with the ligand bound head to tail,
as 1.7 A. The in-place implementation is checked against `obrms` without `-m` in the test
suite.

### Audit outcome for the panel (recorded before docking)
No ion is retained for ESR1_ant, MAPK1, MTORC1 or PPARG (none of their 56 depositions
contains an ion). The TP53 structural zinc lies 29.0-29.7 A from the ligand (where within
30 A), outside the pocket and the search box, and is not retained. The amendment therefore
does not change the receptor of any panel target.

### Additional arm, reported separately
FEN1 (template 5fv7: two catalytic Mg2+ at 2.29 and 2.31 A from the ligand, absent from
the benchmark file) with the same sampling as the panel (40 actives, 300 random and 8
matched inactives per active). It was excluded from the panel rule only by its inactive
count. Its purpose is to exercise the metal path on a real pocket ion. It is excluded from
the macro-average, and its template is reported as validated or not like any other.
Configuration: configs/metal_arm.json.

## Amendment 2 (2026-11-10, before the production run): extended pipeline

Everything below is fixed before any production docking is run. None of it changes the
template rule, the engine set, the primary scores or the interpretation rules above; the
panel result is still defined by the original design. The additions measure what that
result depends on and what a user can do with it.

### Negative-set designs (three, all docked, all reported)
1. random experimentally confirmed inactives (original);
2. property-matched experimentally confirmed inactives (original; drawn from the full inactive
   pool before docking, never matched after the fact);
3. DUD-E-protocol decoys (Mysinger et al. 2012): for each active, the 100 nearest unused
   compounds in standardised property space (mw, clogp, hbd, hba, rotb, tpsa, charge) are taken
   from a pool of inactives of the other LIT-PCBA targets, the 25% most similar (ECFP4
   Tanimoto) to any active are discarded, and the nearest remaining compound is used;
   8 per active. Pool compounds were never tested against this target, so a pool compound
   could in principle be active; decoys are labelled presumed inactive, not confirmed
   inactive. DeepCoy (generative decoys) is not used: it needs a trained model and an
   environment that this machine does not have, and its published gain over DUD-E-style
   decoys on Vina is small (AUC 0.70 to 0.63 in its own paper). This is a deviation from
   the plan, reported as such.

The primary negative set remains the matched set. The decoy arm is secondary. The
difference in ROC-AUC between designs is reported per target with an unpaired-negatives
bootstrap, and described as a design effect only if the interval excludes 0.

### Replication and sampling
* AutoDock-GPU is run with three seeds (42, 43, 44). The primary score is the mean of the
  per-seed lowest energies; per-seed ROC-AUC range, rank correlation between seeds and the
  fraction of the top 20 shared by all seeds are reported.
* gnina is run once for every compound (seed 42) and re-run with two further seeds on a
  subset of 60 compounds (all actives of the target up to 60) to measure its seed variance.
* Bootstrap: 1000 stratified resamples, 95% percentile intervals, seed fixed in the
  configuration.

### Structure quality (flag, never exclude)
Every template is characterised for chain breaks not explained by a numbering gap, missing
backbone atoms, number of chains, resolution, engineered residues, expression tags, altlocs,
hetero groups near the ligand and Ramachandran core/allowed/outside fractions (coarse regions;
a screening flag, not MolProbity). Flags are reported beside the result. They do not alter
template selection, which stays label-free and redocking-based.

### Site and druggability assessment
* fpocket on every template; DoGSiteScorer (public service, public PDB codes) on the five
  panel targets, FEN1 and the controls; P2Rank only if its execution is permitted on this machine.
  FTMap and CryptoSite are not automatable and are not run (a deviation, reported).
* A pocket recovers the ligand site if it contains at least half of the protein residues
  within 4.5 A of the co-crystal ligand. Reported: top-1 and top-3 recovery per detector with
  Wilson intervals, and the best-pocket recovery.
* Four controls, fixed in advance: positive HIV-1 protease/saquinavir (1HXB), EGFR
  kinase/erlotinib (1M17), carbonic anhydrase II/sulfonamide (1CIL); negative PTP1B/inhibitor
  (1BZC), a polar and historically difficult pocket. A detector's druggability score is
  described as informative only if it ranks the negative control below the positives.
  Detector scores are not used to predict docking performance beyond a descriptive
  correlation, following the caution of Michel et al. (2020) that such scores depend on the
  pocket definition.

### Quality gates and ablation
Two gates act on ligand preparation: a frozen-flexible-torsion gate and a partial-charge
gate (cap on the magnitude). Their effect is measured by docking 30 actives and 30 matched
inactives per panel target under four conditions: baseline, charges scaled by 4 without the
cap, charges scaled by 4 with the cap, all torsions frozen. Reported: how many compounds each
gate flags and the mean shift in docking energy against the baseline.

### Docking-box strategies
The co-crystal ligand of every template of the panel targets is redocked with (a) the ligand
box used by the campaign, (b) a 24 A cube at the centroid of the alpha spheres of the
top-ranked fpocket pocket, and (c) one box over the protein bounding box plus 2 A (grid
spacing coarsened only if 256 points per axis would be exceeded). RMSD is in place. Reported:
success rate (RMSD at or below 2 A) per strategy and engine with Wilson intervals. The
campaign's box is the oracle choice; (b) and (c) show what it is worth.

### Interaction analysis and lead selection
* PLIP interaction sets (type, chain, residue) of the co-crystal ligand are the reference;
  the top gnina pose of every active and of the 100 best-scoring compounds of each evaluation
  set is compared with it (recall, precision, Tanimoto, residue-level recall). A pose PLIP
  cannot analyse reproduces nothing and is counted, not dropped.
* ADMET-AI predictions (hERG, AMES, DILI, HIA, BBB, CYP3A4, CYP2D6, solubility), PAINS
  (A, B, C) and the Veber rules (rotatable bonds at most 10, TPSA at most 140 A^2; not
  Lipinski). Tiers, in order: EXCLUDE (any PAINS alert, or AMES probability above 0.5),
  HIGH_RISK (hERG probability above 0.5), MODERATE (Veber violation), PASS.
* Four strategies are compared on the pool of the 100 best consensus-scoring compounds of
  each evaluation set, chosen before any label is used: docking alone; docking after removing
  EXCLUDE and HIGH_RISK; docking plus interaction recall (mean z-score); Pareto ordering on
  docking, interaction recall and ADMET risk. Outcome: fraction of actives among the top 10,
  20 and 50, with bootstrap intervals. A filter that removes actives is penalised by the
  metric. ADMET predictions are ranking aids from a model trained on other compounds; they are
  not claimed to be correct for these compounds.

### Performance
AutoDock-GPU and gnina (GPU) at 8, 16 and 32 ligands; gnina with CPU-only CNN scoring, gnina
with Vina-type scoring on CPU, and AutoDock Vina (one process per ligand, parallel across
ligands, the most favourable CPU arrangement) at 1, 4, 8 and 16 threads on 16 ligands. Three
repeats after a warm-up run; wall time of the engine call only; GPU memory and utilisation
sampled every 0.5 s. Vina is a throughput reference only and does not enter any scoring
comparison. Results are specific to the recorded hardware (RTX 3050 6 GB laptop GPU, 20
logical CPUs) and are not claimed to transfer.

### Reporting
Each deviation above (decoy source, FTMap/CryptoSite, P2Rank if blocked, DoGSite on the
panel, FEN1 and controls only) is stated in the Methods and Limitations. Negative or null outcomes of any added
analysis are reported with the same prominence as positive ones.

## Amendment 3 (2026-11-11, after the screen, before any of the analyses below are computed)

The screening stage of the q1 campaign finished on 2026-11-11 and its outputs are frozen.
Nothing below re-runs it or changes a score already recorded. Four analyses are added, all of
them specified here before being computed, and one targeted diagnostic.

**Anti-tuning rule, governing every item in this amendment.** The primary score of the paper
remains the one registered before docking: the per-engine z-score mean of the AutoDock-GPU
lowest energy and the gnina CNN_VS score. Variants introduced below are reported as a
secondary comparison with paired bootstrap intervals. The best-performing variant is not
promoted to primary, and no variant, stratum, box or template is chosen by inspecting an
activity label or a resulting enrichment. Where a choice is made, it is made on label-free
evidence only (redocking geometry, molecular properties), and the rule making it is stated
here before the result is known.

### A. Template difficulty and what the redocking gate measures
For every template redocked in this project (the five panel targets, FEN1, and the
druggability controls), the in-place redocking RMSD of both engines is related to properties
of that template's co-crystal ligand: heavy-atom count, rotatable bonds, molecular weight,
largest ring size, and a macrocycle flag (largest ring of twelve atoms or more). Reported per
target: template count, the number that ran, the number validated at or below 2 A, and the
median of each property. Across templates: the Spearman correlation of the worst-engine RMSD
with each property, and the validation rate split by the macrocycle flag with Wilson
intervals. A ligand whose mol2 record cannot be sanitised is counted as a template whose
properties are unknown and is excluded from the correlation, not from the counts.
*Pre-specified interpretation:* a target whose templates are predominantly macrocyclic is one
that rigid-receptor docking cannot validate, and its screening result is reported as
unvalidated whatever its enrichment. This is a statement about the benchmark's structures,
not about the engines.

### B. Size-controlled discrimination
Docking scores rise with molecular size, so an enrichment measured against a negative set
that is smaller on average is partly a size measurement. For every target and every negative
set: the heavy-atom distribution of the pooled compounds is cut at its quartiles, and the ROC
AUC of each score is computed within each stratum and pooled as the mean weighted by the
number of active-inactive pairs the stratum contributes. Restricting the comparison to pairs
of equal size class and pooling by pair count are the same operation, so this is one
statistic, not two: the conditional AUC, the probability that an active outranks an inactive
of its own size class. It is reported with a stratified bootstrap interval beside the
unstratified value and the per-stratum values, and the heavy-atom AUC already registered as a
property baseline is reported in the same table. A stratum holding no active or no inactive
contributes no pairs and is reported as empty rather than counted as chance.
*Pre-specified interpretation:* where the size-controlled value is materially below the
unstratified one, the unstratified value is inflated by molecular size, the size-controlled
value is quoted as the discrimination attributable to the docking score, and both appear in
the abstract if the difference changes the paper's conclusion. A size-controlled AUC near 0.5
is reported as no discrimination, not as a negative result to be explained away.

### C. Scores already computed but not yet evaluated, and score fusion
The screen recorded two gnina scores that no table evaluates: the CNN score and the CNN
affinity head. Both are evaluated as standalone scores on every target and negative set,
with the same metrics and intervals as the registered scores. Fusion variants, all computed
from per-compound scores already on disk: the mean of the per-engine ranks; reciprocal-rank
fusion; a z-score mean weighted, per engine, by that engine's redocking quality on the
selected template (weight one over one plus its in-place RMSD); and each single engine alone.
Metrics: ROC AUC, EF1 and BEDROC20 per negative set, each compared with the registered
consensus by paired bootstrap on the same compounds.
*Pre-specified interpretation:* a fusion variant is reported as an improvement only where its
paired interval excludes zero on the property-matched and decoy sets, not on the random set
alone. The registered consensus stays primary regardless of the outcome.

### D. Interaction recovery as a score over the full evaluation sets
PLIP interaction recovery was computed for the actives and for the hundred best-scoring
compounds of each set, which is enough to compare selection strategies within a pool but not
to rank a whole set. Interaction recovery against the co-crystal reference is therefore
computed for every screened compound's top gnina pose, and evaluated both as a standalone
score and fused with the registered consensus, with the metrics and intervals of item C. A
pose PLIP cannot analyse recovers nothing and is ranked last, as elsewhere.
*Pre-specified interpretation:* as item C. Interaction recovery against a single co-crystal
reference rewards compounds that bind like that ligand, which is a bias to state, not a
ranking to trust blindly.

### E. Diagnostic: the one MTORC1 template with a drug-sized ligand
Nine of the eleven MTORC1 templates carry a 29-membered macrocycle of molecular weight 914 to
980. The exception, 4jsx, carries a 432 Da ligand with two rotatable bonds, and still failed
redocking at 6.13 A (AutoDock-GPU) and 6.21 A (gnina) with good scores, no pose among those
generated closer than 5.79 A, and nothing removed in receptor preparation. Its co-crystal
ligand is redocked with box edges of 16, 20, 24 and 30 A centred on the ligand, and with the
fpocket-centred box of the grid-strategy study, both engines, three seeds. Reported: RMSD
against box edge for both engines.
*Pre-specified interpretation:* redocking uses no activity label, so a box chosen because it
recovers a crystal pose is a label-free choice and is permitted. If a box recovers the pose at
or below 2 A, MTORC1 is re-screened with that box, both the original and the re-screened
result are reported, and the change of box is stated in the Methods. If none does, MTORC1
remains unvalidated and is reported as a target this protocol cannot validate. The diagnostic
runs only after the performance study has finished, so that no measured wall time includes it.

### C2. Ligand efficiency (added 2026-11-11, after item B, before being computed)
Item B showed that the AutoDock4 score's apparent retrieval against random inactives is largely
molecular size. The standard medicinal-chemistry response is to divide the score by ligand
size, so this is added as a score variant rather than left for a reader to ask about. Two
normalizations are evaluated for each engine's energy-like score: the energy per heavy atom
(ligand efficiency) and the energy divided by the cube root of the heavy-atom count, which is
the usual size-independent form. Both are computed from values already on disk, and both are
evaluated with the metrics, intervals and promotion rule of item C.
*Pre-specified interpretation:* dividing by size must reduce the retrieval measured against
random inactives, because that is the confound it removes, so a fall there is the expected
behaviour and not evidence against the normalization. The question it is being asked is
whether it helps on the property-matched and decoy sets, where the size gap is already near
zero; an improvement there would have to come from something other than removing a size
difference that is not present. The registered consensus stays primary either way.

## Correction 1 (2026-11-12): residue identity lost for two templates

Found while auditing the interaction-recovery output, not reported by the pipeline. It is
recorded here because it changed numbers that an earlier draft of the manuscript already
contained.

**Cause.** A protein mol2 names each residue in its `SUBSTRUCTURE` records. Those records
normally end `... PRO A 206`, with the chain and the residue number in separate fields. Where a
residue number reaches four digits the chain and the number are written joined, `... GLU B1385`,
one field shorter. Both readers of that section required the longer layout, so every record
written the short way was discarded and its residue became `UNK`, on a chain of its own,
renumbered by position in the file.

The layout is per record, not per file, so the damage was proportional to how many of a
structure's residues are numbered above 999. Of the 129 panel templates, 28 contain at least one
such record. In 12 the only ones affected are crystallographic waters. In 16 real residues are
affected: the FKBP12 complexes of MTORC1, whose mTOR fragment is numbered above 2000 while FKBP12
is not, five ADRB2 constructs, and 3sa0. In two, 4jsx and 4jt5, every residue is affected, because
those structures are the mTOR kinase domain alone and its residues run 1385 to 2549. 4jsx is the
template MTORC1 screened against, which is why the failure surfaced there first and completely.

**What it invalidated.** Interaction recovery for MTORC1: PLIP cannot type an unidentifiable
residue, so the co-crystal reference held zero contacts and all 980 compounds scored as missing
values, which reads as an uninformative target rather than a broken input. Site detection on the
affected templates: the residues lining the ligand were keyed to a numbering the detectors never
saw, so recovery was understated wherever a lining residue was numbered above 999, and was 0.000
for 4jsx and 4jt5, the only four such rows in the study. Structure quality on the same templates:
where numbering was synthetic and consecutive a numbering gap cannot exist, so chain breaks that
a gap does explain were reported as unexplained.

**What it did not touch.** Docking and screening read coordinates, not residue names, and the box
centre comes from the co-crystal ligand, so every docking score, redocking RMSD and enrichment
figure is unaffected. fpocket returned the same pocket count and the same top pocket before and
after the correction, so the grid-strategy study is unaffected.

**A second defect the same cause had hidden.** Crystallographic waters are numbered in the
thousands, so they are written the short way too, and their names were lost wherever that
happened. A water then became a residue of unknown name on a chain of its own, and because a
water carries no backbone atoms it was counted both as an extra chain and as a residue missing
backbone atoms. Waters and free ions are now excluded from the backbone report, which is a
report on the polypeptide; that exclusion, not the parser, is what most of the panel-wide change
below comes from.

**Corrected values.** Structure quality, panel-wide: templates with more than one chain 41 to 25,
templates missing backbone atoms away from the termini 13 to 2, templates with an unexplained
chain break 9 to 1. Site detection, on re-running all 129 templates: 31 of 192 rows changed, the
11 rows that reported 0.000 recovery fell to none, and mean best recovery rose from 0.810 to 0.889
(fpocket, 129 templates) and from 0.854 to 0.918 (DoGSite, 63). DoGSite's rank-1 count rose from
24 to 27 of 63 and its first-three count from 51 to 55. MTORC1 interaction recovery is now defined,
against seven reference contacts at the mTOR hinge and ATP site, residues 2185, 2190, 2237, 2239,
2240, 2245 and 2356 of chain B; the other four targets were re-run with the same code and their
mean recovery changed by at most 0.0002, which is one pose.

No conclusion is reversed, and one is strengthened: the detectors essentially always contain the
ligand site, so their weakness is ranking rather than detection. The earlier draft's statement
that MTORC1's best pocket held only about half of the lining residues was an artefact of this
defect; the figure is 0.92 for fpocket and 0.83 for DoGSite.

**Guards added.** A mol2 whose `SUBSTRUCTURE` records yield no residue at all is now refused
rather than written with every residue unknown, and an empty interaction reference is refused
rather than scored as missing values for every compound. Both were silent before. The parser is
one shared function with tests for both layouts, for the `****` unusual-residue marker in each,
and against the real 4jsx file.

---

## Amendment 4 (2026-11-12, after the variant comparison had been computed): multiplicity

Written after the comparison it governs, which Amendments 1 to 3 were not. It is recorded that
way deliberately, because the change can only tighten a conclusion and never create one.

**What was wrong with the reporting.** Amendment 3 compares twelve score variants against the
registered consensus at each metric and negative set, and counts an improvement where the paired
bootstrap interval excludes zero. Twelve simultaneous comparisons read at a pointwise 95% level
produce a spurious winner at roughly the nominal rate per family; the rule as written had no
control for that.

**The change.** Each comparison now carries a second interval for the difference, Bonferroni
widened to 1 - 0.05/m over the m variants compared within that metric and negative set, computed
from the same bootstrap draws so it costs no extra sampling. Promotion reads the adjusted
interval. The pointwise interval is still reported, because it describes the single comparison a
reader may want, and the two are distinguished in every table.

**Why this cannot have been chosen to suit the result.** The adjusted interval contains the
pointwise one by construction, so a variant that fails the pointwise test fails the adjusted test
as well. The pre-amendment result was already that no variant beat the registered consensus on
both hard negative sets. The amendment therefore cannot promote anything that was not promoted
before, and it cannot demote anything that was, because nothing was. It changes what can be
claimed in future runs, and here it only makes an existing null result harder to dispute.

**Resolution.** An adjusted interval reads a far quantile: at 95% over twelve variants each tail
holds 0.2% of the draws, two replicates in a thousand. The variant comparison is therefore run at
5,000 replicates, about twenty per tail, while every other table keeps the 1,000 replicates fixed
in the original protocol, which need no such depth. Where the draws are too few to resolve the
adjusted interval the tables say so in an `adj_resolved` column, and no variant may be promoted on
an unresolved one.

**Scaffold dependence, disclosed not fixed.** Every bootstrap here resamples compounds, which
treats them as independent. Confirmed actives for one target are often analogues, so the true
information content is below the nominal count and the intervals are, if anything, narrower than a
scaffold-aware resampling would give. The bias favours finding a difference, so it weakens a
positive claim and not a null one. No scaffold-level resampling is implemented; this is stated in
the paper's limitations rather than corrected.

---

## Amendment 5 (2026-11-13, after the performance study had been computed): the baseline comparison

**What the original protocol said.** Single properties were to be scored as if they were docking
scores on every negative set, and a docking result was to be described as exceeding a property
baseline only where the paired comparison supported it. The maximum ECFP4 Tanimoto similarity of
each compound to any active other than itself was listed as an analogue baseline and reported as a
per-target table.

**What was unspecified.** The protocol named the comparison but not its multiplicity treatment,
and it did not say whether the analogue baseline entered the same comparison as the physical
descriptors or stood apart as a descriptive table.

**The change.** Docking is now compared against the best of all eight structure-free baselines,
the seven descriptors and the similarity score together, in the same target and negative set, by
paired bootstrap on the same compounds at 5,000 replicates. The difference carries both a
pointwise and a Bonferroni-widened interval over the eight baselines, under the rule fixed in
Amendment 4. Taking the best of eight is the conservative direction for a claim that docking
exceeds a baseline and the anti-conservative direction for a claim that it falls below one, so
the adjusted interval governs both readings and both are reported.

**Why this cannot have been chosen to suit the result.** Docking did not exceed the best baseline
in any of the 45 target, score and negative-set cells at either interval. Restricting the
comparison to the seven physical descriptors, which is the narrower reading the original protocol
permits, would have left the similarity baseline out of the family and made the adjusted interval
narrower, which would have strengthened rather than weakened the reported conclusion that docking
falls below the baseline in a majority of cells on the matched designs. The wider family chosen
here is the less favourable one for that conclusion.

---

## Amendment 6 (2026-11-21, after every analysis above had been computed): confirmatory arm at full active-set size

**What the original protocol said.** Forty actives were to be drawn per target, or all actives
where fewer than forty exist. That cap was chosen for compute cost and fixes the width of every
interval the registered design can report.

**What was unspecified.** The protocol did not say whether the panel would also be screened
without the cap. The question it leaves open was stated in the manuscript's conclusions before
this arm was run: whether the matched-design result is a property of the panel or only of the
interval width that forty actives permit.

**The change.** The whole screen was repeated with the cap removed, admitting every available
active: 102, 308, 97, 27 and 79 for ESR1_ant, MAPK1, MTORC1, PPARG and TP53, which is 613 actives
against the 187 sampled and 11,921 compounds against 4,679. The negative-set constructions, the
matching ratio, the templates, the engines, the scores and the random seed are unchanged, so the
arms differ in active-set size alone.

**What this amendment can do to the reported result, stated plainly.** Unlike Amendments 4 and 5,
this one is not constrained to withhold a claim. More actives narrow every interval, and a
narrower interval can resolve an effect the registered arm could not. It did: on the two matched
designs the mean interval width falls from 0.200 to 0.134 and the number of registered cells whose
interval excludes 0.5 rises from 1 of 30 to 10 of 30, against 1.5 expected by chance. The
registered arm's reading that the matched-design signal is indistinguishable from chance was
therefore a limit of its resolution, and the confirmatory arm corrects it. The effect it resolves
is small: the macro-average ROC-AUC moves by at most 0.016 in any score and negative-set
combination, so what changes is what can be distinguished from chance and not how well docking
ranks these compounds.

The central comparison moves the other way. Docking exceeds the best of the eight structure-free
baselines in 0 of 45 cells in both arms, and falls resolvably below it in 25 of the 30
matched-design cells against 9 of 30 in the registered arm, with a lower point estimate in all 30
against 27 of 30. The conclusion that docking does not outrank a ligand-similarity baseline on
property-matched confirmed inactives is therefore strengthened by this arm, while the weaker
secondary claim about chance is revised. Both are reported in Section 3.5 of the manuscript, and
the registered arm remains the primary analysis.
