**Table S1. Pose selection**

| target | negative set | AUC lowest energy | AUC largest cluster | difference [95% CI] |
|---|---|---|---|---|
| ESR1_ant | random | 0.65 | 0.65 | -0.003 [-0.015, +0.008] |
| ESR1_ant | matched | 0.52 | 0.52 | -0.004 [-0.015, +0.006] |
| ESR1_ant | decoy | 0.43 | 0.42 | -0.003 [-0.016, +0.009] |
| MAPK1 | random | 0.57 | 0.58 | +0.009 [-0.004, +0.021] |
| MAPK1 | matched | 0.51 | 0.52 | +0.008 [-0.002, +0.017] |
| MAPK1 | decoy | 0.50 | 0.51 | +0.004 [-0.006, +0.013] |
| MTORC1 | random | 0.49 | 0.48 | -0.018 [-0.053, +0.019] |
| MTORC1 | matched | 0.50 | 0.47 | -0.026 [-0.064, +0.005] |
| MTORC1 | decoy | 0.59 | 0.57 | -0.019 [-0.055, +0.014] |
| PPARG | random | 0.77 | 0.78 | +0.003 [-0.011, +0.015] |
| PPARG | matched | 0.58 | 0.58 | -0.001 [-0.030, +0.021] |
| PPARG | decoy | 0.49 | 0.49 | +0.000 [-0.031, +0.025] |
| TP53 | random | 0.59 | 0.58 | -0.004 [-0.032, +0.024] |
| TP53 | matched | 0.56 | 0.56 | -0.005 [-0.027, +0.015] |
| TP53 | decoy | 0.42 | 0.43 | +0.005 [-0.017, +0.025] |

**Table S2. Throughput**

| target | compounds | ligand prep (s/compound) | AutoDock-GPU (s/compound) | gnina (s/compound) |
|---|---|---|---|---|
| ESR1_ant | 2034 | 0.004 | 0.85 | 4.03 |
| MAPK1 | 5536 | 0.008 | 3.99 | 9.61 |
| MTORC1 | 1949 | 0.008 | 10.19 | 9.31 |
| PPARG | 759 | 0.004 | 6.43 | 7.60 |
| TP53 | 1643 | 0.006 | 8.10 | 5.71 |

**Table S3. Failures**

| target | class | n | prep failed | AutoDock-GPU failed | gnina failed | gnina positive energy | chelator needs manual ionisation | chelator deprotonated |
|---|---|---|---|---|---|---|---|---|
| ESR1_ant | active | 102 | 0 | 0 | 0 | 7% | 0 | 0 |
| ESR1_ant | inactive | 1932 | 0 | 0 | 0 | 4% | 22 | 0 |
| MAPK1 | active | 308 | 0 | 0 | 0 | 0% | 3 | 0 |
| MAPK1 | inactive | 5228 | 6 | 0 | 0 | 0% | 28 | 0 |
| MTORC1 | active | 97 | 0 | 0 | 0 | 0% | 2 | 0 |
| MTORC1 | inactive | 1852 | 0 | 0 | 0 | 0% | 2 | 0 |
| PPARG | active | 27 | 0 | 0 | 0 | 0% | 0 | 0 |
| PPARG | inactive | 732 | 0 | 0 | 0 | 1% | 14 | 0 |
| TP53 | active | 79 | 0 | 0 | 0 | 0% | 0 | 0 |
| TP53 | inactive | 1564 | 1 | 0 | 0 | 0% | 33 | 0 |

**Table S4. Metal-ion audit of the templates against their original depositions**

| target | templates | min frame overlap | ions in benchmark files | templates with a pocket ion | pocket ion (nearest) | nearest functional ion within 30 A |
|---|---|---|---|---|---|---|
| ADRB2 | 8 | 1.000 | 0 | 0 |  | MG 10.1 A |
| ALDH1 | 8 | 0.998 | 0 | 0 |  |  |
| ESR1_ago | 15 | 0.996 | 0 | 0 |  |  |
| ESR1_ant | 15 | 1.000 | 0 | 0 |  |  |
| FEN1 | 1 | 1.000 | 0 | 1 | MG 2.3 A | MG 2.3 A |
| GBA | 6 | 1.000 | 0 | 0 |  |  |
| IDH1 | 14 | 1.000 | 0 | 0 |  | MG 19.2 A |
| KAT2A | 3 | 1.000 | 0 | 0 |  |  |
| MAPK1 | 15 | 0.997 | 0 | 0 |  |  |
| MTORC1 | 11 | 1.000 | 0 | 0 |  |  |
| OPRK1 | 1 | 1.000 | 0 | 0 |  |  |
| PKM2 | 9 | 1.000 | 0 | 0 |  | MG 23.2 A |
| PPARG | 15 | 1.000 | 0 | 0 |  |  |
| TP53 | 6 | 1.000 | 0 | 0 |  | ZN 29.0 A |
| VDR | 2 | 1.000 | 0 | 0 |  |  |

**Table S5. Change in ROC-AUC between negative-set designs**

| negatives A | negatives B | score | mean_delta_auc | targets | targets_with_ci_excluding_zero |
|---|---|---|---|---|---|
| random | matched | autodock-gpu | +0.083 | 5 | 3 |
| random | matched | gnina | +0.056 | 5 | 2 |
| random | matched | consensus | +0.091 | 5 | 4 |
| matched | decoy | autodock-gpu | +0.047 | 5 | 4 |
| matched | decoy | gnina | +0.010 | 5 | 3 |
| matched | decoy | consensus | +0.026 | 5 | 4 |
| random | decoy | autodock-gpu | +0.130 | 5 | 5 |
| random | decoy | gnina | +0.065 | 5 | 4 |
| random | decoy | consensus | +0.117 | 5 | 5 |

**Table S6. Seed stability of AutoDock-GPU**

| target | n_seeds | auc_min | auc_max | auc_range | mean_spearman_between_seeds | stable_top_n_fraction | mean_per_compound_sd |
|---|---|---|---|---|---|---|---|
| ESR1_ant | 3 | 0.647 | 0.654 | 0.008 | 0.995 | 0.900 | 0.043 |
| MAPK1 | 3 | 0.570 | 0.573 | 0.003 | 0.984 | 0.600 | 0.086 |
| MTORC1 | 3 | 0.490 | 0.498 | 0.008 | 0.938 | 0.700 | 0.157 |
| PPARG | 3 | 0.771 | 0.775 | 0.004 | 0.997 | 0.950 | 0.049 |
| TP53 | 3 | 0.585 | 0.592 | 0.007 | 0.981 | 0.750 | 0.076 |

**Table S7. Similarity of the negative sets to the actives**

| target | actives | median NN similarity among actives | clusters at 0.4 | analogue AUC vs random | analogue AUC vs matched | analogue AUC vs decoy |
|---|---|---|---|---|---|---|
| ESR1_ant | 102 | 0.29 | 73 | 0.71 | 0.67 | 0.77 |
| MAPK1 | 308 | 0.37 | 222 | 0.70 | 0.64 | 0.77 |
| MTORC1 | 97 | 0.41 | 64 | 0.74 | 0.68 | 0.85 |
| PPARG | 27 | 0.29 | 21 | 0.75 | 0.66 | 0.68 |
| TP53 | 79 | 0.50 | 49 | 0.73 | 0.72 | 0.77 |

**Table S8. Median heavy-atom gap, actives against each negative set**

| target | decoy | matched | random |
|---|---|---|---|
| ESR1_ant | -1.0 | +1.0 | +5.5 |
| MAPK1 | +1.0 | +0.5 | +2.0 |
| MTORC1 | +0.0 | +0.0 | +0.0 |
| PPARG | +0.0 | +2.0 | +7.0 |
| TP53 | -1.0 | +1.0 | +2.0 |

**Table S9. ADMET risk tiers**

| target | actives PASS | actives MODERATE | actives HIGH_RISK | actives EXCLUDE | negatives PASS | negatives MODERATE | negatives HIGH_RISK | negatives EXCLUDE |
|---|---|---|---|---|---|---|---|---|
| ESR1_ant | 42% | 1% | 35% | 22% | 52% | 2% | 28% | 19% |
| MAPK1 | 22% | 0% | 43% | 35% | 35% | 0% | 45% | 20% |
| MTORC1 | 16% | 0% | 61% | 23% | 19% | 0% | 65% | 15% |
| PPARG | 67% | 11% | 7% | 15% | 68% | 3% | 14% | 15% |
| TP53 | 25% | 4% | 9% | 62% | 60% | 2% | 13% | 25% |

**Table S10. Effect of ADMET filters on retrieval**

| negative_set | filter | actives_retained | inactives_retained | actives_lost | auc_before | auc_after | targets |
|---|---|---|---|---|---|---|---|
| random | not EXCLUDE (no PAINS, Ames <= 0.5) | 69% | 82% | 205 | 0.643 | 0.653 | 5 |
| random | not HIGH_RISK/EXCLUDE (also hERG <= 0.5) | 38% | 53% | 442 | 0.643 | 0.659 | 5 |
| random | Veber pass | 95% | 96% | 17 | 0.643 | 0.649 | 5 |
| random | no PAINS alert | 92% | 96% | 50 | 0.643 | 0.647 | 5 |
| matched | not EXCLUDE (no PAINS, Ames <= 0.5) | 69% | 82% | 205 | 0.552 | 0.566 | 5 |
| matched | not HIGH_RISK/EXCLUDE (also hERG <= 0.5) | 38% | 50% | 442 | 0.552 | 0.554 | 5 |
| matched | Veber pass | 95% | 97% | 17 | 0.552 | 0.555 | 5 |
| matched | no PAINS alert | 92% | 95% | 50 | 0.552 | 0.557 | 5 |
| decoy | not EXCLUDE (no PAINS, Ames <= 0.5) | 69% | 81% | 205 | 0.526 | 0.545 | 5 |
| decoy | not HIGH_RISK/EXCLUDE (also hERG <= 0.5) | 38% | 46% | 442 | 0.526 | 0.553 | 5 |
| decoy | Veber pass | 95% | 98% | 17 | 0.526 | 0.529 | 5 |
| decoy | no PAINS alert | 92% | 94% | 50 | 0.526 | 0.531 | 5 |

**Table S11. Hit rate of the lead-selection strategies**

| negative_set | strategy | top_n | hit_rate | sd_across_targets | pool_active_fraction | targets |
|---|---|---|---|---|---|---|
| random | docking | 10 | 0.480 | 0.415 | 0.380 | 5 |
| random | docking | 20 | 0.430 | 0.244 | 0.380 | 5 |
| random | docking | 50 | 0.400 | 0.223 | 0.380 | 5 |
| random | docking+ADMET | 10 | 0.280 | 0.192 | 0.380 | 5 |
| random | docking+ADMET | 20 | 0.315 | 0.165 | 0.380 | 5 |
| random | docking+ADMET | 50 | 0.277 | 0.070 | 0.380 | 5 |
| random | docking+IFP | 10 | 0.520 | 0.327 | 0.380 | 5 |
| random | docking+IFP | 20 | 0.380 | 0.182 | 0.380 | 5 |
| random | docking+IFP | 50 | 0.400 | 0.176 | 0.380 | 5 |
| random | pareto | 10 | 0.400 | 0.245 | 0.380 | 5 |
| random | pareto | 20 | 0.360 | 0.195 | 0.380 | 5 |
| random | pareto | 50 | 0.376 | 0.197 | 0.380 | 5 |
| matched | docking | 10 | 0.160 | 0.114 | 0.122 | 5 |
| matched | docking | 20 | 0.130 | 0.084 | 0.122 | 5 |
| matched | docking | 50 | 0.144 | 0.059 | 0.122 | 5 |
| matched | docking+ADMET | 10 | 0.060 | 0.089 | 0.122 | 5 |
| matched | docking+ADMET | 20 | 0.081 | 0.076 | 0.122 | 5 |
| matched | docking+ADMET | 50 | 0.083 | 0.069 | 0.122 | 5 |
| matched | docking+IFP | 10 | 0.080 | 0.084 | 0.122 | 5 |
| matched | docking+IFP | 20 | 0.140 | 0.065 | 0.122 | 5 |
| matched | docking+IFP | 50 | 0.144 | 0.061 | 0.122 | 5 |
| matched | pareto | 10 | 0.120 | 0.110 | 0.122 | 5 |
| matched | pareto | 20 | 0.120 | 0.067 | 0.122 | 5 |
| matched | pareto | 50 | 0.136 | 0.062 | 0.122 | 5 |
| decoy | docking | 10 | 0.160 | 0.114 | 0.110 | 5 |
| decoy | docking | 20 | 0.130 | 0.084 | 0.110 | 5 |
| decoy | docking | 50 | 0.128 | 0.086 | 0.110 | 5 |
| decoy | docking+ADMET | 10 | 0.080 | 0.130 | 0.110 | 5 |
| decoy | docking+ADMET | 20 | 0.116 | 0.112 | 0.110 | 5 |
| decoy | docking+ADMET | 50 | 0.118 | 0.092 | 0.110 | 5 |
| decoy | docking+IFP | 10 | 0.100 | 0.100 | 0.110 | 5 |
| decoy | docking+IFP | 20 | 0.140 | 0.089 | 0.110 | 5 |
| decoy | docking+IFP | 50 | 0.132 | 0.064 | 0.110 | 5 |
| decoy | pareto | 10 | 0.140 | 0.089 | 0.110 | 5 |
| decoy | pareto | 20 | 0.110 | 0.065 | 0.110 | 5 |
| decoy | pareto | 50 | 0.128 | 0.072 | 0.110 | 5 |

**Table S12. Recovery of the co-crystal interactions by the top gnina pose**

| target | crystal contacts | recall, actives | n, actives | recall, top-ranked negatives | n, top-ranked negatives |
|---|---|---|---|---|---|
| ESR1_ant | 16 | 0.23 | 102 | 0.25 | 1932 |
| MAPK1 | 8 | 0.41 | 308 | 0.40 | 5222 |
| MTORC1 | 7 | 0.46 | 97 | 0.46 | 1852 |
| PPARG | 9 | 0.44 | 27 | 0.37 | 732 |
| TP53 | 5 | 0.23 | 79 | 0.23 | 1563 |

**Table S13. Structure quality of the templates**

| target | templates | median_resolution | with_chain_break | multiple_chains | engineered_mutation | median_rama_core |
|---|---|---|---|---|---|---|
| ADRB2 | 8 | 3.2 | 0 | 3 | 8 | 0.815 |
| ALDH1 | 8 | 1.85 | 0 | 0 | 3 | 0.749 |
| ESR1_ago | 15 | 1.95 | 0 | 0 | 14 | 0.84 |
| ESR1_ant | 15 | 1.9 | 0 | 1 | 9 | 0.821 |
| FEN1 | 1 | 2.84 | 0 | 0 | 0 | 0.775 |
| GBA | 6 | 2.355 | 0 | 0 | 0 | 0.652 |
| IDH1 | 14 | 2.7 | 0 | 1 | 14 | 0.77 |
| KAT2A | 3 | 2.0 | 0 | 0 | 0 | 0.796 |
| MAPK1 | 15 | 2.0 | 1 | 2 | 1 | 0.698 |
| MTORC1 | 11 | 2.2 | 0 | 10 | 0 | 0.784 |
| OPRK1 | 1 | 3.1 | 0 | 0 | 1 | 0.886 |
| PKM2 | 9 | 2.1 | 0 | 8 | 0 | 0.768 |
| PPARG | 15 | 2.05 | 0 | 0 | 0 | 0.826 |
| TP53 | 6 | 1.435 | 0 | 0 | 6 | 0.628 |
| VDR | 2 | 2.985 | 0 | 0 | 2 | 0.82 |

**Table S14. Pocket-detector recovery of the ligand site**

| detector | templates | top-1 recovery | top-3 recovery | mean best-pocket recovery | median pockets |
|---|---|---|---|---|---|
| fpocket | 129 | 37/129 (29%) [22%, 37%] | 61/129 (47%) [39%, 56%] | 0.89 | 22 |
| dogsite | 63 | 27/63 (43%) [31%, 55%] | 55/63 (87%) [77%, 93%] | 0.92 | 13 |
| p2rank | 129 | 99/129 (77%) [69%, 83%] | 123/129 (95%) [90%, 98%] | 0.87 | 8 |

**Table S15. Druggability controls**

| template | control_role | description | method | rank_of_site | site_pocket_druggability | ligand_buried_fraction |
|---|---|---|---|---|---|---|
| 1HXB | positive | HIV-1 protease with saquinavir | fpocket | 4 | 0.06 | 0.57 |
| 1HXB | positive | HIV-1 protease with saquinavir | dogsite | 1 | 0.81 | 0.57 |
| 1M17 | positive | EGFR kinase with erlotinib | fpocket | 2 | 0.32 | 0.89 |
| 1M17 | positive | EGFR kinase with erlotinib | dogsite | 1 | 0.81 | 0.89 |
| 1CIL | positive | carbonic anhydrase II with a sulfonamide (zinc enzyme) | fpocket | 8 | 0.34 | 0.73 |
| 1CIL | positive | carbonic anhydrase II with a sulfonamide (zinc enzyme) | dogsite | 1 | 0.64 | 0.73 |
| 1BZC | negative | PTP1B with an inhibitor: a polar, historically difficult pocket | fpocket | 1 | 0.29 | 0.93 |
| 1BZC | negative | PTP1B with an inhibitor: a polar, historically difficult pocket | dogsite | 1 | 0.74 | 0.93 |
| 1HXB | positive | HIV-1 protease with saquinavir | p2rank | 1 | 0.89 | 0.57 |
| 1M17 | positive | EGFR kinase with erlotinib | p2rank | 1 | 0.56 | 0.89 |
| 1CIL | positive | carbonic anhydrase II with a sulfonamide (zinc enzyme) | p2rank | 1 | 0.93 | 0.73 |
| 1BZC | negative | PTP1B with an inhibitor: a polar, historically difficult pocket | p2rank | 1 | 0.43 | 0.93 |

**Table S16. Quality-gate ablation**

| condition | n | mean_dG | median_dG | min_dG | max_dG | non_binding_dG_ge_0 | dG_below_minus20 | n_charge_flagged | n_torsion_frozen_flagged | mean_shift_vs_baseline | shift_ci_low | shift_ci_high | target | template | n_ligands |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| baseline | 60 | -7.617 | -7.625 | -12.71 | -3.72 | 0 | 0 | 0 | 0 | 0.0 | 0.0 | 0.0 | ESR1_ant | 6b0f | 60 |
| charges_x4_uncapped | 60 | -3.588 | -3.38 | -9.85 | 4.5 | 4 | 0 | 58 | 0 | 4.028 | 3.291 | 4.811 | ESR1_ant | 6b0f | 60 |
| charges_x4_capped | 60 | -3.789 | -4.11 | -8.09 | 3.02 | 2 | 0 | 0 | 0 | 3.828 | 3.152 | 4.531 | ESR1_ant | 6b0f | 60 |
| torsions_frozen | 60 | -7.397 | -7.53 | -12.86 | -3.67 | 0 | 0 | 0 | 40 | 0.22 | -0.051 | 0.493 | ESR1_ant | 6b0f | 60 |
| baseline | 60 | -8.17 | -8.345 | -9.81 | -5.35 | 0 | 0 | 0 | 0 | 0.0 | 0.0 | 0.0 | MAPK1 | 4qte | 60 |
| charges_x4_uncapped | 60 | -3.515 | -3.74 | -5.82 | -1.54 | 0 | 0 | 59 | 0 | 4.655 | 4.349 | 4.97 | MAPK1 | 4qte | 60 |
| charges_x4_capped | 60 | -3.652 | -3.915 | -5.83 | -1.21 | 0 | 0 | 0 | 0 | 4.518 | 4.157 | 4.88 | MAPK1 | 4qte | 60 |
| torsions_frozen | 60 | -7.193 | -7.27 | -9.27 | -5.02 | 0 | 0 | 0 | 50 | 0.977 | 0.748 | 1.203 | MAPK1 | 4qte | 60 |
| baseline | 60 | -9.016 | -9.075 | -11.07 | -6.83 | 0 | 0 | 0 | 0 | 0.0 | 0.0 | 0.0 | MTORC1 | 4jsx | 60 |
| charges_x4_uncapped | 60 | -3.68 | -3.275 | -11.73 | -0.57 | 0 | 0 | 60 | 0 | 5.336 | 4.916 | 5.72 | MTORC1 | 4jsx | 60 |
| charges_x4_capped | 60 | -4.999 | -4.3 | -10.21 | -1.42 | 0 | 0 | 0 | 0 | 4.017 | 3.558 | 4.47 | MTORC1 | 4jsx | 60 |
| torsions_frozen | 60 | -8.017 | -7.98 | -10.17 | -6.47 | 0 | 0 | 0 | 58 | 0.999 | 0.764 | 1.246 | MTORC1 | 4jsx | 60 |
| baseline | 57 | -7.777 | -7.98 | -10.66 | -4.59 | 0 | 0 | 0 | 0 | 0.0 | 0.0 | 0.0 | PPARG | 2p4y | 57 |
| charges_x4_uncapped | 57 | -4.031 | -4.87 | -8.22 | 0.77 | 2 | 0 | 57 | 0 | 3.746 | 3.146 | 4.35 | PPARG | 2p4y | 57 |
| charges_x4_capped | 57 | -4.086 | -4.25 | -7.54 | 0.39 | 1 | 0 | 0 | 0 | 3.691 | 3.238 | 4.158 | PPARG | 2p4y | 57 |
| torsions_frozen | 57 | -7.91 | -8.06 | -10.02 | -5.51 | 0 | 0 | 0 | 51 | -0.133 | -0.317 | 0.068 | PPARG | 2p4y | 57 |
| baseline | 60 | -6.241 | -6.215 | -9.83 | -4.08 | 0 | 0 | 0 | 0 | 0.0 | 0.0 | 0.0 | TP53 | 4agq | 60 |
| charges_x4_uncapped | 60 | -2.672 | -2.68 | -7.85 | 1.44 | 6 | 0 | 57 | 0 | 3.57 | 3.067 | 4.019 | TP53 | 4agq | 60 |
| charges_x4_capped | 60 | -3.474 | -3.55 | -7.73 | -0.13 | 0 | 0 | 0 | 0 | 2.767 | 2.422 | 3.101 | TP53 | 4agq | 60 |
| torsions_frozen | 60 | -6.437 | -6.5 | -8.35 | -4.24 | 0 | 0 | 0 | 28 | -0.196 | -0.343 | -0.039 | TP53 | 4agq | 60 |

**Table S17. Template difficulty by target**

| target | templates | redocked | validated (<=2 A) | macrocyclic ligands | median heavy atoms | median MW | best RMSD (A) |
|---|---|---|---|---|---|---|---|
| ESR1_ant | 15 | 13 | 10 | 0 | 33 | 428 | 0.48 |
| MAPK1 | 15 | 7 | 4 | 1 | 26 | 377 | 0.60 |
| MTORC1 | 11 | 10 | 0 | 9 | 65 | 914 | 6.21 |
| PPARG | 15 | 10 | 2 | 0 | 32 | 439 | 0.96 |
| TP53 | 6 | 4 | 1 | 0 | 22 | 330 | 1.78 |

**Table S18. Which co-crystal ligand property predicts redocking failure**

| ligand property | Spearman rho vs RMSD | p (permutation) | templates |
|---|---|---|---|
| largest ring size | +0.490 | 0.002 | 42 |
| molecular weight | +0.320 | 0.049 | 42 |
| heavy atoms | +0.299 | 0.060 | 42 |
| rotatable bonds | +0.207 | 0.190 | 42 |
| number of rings | -0.045 | 0.778 | 42 |

**Table S19. Redocking validation rate by macrocyclicity**

| template group | templates | validated | validation rate | 95% CI |
|---|---|---|---|---|
| non-macrocyclic | 42 | 16 | 38.1% | 25.0% to 53.2% |
| macrocyclic | 10 | 0 | 0.0% | 0.0% to 27.8% |
| ligand not sanitizable | 10 | 1 | 10.0% |  |

**Table S20. Score variants by ROC-AUC against each matched design**

| negative set | score | roc_auc | vs registered consensus | targets better | targets worse | better unadjusted | targets |
|---|---|---|---|---|---|---|---|
| property-matched inactives | fusion: reciprocal rank | 0.559 | +0.008 | 1 | 0 | 1 | 5 |
| property-matched inactives | gnina alone | 0.556 | +0.004 | 0 | 0 | 0 | 5 |
| property-matched inactives | fusion: redocking-weighted z | 0.553 | +0.002 | 0 | 0 | 0 | 5 |
| property-matched inactives | fusion: rank mean | 0.552 | +0.000 | 0 | 0 | 0 | 5 |
| property-matched inactives | gnina (CNN score) | 0.540 | -0.012 | 0 | 0 | 0 | 5 |
| property-matched inactives | gnina (CNN affinity) | 0.535 | -0.016 | 0 | 0 | 0 | 5 |
| property-matched inactives | gnina (Vina-type) per heavy atom^(1/3) | 0.533 | -0.018 | 0 | 0 | 0 | 5 |
| property-matched inactives | AutoDock-GPU alone | 0.533 | -0.019 | 0 | 0 | 0 | 5 |
| property-matched inactives | AutoDock-GPU per heavy atom^(1/3) | 0.528 | -0.023 | 0 | 0 | 0 | 5 |
| property-matched inactives | interaction recall | 0.506 | -0.045 | 0 | 0 | 0 | 5 |
| property-matched inactives | gnina (Vina-type) per heavy atom | 0.500 | -0.051 | 0 | 0 | 0 | 5 |
| property-matched inactives | AutoDock-GPU per heavy atom | 0.493 | -0.059 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | gnina (CNN score) | 0.562 | +0.036 | 0 | 0 | 1 | 5 |
| DUD-E-protocol decoys | gnina alone | 0.546 | +0.020 | 0 | 0 | 1 | 5 |
| DUD-E-protocol decoys | fusion: reciprocal rank | 0.535 | +0.009 | 1 | 0 | 1 | 5 |
| DUD-E-protocol decoys | fusion: redocking-weighted z | 0.529 | +0.003 | 0 | 0 | 1 | 5 |
| DUD-E-protocol decoys | gnina (Vina-type) per heavy atom | 0.526 | -0.000 | 0 | 0 | 1 | 5 |
| DUD-E-protocol decoys | fusion: rank mean | 0.523 | -0.003 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | gnina (Vina-type) per heavy atom^(1/3) | 0.516 | -0.010 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | gnina (CNN affinity) | 0.514 | -0.012 | 1 | 0 | 1 | 5 |
| DUD-E-protocol decoys | AutoDock-GPU per heavy atom | 0.504 | -0.022 | 0 | 1 | 0 | 5 |
| DUD-E-protocol decoys | interaction recall | 0.493 | -0.033 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | AutoDock-GPU per heavy atom^(1/3) | 0.486 | -0.040 | 0 | 2 | 0 | 5 |
| DUD-E-protocol decoys | AutoDock-GPU alone | 0.486 | -0.040 | 0 | 0 | 0 | 5 |

**Table S21. Score variants by EF1 against each matched design**

| negative set | score | ef1 | vs registered consensus | targets better | targets worse | better unadjusted | targets |
|---|---|---|---|---|---|---|---|
| property-matched inactives | gnina (CNN score) | 2.43 | +0.67 | 0 | 0 | 0 | 5 |
| property-matched inactives | gnina alone | 1.93 | +0.17 | 0 | 0 | 0 | 5 |
| property-matched inactives | fusion: redocking-weighted z | 1.77 | +0.00 | 0 | 0 | 0 | 5 |
| property-matched inactives | fusion: reciprocal rank | 1.68 | -0.08 | 0 | 0 | 0 | 5 |
| property-matched inactives | fusion: rank mean | 1.39 | -0.38 | 0 | 0 | 0 | 5 |
| property-matched inactives | AutoDock-GPU alone | 1.31 | -0.45 | 0 | 0 | 0 | 5 |
| property-matched inactives | gnina (CNN affinity) | 1.28 | -0.48 | 0 | 0 | 0 | 5 |
| property-matched inactives | AutoDock-GPU per heavy atom^(1/3) | 1.17 | -0.59 | 0 | 0 | 0 | 5 |
| property-matched inactives | AutoDock-GPU per heavy atom | 1.06 | -0.70 | 0 | 0 | 0 | 5 |
| property-matched inactives | gnina (Vina-type) per heavy atom | 0.86 | -0.91 | 0 | 0 | 0 | 5 |
| property-matched inactives | gnina (Vina-type) per heavy atom^(1/3) | 0.77 | -1.00 | 0 | 0 | 0 | 5 |
| property-matched inactives | interaction recall | 0.63 | -1.14 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | gnina (CNN score) | 2.41 | +1.25 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | AutoDock-GPU per heavy atom | 1.64 | +0.48 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | gnina (Vina-type) per heavy atom | 1.57 | +0.41 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | AutoDock-GPU per heavy atom^(1/3) | 1.37 | +0.21 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | AutoDock-GPU alone | 1.22 | +0.06 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | fusion: redocking-weighted z | 1.17 | +0.00 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | fusion: reciprocal rank | 1.17 | +0.00 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | gnina alone | 1.11 | -0.06 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | fusion: rank mean | 1.03 | -0.14 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | gnina (CNN affinity) | 0.88 | -0.28 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | gnina (Vina-type) per heavy atom^(1/3) | 0.77 | -0.40 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | interaction recall | 0.65 | -0.51 | 0 | 0 | 0 | 5 |
