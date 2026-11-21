**Table 1. Panel, template validation and compounds docked**

| target | actives | random inactives | matched inactives | decoys | templates validated | pocket ions | deposition checked | selected | RMSD AutoDock-GPU (A) | RMSD gnina (A) | validated |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ESR1_ant | 40 | 300 | 320 | 320 | 10/15 | none | yes | 6b0f | 0.48 | 0.40 | yes |
| MAPK1 | 40 | 300 | 320 | 320 | 4/15 | none | yes | 4qte | 0.60 | 0.44 | yes |
| MTORC1 | 40 | 300 | 320 | 320 | 0/11 | none | yes | 4jsx | 6.13 | 6.21 | NO |
| PPARG | 27 | 300 | 216 | 216 | 2/15 | none | yes | 2p4y | 0.96 | 0.45 | yes |
| TP53 | 40 | 300 | 320 | 320 | 1/6 | none | yes | 4agq | 1.78 | 1.37 | yes |

**Table 2. ROC-AUC against each negative set**

| negative set | target | autodock-gpu | gnina | consensus |
|---|---|---|---|---|
| random inactives | ESR1_ant | 0.68 [0.59, 0.77] | 0.61 [0.50, 0.71] | 0.68 [0.58, 0.77] |
| random inactives | MAPK1 | 0.56 [0.47, 0.66] | 0.58 [0.48, 0.68] | 0.59 [0.49, 0.70] |
| random inactives | MTORC1 | 0.53 [0.43, 0.63] | 0.47 [0.39, 0.56] | 0.51 [0.42, 0.59] |
| random inactives | PPARG | 0.77 [0.69, 0.84] | 0.78 [0.68, 0.88] | 0.81 [0.72, 0.88] |
| random inactives | TP53 | 0.60 [0.49, 0.70] | 0.56 [0.46, 0.65] | 0.62 [0.52, 0.73] |
| property-matched inactives | ESR1_ant | 0.49 [0.39, 0.58] | 0.49 [0.39, 0.60] | 0.50 [0.40, 0.59] |
| property-matched inactives | MAPK1 | 0.51 [0.40, 0.61] | 0.56 [0.46, 0.66] | 0.56 [0.45, 0.65] |
| property-matched inactives | MTORC1 | 0.51 [0.42, 0.60] | 0.47 [0.38, 0.56] | 0.47 [0.39, 0.55] |
| property-matched inactives | PPARG | 0.58 [0.47, 0.69] | 0.62 [0.51, 0.73] | 0.60 [0.48, 0.71] |
| property-matched inactives | TP53 | 0.56 [0.46, 0.67] | 0.55 [0.46, 0.64] | 0.58 [0.49, 0.68] |
| DUD-E-protocol decoys | ESR1_ant | 0.43 [0.33, 0.52] | 0.54 [0.44, 0.64] | 0.49 [0.38, 0.59] |
| DUD-E-protocol decoys | MAPK1 | 0.50 [0.41, 0.61] | 0.57 [0.46, 0.67] | 0.56 [0.46, 0.66] |
| DUD-E-protocol decoys | MTORC1 | 0.59 [0.50, 0.67] | 0.53 [0.44, 0.63] | 0.59 [0.50, 0.67] |
| DUD-E-protocol decoys | PPARG | 0.48 [0.37, 0.60] | 0.55 [0.42, 0.66] | 0.52 [0.40, 0.64] |
| DUD-E-protocol decoys | TP53 | 0.47 [0.37, 0.57] | 0.51 [0.42, 0.60] | 0.50 [0.40, 0.59] |

**Table 3. ROC-AUC of single properties used as scores**

| target | decoy: mw | decoy: n_heavy | decoy: tpsa | matched: mw | matched: n_heavy | matched: tpsa | random: mw | random: n_heavy | random: tpsa |
|---|---|---|---|---|---|---|---|---|---|
| ESR1_ant | 0.5 | 0.45 | 0.5 | 0.54 | 0.52 | 0.51 | 0.7 | 0.69 | 0.53 |
| MAPK1 | 0.5 | 0.54 | 0.5 | 0.5 | 0.5 | 0.48 | 0.54 | 0.55 | 0.42 |
| MTORC1 | 0.5 | 0.48 | 0.5 | 0.5 | 0.47 | 0.48 | 0.5 | 0.49 | 0.44 |
| PPARG | 0.5 | 0.5 | 0.5 | 0.54 | 0.57 | 0.55 | 0.74 | 0.77 | 0.59 |
| TP53 | 0.49 | 0.47 | 0.51 | 0.53 | 0.54 | 0.52 | 0.6 | 0.59 | 0.58 |

**Table 4. Size-controlled ROC-AUC against each negative set**

| negative set | score | ROC-AUC | size-controlled AUC | change | SD of change | targets |
|---|---|---|---|---|---|---|
| random inactives | fusion: rank mean | 0.644 | 0.593 | -0.051 | 0.045 | 5 |
| random inactives | gnina | 0.600 | 0.593 | -0.008 | 0.032 | 5 |
| random inactives | fusion: reciprocal rank | 0.644 | 0.591 | -0.053 | 0.062 | 5 |
| random inactives | fusion: redocking-weighted z | 0.641 | 0.590 | -0.052 | 0.046 | 5 |
| random inactives | consensus | 0.642 | 0.587 | -0.054 | 0.049 | 5 |
| random inactives | gnina (CNN score) | 0.537 | 0.584 | +0.047 | 0.037 | 5 |
| random inactives | gnina (Vina-type) per heavy atom^(1/3) | 0.581 | 0.571 | -0.010 | 0.022 | 5 |
| random inactives | gnina (Vina-type) | 0.614 | 0.567 | -0.046 | 0.054 | 5 |
| random inactives | gnina (Vina-type) per heavy atom | 0.453 | 0.566 | +0.113 | 0.105 | 5 |
| random inactives | gnina (CNN affinity) | 0.619 | 0.561 | -0.057 | 0.059 | 5 |
| random inactives | AutoDock-GPU per heavy atom^(1/3) | 0.607 | 0.546 | -0.061 | 0.076 | 5 |
| random inactives | autodock-gpu | 0.630 | 0.545 | -0.085 | 0.103 | 5 |
| random inactives | AutoDock-GPU per heavy atom | 0.447 | 0.540 | +0.093 | 0.086 | 5 |
| random inactives | autodock-gpu (largest cluster) | 0.627 | 0.538 | -0.089 | 0.095 | 5 |
| random inactives | interaction recall | 0.520 | 0.496 | -0.024 | 0.025 | 5 |
| property-matched inactives | fusion: reciprocal rank | 0.545 | 0.554 | +0.009 | 0.020 | 5 |
| property-matched inactives | gnina | 0.539 | 0.542 | +0.003 | 0.013 | 5 |
| property-matched inactives | fusion: redocking-weighted z | 0.542 | 0.542 | -0.000 | 0.023 | 5 |
| property-matched inactives | fusion: rank mean | 0.540 | 0.540 | +0.000 | 0.020 | 5 |
| property-matched inactives | consensus | 0.541 | 0.540 | -0.001 | 0.023 | 5 |
| property-matched inactives | gnina (CNN score) | 0.526 | 0.535 | +0.008 | 0.008 | 5 |
| property-matched inactives | gnina (Vina-type) per heavy atom^(1/3) | 0.535 | 0.532 | -0.003 | 0.015 | 5 |
| property-matched inactives | gnina (Vina-type) | 0.533 | 0.528 | -0.005 | 0.026 | 5 |
| property-matched inactives | gnina (Vina-type) per heavy atom | 0.508 | 0.525 | +0.017 | 0.034 | 5 |
| property-matched inactives | AutoDock-GPU per heavy atom^(1/3) | 0.532 | 0.524 | -0.008 | 0.017 | 5 |
| property-matched inactives | gnina (CNN affinity) | 0.521 | 0.518 | -0.003 | 0.023 | 5 |
| property-matched inactives | autodock-gpu | 0.529 | 0.517 | -0.012 | 0.030 | 5 |
| property-matched inactives | AutoDock-GPU per heavy atom | 0.502 | 0.516 | +0.014 | 0.027 | 5 |
| property-matched inactives | autodock-gpu (largest cluster) | 0.522 | 0.508 | -0.014 | 0.027 | 5 |
| property-matched inactives | interaction recall | 0.485 | 0.485 | +0.000 | 0.011 | 5 |
| DUD-E-protocol decoys | gnina | 0.540 | 0.550 | +0.010 | 0.012 | 5 |
| DUD-E-protocol decoys | gnina (CNN score) | 0.546 | 0.550 | +0.004 | 0.013 | 5 |
| DUD-E-protocol decoys | fusion: reciprocal rank | 0.531 | 0.549 | +0.018 | 0.030 | 5 |
| DUD-E-protocol decoys | fusion: redocking-weighted z | 0.532 | 0.545 | +0.013 | 0.024 | 5 |
| DUD-E-protocol decoys | consensus | 0.529 | 0.542 | +0.013 | 0.024 | 5 |
| DUD-E-protocol decoys | fusion: rank mean | 0.526 | 0.540 | +0.014 | 0.026 | 5 |
| DUD-E-protocol decoys | gnina (CNN affinity) | 0.516 | 0.536 | +0.019 | 0.020 | 5 |
| DUD-E-protocol decoys | gnina (Vina-type) per heavy atom^(1/3) | 0.535 | 0.535 | +0.000 | 0.007 | 5 |
| DUD-E-protocol decoys | gnina (Vina-type) per heavy atom | 0.530 | 0.534 | +0.004 | 0.025 | 5 |
| DUD-E-protocol decoys | gnina (Vina-type) | 0.530 | 0.534 | +0.004 | 0.020 | 5 |
| DUD-E-protocol decoys | autodock-gpu | 0.493 | 0.507 | +0.014 | 0.022 | 5 |
| DUD-E-protocol decoys | AutoDock-GPU per heavy atom^(1/3) | 0.499 | 0.505 | +0.006 | 0.021 | 5 |
| DUD-E-protocol decoys | autodock-gpu (largest cluster) | 0.487 | 0.494 | +0.006 | 0.023 | 5 |
| DUD-E-protocol decoys | AutoDock-GPU per heavy atom | 0.504 | 0.493 | -0.011 | 0.019 | 5 |
| DUD-E-protocol decoys | interaction recall | 0.486 | 0.490 | +0.004 | 0.015 | 5 |

**Table 5. Docking against the best structure-free baseline**

| target | negatives | score | best baseline | baseline AUC | score AUC | difference | 95% CI | adjusted CI | verdict |
|---|---|---|---|---|---|---|---|---|---|
| ESR1_ant | DUD-E-protocol decoys | autodock-gpu | nn_sim | 0.745 | 0.427 | -0.319 | [-0.465, -0.162] | [-0.518, -0.097] | below |
| ESR1_ant | DUD-E-protocol decoys | consensus | nn_sim | 0.745 | 0.489 | -0.257 | [-0.389, -0.116] | [-0.435, -0.057] | below |
| ESR1_ant | DUD-E-protocol decoys | gnina | nn_sim | 0.745 | 0.542 | -0.203 | [-0.318, -0.085] | [-0.362, -0.036] | below |
| ESR1_ant | property-matched inactives | autodock-gpu | nn_sim | 0.603 | 0.487 | -0.117 | [-0.266, +0.044] | [-0.332, +0.100] | not resolved |
| ESR1_ant | property-matched inactives | consensus | nn_sim | 0.603 | 0.495 | -0.108 | [-0.239, +0.029] | [-0.288, +0.082] | not resolved |
| ESR1_ant | property-matched inactives | gnina | nn_sim | 0.603 | 0.491 | -0.113 | [-0.228, +0.003] | [-0.266, +0.046] | not resolved |
| ESR1_ant | random inactives | autodock-gpu | clogp | 0.741 | 0.68 | -0.061 | [-0.147, +0.028] | [-0.190, +0.062] | not resolved |
| ESR1_ant | random inactives | consensus | clogp | 0.741 | 0.681 | -0.06 | [-0.152, +0.034] | [-0.185, +0.066] | not resolved |
| ESR1_ant | random inactives | gnina | clogp | 0.741 | 0.606 | -0.135 | [-0.258, -0.010] | [-0.307, +0.038] | not resolved |
| MAPK1 | DUD-E-protocol decoys | autodock-gpu | nn_sim | 0.621 | 0.505 | -0.116 | [-0.219, -0.015] | [-0.254, +0.019] | not resolved |
| MAPK1 | DUD-E-protocol decoys | consensus | nn_sim | 0.621 | 0.558 | -0.063 | [-0.193, +0.063] | [-0.245, +0.106] | not resolved |
| MAPK1 | DUD-E-protocol decoys | gnina | nn_sim | 0.621 | 0.568 | -0.054 | [-0.213, +0.098] | [-0.267, +0.151] | not resolved |
| MAPK1 | property-matched inactives | autodock-gpu | clogp | 0.516 | 0.506 | -0.01 | [-0.125, +0.101] | [-0.164, +0.148] | not resolved |
| MAPK1 | property-matched inactives | consensus | clogp | 0.516 | 0.556 | 0.04 | [-0.080, +0.154] | [-0.118, +0.207] | not resolved |
| MAPK1 | property-matched inactives | gnina | clogp | 0.516 | 0.558 | 0.042 | [-0.091, +0.166] | [-0.141, +0.208] | not resolved |
| MAPK1 | random inactives | autodock-gpu | clogp | 0.697 | 0.565 | -0.132 | [-0.235, -0.029] | [-0.278, +0.017] | not resolved |
| MAPK1 | random inactives | consensus | clogp | 0.697 | 0.591 | -0.106 | [-0.218, +0.007] | [-0.270, +0.051] | not resolved |
| MAPK1 | random inactives | gnina | clogp | 0.697 | 0.58 | -0.116 | [-0.239, +0.002] | [-0.291, +0.052] | not resolved |
| MTORC1 | DUD-E-protocol decoys | autodock-gpu | nn_sim | 0.811 | 0.586 | -0.225 | [-0.330, -0.119] | [-0.370, -0.078] | below |
| MTORC1 | DUD-E-protocol decoys | consensus | nn_sim | 0.811 | 0.587 | -0.224 | [-0.325, -0.122] | [-0.369, -0.078] | below |
| MTORC1 | DUD-E-protocol decoys | gnina | nn_sim | 0.811 | 0.534 | -0.277 | [-0.394, -0.146] | [-0.442, -0.099] | below |
| MTORC1 | property-matched inactives | autodock-gpu | nn_sim | 0.633 | 0.507 | -0.126 | [-0.242, -0.011] | [-0.290, +0.042] | not resolved |
| MTORC1 | property-matched inactives | consensus | nn_sim | 0.633 | 0.474 | -0.159 | [-0.268, -0.046] | [-0.308, -0.007] | below |
| MTORC1 | property-matched inactives | gnina | nn_sim | 0.633 | 0.47 | -0.163 | [-0.289, -0.025] | [-0.337, +0.025] | not resolved |
| MTORC1 | random inactives | autodock-gpu | nn_sim | 0.685 | 0.533 | -0.151 | [-0.271, -0.024] | [-0.322, +0.022] | not resolved |
| MTORC1 | random inactives | consensus | nn_sim | 0.685 | 0.507 | -0.178 | [-0.295, -0.059] | [-0.340, -0.011] | below |
| MTORC1 | random inactives | gnina | nn_sim | 0.685 | 0.473 | -0.212 | [-0.352, -0.069] | [-0.397, -0.013] | below |
| PPARG | DUD-E-protocol decoys | autodock-gpu | nn_sim | 0.698 | 0.482 | -0.216 | [-0.340, -0.086] | [-0.390, -0.046] | below |
| PPARG | DUD-E-protocol decoys | consensus | nn_sim | 0.698 | 0.516 | -0.182 | [-0.323, -0.031] | [-0.391, +0.022] | not resolved |
| PPARG | DUD-E-protocol decoys | gnina | nn_sim | 0.698 | 0.546 | -0.151 | [-0.318, +0.024] | [-0.401, +0.087] | not resolved |
| PPARG | property-matched inactives | autodock-gpu | nn_sim | 0.656 | 0.581 | -0.075 | [-0.194, +0.048] | [-0.236, +0.087] | not resolved |
| PPARG | property-matched inactives | consensus | nn_sim | 0.656 | 0.599 | -0.057 | [-0.195, +0.083] | [-0.242, +0.133] | not resolved |
| PPARG | property-matched inactives | gnina | nn_sim | 0.656 | 0.624 | -0.032 | [-0.195, +0.139] | [-0.258, +0.196] | not resolved |
| PPARG | random inactives | autodock-gpu | n_heavy | 0.767 | 0.773 | 0.005 | [-0.055, +0.062] | [-0.085, +0.083] | not resolved |
| PPARG | random inactives | consensus | n_heavy | 0.767 | 0.809 | 0.042 | [-0.032, +0.109] | [-0.067, +0.133] | not resolved |
| PPARG | random inactives | gnina | n_heavy | 0.767 | 0.783 | 0.016 | [-0.094, +0.118] | [-0.143, +0.152] | not resolved |
| TP53 | DUD-E-protocol decoys | autodock-gpu | nn_sim | 0.665 | 0.466 | -0.199 | [-0.337, -0.056] | [-0.395, +0.008] | not resolved |
| TP53 | DUD-E-protocol decoys | consensus | nn_sim | 0.665 | 0.497 | -0.168 | [-0.284, -0.054] | [-0.321, -0.007] | below |
| TP53 | DUD-E-protocol decoys | gnina | nn_sim | 0.665 | 0.509 | -0.156 | [-0.272, -0.036] | [-0.312, +0.016] | not resolved |
| TP53 | property-matched inactives | autodock-gpu | nn_sim | 0.563 | 0.562 | -0.001 | [-0.134, +0.136] | [-0.178, +0.187] | not resolved |
| TP53 | property-matched inactives | consensus | nn_sim | 0.563 | 0.581 | 0.018 | [-0.092, +0.128] | [-0.136, +0.176] | not resolved |
| TP53 | property-matched inactives | gnina | nn_sim | 0.563 | 0.555 | -0.008 | [-0.118, +0.102] | [-0.155, +0.143] | not resolved |
| TP53 | random inactives | autodock-gpu | mw | 0.603 | 0.597 | -0.007 | [-0.097, +0.085] | [-0.138, +0.121] | not resolved |
| TP53 | random inactives | consensus | mw | 0.603 | 0.621 | 0.018 | [-0.094, +0.126] | [-0.138, +0.160] | not resolved |
| TP53 | random inactives | gnina | mw | 0.603 | 0.56 | -0.043 | [-0.181, +0.092] | [-0.237, +0.142] | not resolved |

**Table 6. Resampling sensitivity: compound, matched-group and scaffold-aware bootstraps**

| target | score | auc | registered_low | registered_high | grouped_low | grouped_high | scaffold_low | scaffold_high | n_groups | n_clusters | registered_excludes_chance | registered_width | grouped_excludes_chance | grouped_width | scaffold_excludes_chance | scaffold_width |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ESR1_ant | autodock-gpu | 0.487 | 0.392 | 0.581 | 0.43 | 0.54 | 0.432 | 0.542 | 40 | 33 | False | 0.189 | False | 0.111 | False | 0.11 |
| ESR1_ant | gnina | 0.491 | 0.382 | 0.598 | 0.392 | 0.585 | 0.378 | 0.586 | 40 | 33 | False | 0.215 | False | 0.193 | False | 0.207 |
| ESR1_ant | consensus | 0.495 | 0.396 | 0.589 | 0.426 | 0.559 | 0.424 | 0.558 | 40 | 33 | False | 0.193 | False | 0.133 | False | 0.135 |
| MAPK1 | autodock-gpu | 0.506 | 0.405 | 0.608 | 0.441 | 0.567 | 0.435 | 0.575 | 40 | 37 | False | 0.203 | False | 0.126 | False | 0.14 |
| MAPK1 | gnina | 0.558 | 0.462 | 0.658 | 0.477 | 0.639 | 0.474 | 0.638 | 40 | 37 | False | 0.196 | False | 0.162 | False | 0.163 |
| MAPK1 | consensus | 0.556 | 0.454 | 0.657 | 0.475 | 0.63 | 0.47 | 0.642 | 40 | 37 | False | 0.203 | False | 0.155 | False | 0.173 |
| MTORC1 | autodock-gpu | 0.507 | 0.411 | 0.601 | 0.448 | 0.569 | 0.441 | 0.574 | 40 | 32 | False | 0.19 | False | 0.122 | False | 0.133 |
| MTORC1 | gnina | 0.47 | 0.381 | 0.558 | 0.391 | 0.549 | 0.391 | 0.542 | 40 | 32 | False | 0.177 | False | 0.158 | False | 0.152 |
| MTORC1 | consensus | 0.474 | 0.39 | 0.555 | 0.395 | 0.552 | 0.406 | 0.542 | 40 | 32 | False | 0.166 | False | 0.156 | False | 0.136 |
| PPARG | autodock-gpu | 0.581 | 0.466 | 0.684 | 0.495 | 0.672 | 0.481 | 0.672 | 27 | 21 | False | 0.218 | False | 0.177 | False | 0.191 |
| PPARG | gnina | 0.624 | 0.505 | 0.737 | 0.507 | 0.741 | 0.498 | 0.749 | 27 | 21 | True | 0.232 | True | 0.235 | False | 0.251 |
| PPARG | consensus | 0.599 | 0.485 | 0.711 | 0.496 | 0.705 | 0.481 | 0.712 | 27 | 21 | False | 0.226 | False | 0.209 | False | 0.231 |
| TP53 | autodock-gpu | 0.562 | 0.463 | 0.665 | 0.481 | 0.639 | 0.471 | 0.653 | 40 | 33 | False | 0.202 | False | 0.158 | False | 0.182 |
| TP53 | gnina | 0.555 | 0.464 | 0.643 | 0.473 | 0.643 | 0.461 | 0.648 | 40 | 33 | False | 0.179 | False | 0.17 | False | 0.187 |
| TP53 | consensus | 0.581 | 0.486 | 0.678 | 0.506 | 0.664 | 0.492 | 0.665 | 40 | 33 | False | 0.192 | True | 0.158 | False | 0.173 |

**Table 7. Retrospective and deployable ligand-similarity baselines**

| target | negative_set | baseline | n_reference_actives | roc_auc | ci_low | ci_high | excludes_chance | autodock-gpu_auc | autodock-gpu_minus_baseline | autodock-gpu_beats | autodock-gpu_below | gnina_auc | gnina_minus_baseline | gnina_beats | gnina_below | consensus_auc | consensus_minus_baseline | consensus_beats | consensus_below | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ESR1_ant | random | similarity (oracle, evaluation-set actives) | 39 | 0.685 | 0.582 | 0.774 | True | 0.68 | -0.0045 | False | False | 0.606 | -0.079 | False | False | 0.681 | -0.0035 | False | False |  |
| ESR1_ant | random | similarity (deployable, held-out actives) | 62 | 0.732 | 0.637 | 0.819 | True | 0.68 | -0.0519 | False | False | 0.606 | -0.1263 | False | False | 0.681 | -0.0508 | False | False |  |
| ESR1_ant | matched | similarity (oracle, evaluation-set actives) | 39 | 0.603 | 0.491 | 0.704 | False | 0.487 | -0.1166 | False | False | 0.491 | -0.1128 | False | False | 0.495 | -0.1081 | False | False |  |
| ESR1_ant | matched | similarity (deployable, held-out actives) | 62 | 0.686 | 0.578 | 0.78 | True | 0.487 | -0.1996 | False | True | 0.491 | -0.1959 | False | True | 0.495 | -0.1912 | False | True |  |
| ESR1_ant | decoy | similarity (oracle, evaluation-set actives) | 39 | 0.745 | 0.638 | 0.839 | True | 0.427 | -0.3187 | False | True | 0.542 | -0.2032 | False | True | 0.489 | -0.2567 | False | True |  |
| ESR1_ant | decoy | similarity (deployable, held-out actives) | 62 | 0.739 | 0.633 | 0.832 | True | 0.427 | -0.3125 | False | True | 0.542 | -0.1971 | False | True | 0.489 | -0.2505 | False | True |  |
| MAPK1 | random | similarity (oracle, evaluation-set actives) | 39 | 0.568 | 0.453 | 0.676 | False | 0.565 | -0.0034 | False | False | 0.58 | 0.0119 | False | False | 0.591 | 0.0227 | False | False |  |
| MAPK1 | random | similarity (deployable, held-out actives) | 268 | 0.71 | 0.616 | 0.803 | True | 0.565 | -0.1457 | False | True | 0.58 | -0.1303 | False | False | 0.591 | -0.1195 | False | False |  |
| MAPK1 | matched | similarity (oracle, evaluation-set actives) | 39 | 0.499 | 0.383 | 0.615 | False | 0.506 | 0.007 | False | False | 0.558 | 0.0587 | False | False | 0.556 | 0.0564 | False | False |  |
| MAPK1 | matched | similarity (deployable, held-out actives) | 268 | 0.665 | 0.565 | 0.759 | True | 0.506 | -0.159 | False | True | 0.558 | -0.1073 | False | False | 0.556 | -0.1096 | False | False |  |
| MAPK1 | decoy | similarity (oracle, evaluation-set actives) | 39 | 0.621 | 0.502 | 0.731 | True | 0.505 | -0.1164 | False | False | 0.568 | -0.0537 | False | False | 0.558 | -0.0634 | False | False |  |
| MAPK1 | decoy | similarity (deployable, held-out actives) | 268 | 0.724 | 0.629 | 0.813 | True | 0.505 | -0.2187 | False | True | 0.568 | -0.1561 | False | False | 0.558 | -0.1658 | False | True |  |
| MTORC1 | random | similarity (oracle, evaluation-set actives) | 39 | 0.685 | 0.562 | 0.794 | True | 0.533 | -0.1513 | False | False | 0.473 | -0.2117 | False | True | 0.507 | -0.178 | False | True |  |
| MTORC1 | random | similarity (deployable, held-out actives) | 57 | 0.708 | 0.601 | 0.812 | True | 0.533 | -0.1747 | False | True | 0.473 | -0.235 | False | True | 0.507 | -0.2014 | False | True |  |
| MTORC1 | matched | similarity (oracle, evaluation-set actives) | 39 | 0.633 | 0.527 | 0.738 | True | 0.507 | -0.1261 | False | False | 0.47 | -0.1626 | False | False | 0.474 | -0.1593 | False | True |  |
| MTORC1 | matched | similarity (deployable, held-out actives) | 57 | 0.705 | 0.597 | 0.8 | True | 0.507 | -0.1986 | False | True | 0.47 | -0.235 | False | True | 0.474 | -0.2318 | False | True |  |
| MTORC1 | decoy | similarity (oracle, evaluation-set actives) | 39 | 0.811 | 0.714 | 0.896 | True | 0.586 | -0.2246 | False | True | 0.534 | -0.2769 | False | True | 0.587 | -0.2237 | False | True |  |
| MTORC1 | decoy | similarity (deployable, held-out actives) | 57 | 0.813 | 0.72 | 0.894 | True | 0.586 | -0.2273 | False | True | 0.534 | -0.2796 | False | True | 0.587 | -0.2264 | False | True |  |
| PPARG | random | similarity (oracle, evaluation-set actives) | 26 | 0.75 | 0.624 | 0.859 | True | 0.773 | 0.0231 | False | False | 0.783 | 0.0333 | False | False | 0.809 | 0.0594 | False | False |  |
| PPARG | random | similarity (deployable, held-out actives) | 0 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | no held-out actives |
| PPARG | matched | similarity (oracle, evaluation-set actives) | 26 | 0.656 | 0.519 | 0.785 | True | 0.581 | -0.0748 | False | False | 0.624 | -0.0318 | False | False | 0.599 | -0.0568 | False | False |  |
| PPARG | matched | similarity (deployable, held-out actives) | 0 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | no held-out actives |
| PPARG | decoy | similarity (oracle, evaluation-set actives) | 26 | 0.698 | 0.555 | 0.828 | True | 0.482 | -0.2161 | False | True | 0.546 | -0.1514 | False | False | 0.516 | -0.1818 | False | True |  |
| PPARG | decoy | similarity (deployable, held-out actives) | 0 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | no held-out actives |
| TP53 | random | similarity (oracle, evaluation-set actives) | 39 | 0.578 | 0.449 | 0.688 | False | 0.597 | 0.0191 | False | False | 0.56 | -0.0174 | False | False | 0.621 | 0.0435 | False | False |  |
| TP53 | random | similarity (deployable, held-out actives) | 39 | 0.659 | 0.542 | 0.764 | True | 0.597 | -0.0622 | False | False | 0.56 | -0.0987 | False | False | 0.621 | -0.0378 | False | False |  |
| TP53 | matched | similarity (oracle, evaluation-set actives) | 39 | 0.563 | 0.432 | 0.67 | False | 0.562 | -0.0006 | False | False | 0.555 | -0.0082 | False | False | 0.581 | 0.0182 | False | False |  |
| TP53 | matched | similarity (deployable, held-out actives) | 39 | 0.635 | 0.513 | 0.741 | True | 0.562 | -0.0724 | False | False | 0.555 | -0.08 | False | False | 0.581 | -0.0536 | False | False |  |
| TP53 | decoy | similarity (oracle, evaluation-set actives) | 39 | 0.665 | 0.522 | 0.775 | True | 0.466 | -0.199 | False | True | 0.509 | -0.1557 | False | True | 0.497 | -0.1683 | False | True |  |
| TP53 | decoy | similarity (deployable, held-out actives) | 39 | 0.688 | 0.562 | 0.791 | True | 0.466 | -0.2216 | False | True | 0.509 | -0.1783 | False | True | 0.497 | -0.191 | False | True |  |

**Table 8. Docking-box strategies**

| strategy | engine | n | successes | success_rate | ci_low | ci_high | median_rmsd | median_wall_s |
|---|---|---|---|---|---|---|---|---|
| ligand | autodock-gpu | 62 | 18 | 0.29 | 0.192 | 0.413 | 4.969 | 24.88 |
| ligand | gnina | 62 | 26 | 0.419 | 0.305 | 0.543 | 1.42 | 24.88 |
| pocket | autodock-gpu | 62 | 10 | 0.161 | 0.09 | 0.272 | 18.529 | 24.76 |
| pocket | gnina | 62 | 14 | 0.226 | 0.14 | 0.344 | 19.54 | 24.76 |
| p2rank | autodock-gpu | 62 | 17 | 0.274 | 0.179 | 0.396 | 4.92 | 19.736 |
| p2rank | gnina | 62 | 23 | 0.371 | 0.262 | 0.495 | 1.855 | 19.736 |
| blind | autodock-gpu | 62 | 13 | 0.21 | 0.127 | 0.326 | 5.483 | 55.59 |
| blind | gnina | 62 | 20 | 0.323 | 0.22 | 0.446 | 1.994 | 55.59 |

**Table 9. Docking-box strategies, paired per template (exact McNemar)**

| engine | strategy_a | strategy_b | scope | n_paired | a_only | b_only | discordant | p_mcnemar_exact | favours |
|---|---|---|---|---|---|---|---|---|---|
| autodock-gpu | ligand | pocket | pairwise | 43 | 8 | 0 | 8 | 0.0078 | ligand |
| autodock-gpu | ligand | p2rank | pairwise | 44 | 4 | 3 | 7 | 1.0 | not resolved |
| autodock-gpu | ligand | blind | pairwise | 39 | 3 | 0 | 3 | 0.25 | not resolved |
| autodock-gpu | pocket | p2rank | pairwise | 43 | 2 | 9 | 11 | 0.0654 | not resolved |
| autodock-gpu | pocket | blind | pairwise | 39 | 3 | 7 | 10 | 0.3438 | not resolved |
| autodock-gpu | p2rank | blind | pairwise | 39 | 5 | 1 | 6 | 0.2188 | not resolved |
| gnina | ligand | pocket | pairwise | 43 | 12 | 1 | 13 | 0.0034 | ligand |
| gnina | ligand | p2rank | pairwise | 44 | 3 | 0 | 3 | 0.25 | not resolved |
| gnina | ligand | blind | pairwise | 39 | 3 | 0 | 3 | 0.25 | not resolved |
| gnina | pocket | p2rank | pairwise | 43 | 3 | 12 | 15 | 0.0352 | p2rank |
| gnina | pocket | blind | pairwise | 39 | 2 | 10 | 12 | 0.0386 | blind |
| gnina | p2rank | blind | pairwise | 39 | 3 | 0 | 3 | 0.25 | not resolved |
| autodock-gpu | ligand | pocket | listwise | 39 | 7 | 0 | 7 | 0.0156 | ligand |
| autodock-gpu | ligand | p2rank | listwise | 39 | 2 | 3 | 5 | 1.0 | not resolved |
| autodock-gpu | ligand | blind | listwise | 39 | 3 | 0 | 3 | 0.25 | not resolved |
| autodock-gpu | pocket | p2rank | listwise | 39 | 1 | 9 | 10 | 0.0215 | p2rank |
| autodock-gpu | pocket | blind | listwise | 39 | 3 | 7 | 10 | 0.3438 | not resolved |
| autodock-gpu | p2rank | blind | listwise | 39 | 5 | 1 | 6 | 0.2188 | not resolved |
| gnina | ligand | pocket | listwise | 39 | 12 | 1 | 13 | 0.0034 | ligand |
| gnina | ligand | p2rank | listwise | 39 | 0 | 0 | 0 | 1.0 | not resolved |
| gnina | ligand | blind | listwise | 39 | 3 | 0 | 3 | 0.25 | not resolved |
| gnina | pocket | p2rank | listwise | 39 | 1 | 12 | 13 | 0.0034 | p2rank |
| gnina | pocket | blind | listwise | 39 | 2 | 10 | 12 | 0.0386 | blind |
| gnina | p2rank | blind | listwise | 39 | 3 | 0 | 3 | 0.25 | not resolved |
