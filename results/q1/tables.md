## Main-text tables

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

## Additional file: supplementary tables

**Table S1. Pose selection**

| target | negative set | AUC lowest energy | AUC largest cluster | difference [95% CI] |
|---|---|---|---|---|
| ESR1_ant | random | 0.68 | 0.68 | -0.002 [-0.018, +0.011] |
| ESR1_ant | matched | 0.49 | 0.48 | -0.005 [-0.017, +0.006] |
| ESR1_ant | decoy | 0.43 | 0.43 | +0.001 [-0.014, +0.015] |
| MAPK1 | random | 0.56 | 0.56 | -0.000 [-0.030, +0.026] |
| MAPK1 | matched | 0.51 | 0.51 | -0.000 [-0.029, +0.024] |
| MAPK1 | decoy | 0.50 | 0.49 | -0.013 [-0.042, +0.013] |
| MTORC1 | random | 0.53 | 0.51 | -0.025 [-0.088, +0.030] |
| MTORC1 | matched | 0.51 | 0.47 | -0.041 [-0.100, +0.014] |
| MTORC1 | decoy | 0.59 | 0.56 | -0.026 [-0.082, +0.032] |
| PPARG | random | 0.77 | 0.78 | +0.003 [-0.011, +0.015] |
| PPARG | matched | 0.58 | 0.58 | -0.001 [-0.030, +0.021] |
| PPARG | decoy | 0.48 | 0.48 | +0.000 [-0.028, +0.023] |
| TP53 | random | 0.60 | 0.61 | +0.010 [-0.017, +0.034] |
| TP53 | matched | 0.56 | 0.57 | +0.011 [-0.010, +0.030] |
| TP53 | decoy | 0.47 | 0.47 | +0.008 [-0.019, +0.031] |

**Table S2. Throughput**

| target | compounds | ligand prep (s/compound) | AutoDock-GPU (s/compound) | gnina (s/compound) |
|---|---|---|---|---|
| ESR1_ant | 980 | 0.004 | 0.80 | 3.91 |
| MAPK1 | 980 | 0.018 | 1.35 | 5.57 |
| MTORC1 | 980 | 0.010 | 2.01 | 6.87 |
| PPARG | 759 | 0.005 | 1.11 | 4.03 |
| TP53 | 980 | 0.006 | 1.42 | 3.37 |

**Table S3. Throughput, CPU against GPU**

| case | n | threads | wall_mean_s | wall_sd_s | n_ok | peak_vram_mib | mean_gpu_util | repeats | s_per_compound | compounds_per_hour | speedup_vs_ref |
|---|---|---|---|---|---|---|---|---|---|---|---|
| autodock-gpu | 8 | 16 | 5.036 | 0.061 | 8 | 295.0 | 80.2333 | 3 | 0.629 | 5719.113 |  |
| autodock-gpu | 16 | 16 | 11.968 | 0.042 | 16 | 295.0 | 89.1364 | 3 | 0.748 | 4812.779 | 8.998 |
| autodock-gpu | 32 | 16 | 36.36 | 0.124 | 32 | 295.0 | 93.6776 | 3 | 1.136 | 3168.347 |  |
| gnina (CNN, GPU) | 8 | 16 | 30.205 | 0.207 | 8 | 175.0 | 15.2711 | 3 | 3.776 | 953.473 |  |
| gnina (CNN, GPU) | 16 | 16 | 64.05 | 0.217 | 16 | 175.0 | 14.4526 | 3 | 4.003 | 899.293 | 1.681 |
| gnina (CNN, GPU) | 32 | 16 | 127.89 | 0.385 | 32 | 175.0 | 15.4584 | 3 | 3.997 | 900.776 |  |
| gnina (CNN, CPU) | 16 | 1 | 656.785 | 42.926 | 16 | 81.0 | 0.0 | 3 | 41.049 | 87.7 | 0.164 |
| gnina (CNN, CPU) | 16 | 4 | 481.694 | 24.455 | 16 | 81.0 | 0.0 | 3 | 30.106 | 119.578 | 0.224 |
| gnina (CNN, CPU) | 16 | 8 | 660.219 | 6.884 | 16 | 81.0 | 0.0 | 3 | 41.264 | 87.244 | 0.163 |
| gnina (CNN, CPU) | 16 | 16 | 1048.799 | 5.786 | 16 | 81.0 | 0.0 | 3 | 65.55 | 54.92 | 0.103 |
| gnina (Vina, CPU) | 16 | 1 | 247.41 | 0.416 | 16 | 81.0 | 0.0008 | 3 | 15.463 | 232.812 | 0.435 |
| gnina (Vina, CPU) | 16 | 4 | 76.652 | 0.479 | 16 | 81.0 | 0.0097 | 3 | 4.791 | 751.451 | 1.405 |
| gnina (Vina, CPU) | 16 | 8 | 51.483 | 0.06 | 16 | 81.0 | 0.0074 | 3 | 3.218 | 1118.817 | 2.092 |
| gnina (Vina, CPU) | 16 | 16 | 51.733 | 0.47 | 16 | 81.0 | 0.0037 | 3 | 3.233 | 1113.409 | 2.082 |
| vina (CPU) | 16 | 1 | 305.157 | 1.151 | 16 | 0.0 | 0.0 | 3 | 19.072 | 188.755 | 0.353 |
| vina (CPU) | 16 | 4 | 126.906 | 0.424 | 16 | 0.0 | 0.0 | 3 | 7.932 | 453.879 | 0.849 |
| vina (CPU) | 16 | 8 | 113.956 | 0.454 | 16 | 0.0 | 0.0 | 3 | 7.122 | 505.458 | 0.945 |
| vina (CPU) | 16 | 16 | 107.684 | 0.362 | 16 | 0.0 | 0.0 | 3 | 6.73 | 534.9 | 1.0 |

**Table S4. Failures**

| target | class | n | prep failed | AutoDock-GPU failed | gnina failed | gnina positive energy | chelator needs manual ionisation | chelator deprotonated |
|---|---|---|---|---|---|---|---|---|
| ESR1_ant | active | 40 | 0 | 0 | 0 | 5% | 0 | 0 |
| ESR1_ant | inactive | 940 | 0 | 0 | 0 | 4% | 13 | 0 |
| MAPK1 | active | 40 | 0 | 0 | 0 | 0% | 1 | 0 |
| MAPK1 | inactive | 940 | 2 | 0 | 0 | 0% | 5 | 0 |
| MTORC1 | active | 40 | 0 | 0 | 0 | 0% | 0 | 0 |
| MTORC1 | inactive | 940 | 0 | 0 | 0 | 0% | 0 | 0 |
| PPARG | active | 27 | 0 | 0 | 0 | 0% | 0 | 0 |
| PPARG | inactive | 732 | 0 | 0 | 0 | 1% | 13 | 0 |
| TP53 | active | 40 | 0 | 0 | 0 | 0% | 0 | 0 |
| TP53 | inactive | 940 | 1 | 0 | 0 | 0% | 26 | 0 |

**Table S5. Metal-ion audit of the templates against their original depositions**

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

**Table S6. Change in ROC-AUC between negative-set designs**

| negatives A | negatives B | score | mean_delta_auc | targets | targets_with_ci_excluding_zero |
|---|---|---|---|---|---|
| random | matched | autodock-gpu | +0.101 | 5 | 3 |
| random | matched | gnina | +0.061 | 5 | 2 |
| random | matched | consensus | +0.101 | 5 | 2 |
| matched | decoy | autodock-gpu | +0.036 | 5 | 4 |
| matched | decoy | gnina | -0.000 | 5 | 3 |
| matched | decoy | consensus | +0.012 | 5 | 3 |
| random | decoy | autodock-gpu | +0.137 | 5 | 5 |
| random | decoy | gnina | +0.061 | 5 | 4 |
| random | decoy | consensus | +0.113 | 5 | 4 |

**Table S7. Seed stability of AutoDock-GPU**

| target | n_seeds | auc_min | auc_max | auc_range | mean_spearman_between_seeds | stable_top_n_fraction | mean_per_compound_sd |
|---|---|---|---|---|---|---|---|
| ESR1_ant | 3 | 0.673 | 0.684 | 0.011 | 0.994 | 0.950 | 0.045 |
| MAPK1 | 3 | 0.555 | 0.570 | 0.015 | 0.984 | 0.650 | 0.090 |
| MTORC1 | 3 | 0.527 | 0.536 | 0.009 | 0.936 | 0.600 | 0.157 |
| PPARG | 3 | 0.771 | 0.775 | 0.004 | 0.997 | 0.950 | 0.049 |
| TP53 | 3 | 0.593 | 0.598 | 0.005 | 0.983 | 0.750 | 0.073 |

**Table S8. Similarity of the negative sets to the actives**

| target | actives | median NN similarity among actives | clusters at 0.4 | analogue AUC vs random | analogue AUC vs matched | analogue AUC vs decoy |
|---|---|---|---|---|---|---|
| ESR1_ant | 40 | 0.24 | 33 | 0.68 | 0.60 | 0.75 |
| MAPK1 | 40 | 0.24 | 37 | 0.57 | 0.50 | 0.62 |
| MTORC1 | 40 | 0.37 | 32 | 0.68 | 0.63 | 0.81 |
| PPARG | 27 | 0.29 | 21 | 0.75 | 0.66 | 0.70 |
| TP53 | 40 | 0.22 | 33 | 0.58 | 0.56 | 0.67 |

**Table S9. Median heavy-atom gap, actives against each negative set**

| target | decoy | matched | random |
|---|---|---|---|
| ESR1_ant | -1.0 | +0.5 | +5.5 |
| MAPK1 | +0.0 | +0.0 | +1.0 |
| MTORC1 | +0.0 | +0.0 | +0.0 |
| PPARG | +0.0 | +2.0 | +7.0 |
| TP53 | -1.0 | +1.0 | +2.0 |

**Table S10. ADMET risk tiers**

| target | actives PASS | actives MODERATE | actives HIGH_RISK | actives EXCLUDE | negatives PASS | negatives MODERATE | negatives HIGH_RISK | negatives EXCLUDE |
|---|---|---|---|---|---|---|---|---|
| ESR1_ant | 45% | 2% | 35% | 18% | 57% | 2% | 22% | 19% |
| MAPK1 | 28% | 0% | 32% | 40% | 38% | 1% | 40% | 21% |
| MTORC1 | 28% | 0% | 52% | 20% | 22% | 0% | 63% | 14% |
| PPARG | 67% | 11% | 7% | 15% | 68% | 3% | 13% | 16% |
| TP53 | 25% | 2% | 10% | 62% | 61% | 1% | 13% | 25% |

**Table S11. Effect of ADMET filters on retrieval**

| negative_set | filter | actives_retained | inactives_retained | actives_lost | auc_before | auc_after | targets |
|---|---|---|---|---|---|---|---|
| random | not EXCLUDE (no PAINS, Ames <= 0.5) | 69% | 82% | 60 | 0.642 | 0.653 | 5 |
| random | not HIGH_RISK/EXCLUDE (also hERG <= 0.5) | 42% | 53% | 114 | 0.642 | 0.651 | 5 |
| random | Veber pass | 96% | 96% | 6 | 0.642 | 0.645 | 5 |
| random | no PAINS alert | 91% | 96% | 17 | 0.642 | 0.647 | 5 |
| matched | not EXCLUDE (no PAINS, Ames <= 0.5) | 69% | 83% | 60 | 0.541 | 0.552 | 5 |
| matched | not HIGH_RISK/EXCLUDE (also hERG <= 0.5) | 42% | 53% | 114 | 0.541 | 0.532 | 5 |
| matched | Veber pass | 96% | 98% | 6 | 0.541 | 0.542 | 5 |
| matched | no PAINS alert | 91% | 95% | 17 | 0.541 | 0.546 | 5 |
| decoy | not EXCLUDE (no PAINS, Ames <= 0.5) | 69% | 79% | 60 | 0.529 | 0.557 | 5 |
| decoy | not HIGH_RISK/EXCLUDE (also hERG <= 0.5) | 42% | 47% | 114 | 0.529 | 0.564 | 5 |
| decoy | Veber pass | 96% | 99% | 6 | 0.529 | 0.528 | 5 |
| decoy | no PAINS alert | 91% | 94% | 17 | 0.529 | 0.538 | 5 |

**Table S12. Hit rate of the lead-selection strategies**

| negative_set | strategy | top_n | hit_rate | sd_across_targets | pool_active_fraction | targets |
|---|---|---|---|---|---|---|
| random | docking | 10 | 0.220 | 0.228 | 0.182 | 5 |
| random | docking | 20 | 0.250 | 0.166 | 0.182 | 5 |
| random | docking | 50 | 0.220 | 0.094 | 0.182 | 5 |
| random | docking+ADMET | 10 | 0.220 | 0.179 | 0.182 | 5 |
| random | docking+ADMET | 20 | 0.211 | 0.154 | 0.182 | 5 |
| random | docking+ADMET | 50 | 0.148 | 0.097 | 0.182 | 5 |
| random | docking+IFP | 10 | 0.325 | 0.171 | 0.203 | 4 |
| random | docking+IFP | 20 | 0.263 | 0.075 | 0.203 | 4 |
| random | docking+IFP | 50 | 0.240 | 0.073 | 0.203 | 4 |
| random | pareto | 10 | 0.300 | 0.163 | 0.203 | 4 |
| random | pareto | 20 | 0.275 | 0.132 | 0.203 | 4 |
| random | pareto | 50 | 0.245 | 0.077 | 0.203 | 4 |
| matched | docking | 10 | 0.140 | 0.089 | 0.136 | 5 |
| matched | docking | 20 | 0.130 | 0.045 | 0.136 | 5 |
| matched | docking | 50 | 0.132 | 0.046 | 0.136 | 5 |
| matched | docking+ADMET | 10 | 0.080 | 0.084 | 0.136 | 5 |
| matched | docking+ADMET | 20 | 0.100 | 0.061 | 0.136 | 5 |
| matched | docking+ADMET | 50 | 0.100 | 0.051 | 0.136 | 5 |
| matched | docking+IFP | 10 | 0.125 | 0.050 | 0.155 | 4 |
| matched | docking+IFP | 20 | 0.213 | 0.025 | 0.155 | 4 |
| matched | docking+IFP | 50 | 0.160 | 0.037 | 0.155 | 4 |
| matched | pareto | 10 | 0.150 | 0.058 | 0.155 | 4 |
| matched | pareto | 20 | 0.175 | 0.050 | 0.155 | 4 |
| matched | pareto | 50 | 0.140 | 0.028 | 0.155 | 4 |
| decoy | docking | 10 | 0.140 | 0.089 | 0.136 | 5 |
| decoy | docking | 20 | 0.110 | 0.074 | 0.136 | 5 |
| decoy | docking | 50 | 0.132 | 0.059 | 0.136 | 5 |
| decoy | docking+ADMET | 10 | 0.160 | 0.167 | 0.136 | 5 |
| decoy | docking+ADMET | 20 | 0.143 | 0.092 | 0.136 | 5 |
| decoy | docking+ADMET | 50 | 0.147 | 0.061 | 0.136 | 5 |
| decoy | docking+IFP | 10 | 0.175 | 0.150 | 0.138 | 4 |
| decoy | docking+IFP | 20 | 0.163 | 0.075 | 0.138 | 4 |
| decoy | docking+IFP | 50 | 0.135 | 0.041 | 0.138 | 4 |
| decoy | pareto | 10 | 0.175 | 0.096 | 0.138 | 4 |
| decoy | pareto | 20 | 0.175 | 0.087 | 0.138 | 4 |
| decoy | pareto | 50 | 0.125 | 0.034 | 0.138 | 4 |

**Table S13. Recovery of the co-crystal interactions by the top gnina pose**

| target | crystal contacts | recall, actives | n, actives | recall, top-ranked negatives | n, top-ranked negatives |
|---|---|---|---|---|---|
| ESR1_ant | 16 | 0.21 | 40 | 0.24 | 940 |
| MAPK1 | 8 | 0.40 | 40 | 0.40 | 938 |
| MTORC1 | 7 | 0.43 | 40 | 0.47 | 940 |
| PPARG | 9 | 0.44 | 27 | 0.36 | 732 |
| TP53 | 5 | 0.23 | 40 | 0.23 | 939 |

**Table S14. Structure quality of the templates**

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

**Table S15. Pocket-detector recovery of the ligand site**

| detector | templates | top-1 recovery | top-3 recovery | mean best-pocket recovery | median pockets |
|---|---|---|---|---|---|
| fpocket | 129 | 37/129 (29%) [22%, 37%] | 61/129 (47%) [39%, 56%] | 0.89 | 22 |
| dogsite | 63 | 27/63 (43%) [31%, 55%] | 55/63 (87%) [77%, 93%] | 0.92 | 13 |
| p2rank | 129 | 99/129 (77%) [69%, 83%] | 123/129 (95%) [90%, 98%] | 0.87 | 8 |

**Table S16. Druggability controls**

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

**Table S17. Quality-gate ablation**

| condition | n | mean_dG | median_dG | min_dG | max_dG | non_binding_dG_ge_0 | dG_below_minus20 | n_charge_flagged | n_torsion_frozen_flagged | mean_shift_vs_baseline | shift_ci_low | shift_ci_high | target | template | n_ligands |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| baseline | 60 | -7.403 | -7.65 | -10.5 | -2.77 | 0 | 0 | 0 | 0 | 0.0 | 0.0 | 0.0 | ESR1_ant | 6b0f | 60 |
| charges_x4_uncapped | 60 | -4.01 | -3.66 | -8.51 | 0.03 | 1 | 0 | 60 | 0 | 3.392 | 2.818 | 3.935 | ESR1_ant | 6b0f | 60 |
| charges_x4_capped | 60 | -4.105 | -3.855 | -8.69 | -0.92 | 0 | 0 | 0 | 0 | 3.298 | 2.854 | 3.757 | ESR1_ant | 6b0f | 60 |
| torsions_frozen | 60 | -7.681 | -8.01 | -10.68 | -4.26 | 0 | 0 | 0 | 32 | -0.278 | -0.531 | -0.012 | ESR1_ant | 6b0f | 60 |
| baseline | 60 | -7.975 | -8.155 | -10.54 | -5.02 | 0 | 0 | 0 | 0 | 0.0 | 0.0 | 0.0 | MAPK1 | 4qte | 60 |
| charges_x4_uncapped | 60 | -3.474 | -3.41 | -7.91 | -0.05 | 0 | 0 | 60 | 0 | 4.501 | 4.166 | 4.83 | MAPK1 | 4qte | 60 |
| charges_x4_capped | 60 | -3.733 | -3.53 | -6.76 | -0.58 | 0 | 0 | 0 | 0 | 4.242 | 3.916 | 4.574 | MAPK1 | 4qte | 60 |
| torsions_frozen | 60 | -6.978 | -7.085 | -9.02 | -5.04 | 0 | 0 | 0 | 51 | 0.997 | 0.794 | 1.201 | MAPK1 | 4qte | 60 |
| baseline | 60 | -9.156 | -9.205 | -11.1 | -6.3 | 0 | 0 | 0 | 0 | 0.0 | 0.0 | 0.0 | MTORC1 | 4jsx | 60 |
| charges_x4_uncapped | 60 | -4.083 | -3.495 | -9.37 | 0.32 | 1 | 0 | 60 | 0 | 5.073 | 4.55 | 5.569 | MTORC1 | 4jsx | 60 |
| charges_x4_capped | 60 | -5.122 | -4.66 | -10.21 | -0.33 | 0 | 0 | 0 | 0 | 4.034 | 3.547 | 4.516 | MTORC1 | 4jsx | 60 |
| torsions_frozen | 60 | -8.199 | -8.13 | -9.81 | -6.12 | 0 | 0 | 0 | 56 | 0.957 | 0.711 | 1.175 | MTORC1 | 4jsx | 60 |
| baseline | 57 | -7.777 | -7.98 | -10.66 | -4.59 | 0 | 0 | 0 | 0 | 0.0 | 0.0 | 0.0 | PPARG | 2p4y | 57 |
| charges_x4_uncapped | 57 | -4.031 | -4.87 | -8.22 | 0.77 | 2 | 0 | 57 | 0 | 3.746 | 3.146 | 4.35 | PPARG | 2p4y | 57 |
| charges_x4_capped | 57 | -4.086 | -4.25 | -7.54 | 0.39 | 1 | 0 | 0 | 0 | 3.691 | 3.238 | 4.158 | PPARG | 2p4y | 57 |
| torsions_frozen | 57 | -7.91 | -8.06 | -10.02 | -5.51 | 0 | 0 | 0 | 51 | -0.133 | -0.317 | 0.068 | PPARG | 2p4y | 57 |
| baseline | 60 | -6.555 | -6.735 | -8.71 | -3.53 | 0 | 0 | 0 | 0 | 0.0 | 0.0 | 0.0 | TP53 | 4agq | 60 |
| charges_x4_uncapped | 60 | -2.641 | -2.29 | -6.7 | 0.21 | 2 | 0 | 55 | 0 | 3.915 | 3.403 | 4.401 | TP53 | 4agq | 60 |
| charges_x4_capped | 60 | -3.519 | -3.455 | -6.84 | -0.92 | 0 | 0 | 0 | 0 | 3.037 | 2.653 | 3.419 | TP53 | 4agq | 60 |
| torsions_frozen | 60 | -6.365 | -6.365 | -7.93 | -4.01 | 0 | 0 | 0 | 33 | 0.19 | 0.031 | 0.348 | TP53 | 4agq | 60 |

**Table S18. Template difficulty by target**

| target | templates | redocked | validated (<=2 A) | macrocyclic ligands | median heavy atoms | median MW | best RMSD (A) |
|---|---|---|---|---|---|---|---|
| ESR1_ant | 15 | 13 | 10 | 0 | 33 | 428 | 0.48 |
| MAPK1 | 15 | 7 | 4 | 1 | 26 | 377 | 0.60 |
| MTORC1 | 11 | 10 | 0 | 9 | 65 | 914 | 6.21 |
| PPARG | 15 | 10 | 2 | 0 | 32 | 439 | 0.96 |
| TP53 | 6 | 4 | 1 | 0 | 22 | 330 | 1.78 |

**Table S19. Which co-crystal ligand property predicts redocking failure**

| ligand property | Spearman rho vs RMSD | p (permutation) | templates |
|---|---|---|---|
| largest ring size | +0.490 | 0.002 | 42 |
| molecular weight | +0.320 | 0.049 | 42 |
| heavy atoms | +0.299 | 0.060 | 42 |
| rotatable bonds | +0.207 | 0.190 | 42 |
| number of rings | -0.045 | 0.778 | 42 |

**Table S20. Redocking validation rate by macrocyclicity**

| template group | templates | validated | validation rate | 95% CI |
|---|---|---|---|---|
| non-macrocyclic | 42 | 16 | 38.1% | 25.0% to 53.2% |
| macrocyclic | 10 | 0 | 0.0% | 0.0% to 27.8% |
| ligand not sanitizable | 10 | 1 | 10.0% |  |

**Table S21. Score variants by ROC-AUC against each matched design**

| negative set | score | roc_auc | vs registered consensus | targets better | targets worse | better unadjusted | targets |
|---|---|---|---|---|---|---|---|
| property-matched inactives | fusion: reciprocal rank | 0.545 | +0.004 | 0 | 0 | 0 | 5 |
| property-matched inactives | fusion: redocking-weighted z | 0.542 | +0.001 | 0 | 0 | 0 | 5 |
| property-matched inactives | fusion: rank mean | 0.540 | -0.001 | 0 | 0 | 0 | 5 |
| property-matched inactives | gnina alone | 0.539 | -0.001 | 0 | 0 | 0 | 5 |
| property-matched inactives | gnina (Vina-type) per heavy atom^(1/3) | 0.535 | -0.006 | 0 | 0 | 0 | 5 |
| property-matched inactives | AutoDock-GPU per heavy atom^(1/3) | 0.532 | -0.009 | 0 | 0 | 0 | 5 |
| property-matched inactives | AutoDock-GPU alone | 0.529 | -0.012 | 0 | 0 | 0 | 5 |
| property-matched inactives | gnina (CNN score) | 0.526 | -0.015 | 0 | 0 | 0 | 5 |
| property-matched inactives | gnina (CNN affinity) | 0.521 | -0.020 | 0 | 0 | 0 | 5 |
| property-matched inactives | gnina (Vina-type) per heavy atom | 0.508 | -0.033 | 0 | 0 | 0 | 5 |
| property-matched inactives | AutoDock-GPU per heavy atom | 0.502 | -0.039 | 0 | 0 | 0 | 5 |
| property-matched inactives | interaction recall | 0.485 | -0.056 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | gnina (CNN score) | 0.546 | +0.017 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | gnina alone | 0.540 | +0.011 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | gnina (Vina-type) per heavy atom^(1/3) | 0.535 | +0.005 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | fusion: redocking-weighted z | 0.532 | +0.002 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | fusion: reciprocal rank | 0.531 | +0.002 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | gnina (Vina-type) per heavy atom | 0.530 | +0.001 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | fusion: rank mean | 0.526 | -0.003 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | gnina (CNN affinity) | 0.516 | -0.013 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | AutoDock-GPU per heavy atom | 0.504 | -0.025 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | AutoDock-GPU per heavy atom^(1/3) | 0.499 | -0.030 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | AutoDock-GPU alone | 0.493 | -0.036 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | interaction recall | 0.486 | -0.043 | 0 | 0 | 0 | 5 |

**Table S22. Score variants by EF1 against each matched design**

| negative set | score | ef1 | vs registered consensus | targets better | targets worse | better unadjusted | targets |
|---|---|---|---|---|---|---|---|
| property-matched inactives | gnina (CNN score) | 3.45 | +1.50 | 0 | 0 | 0 | 5 |
| property-matched inactives | gnina alone | 2.40 | +0.45 | 0 | 0 | 0 | 5 |
| property-matched inactives | fusion: reciprocal rank | 1.95 | +0.00 | 0 | 0 | 0 | 5 |
| property-matched inactives | fusion: redocking-weighted z | 1.95 | +0.00 | 0 | 0 | 0 | 5 |
| property-matched inactives | AutoDock-GPU alone | 1.80 | -0.15 | 0 | 0 | 0 | 5 |
| property-matched inactives | fusion: rank mean | 1.50 | -0.45 | 0 | 0 | 0 | 5 |
| property-matched inactives | gnina (Vina-type) per heavy atom^(1/3) | 1.35 | -0.60 | 0 | 0 | 0 | 5 |
| property-matched inactives | AutoDock-GPU per heavy atom | 0.90 | -1.05 | 0 | 0 | 0 | 5 |
| property-matched inactives | interaction recall | 0.78 | -1.17 | 0 | 0 | 0 | 5 |
| property-matched inactives | gnina (CNN affinity) | 0.60 | -1.35 | 0 | 0 | 0 | 5 |
| property-matched inactives | AutoDock-GPU per heavy atom^(1/3) | 0.45 | -1.50 | 0 | 0 | 0 | 5 |
| property-matched inactives | gnina (Vina-type) per heavy atom | 0.45 | -1.50 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | gnina (CNN score) | 3.45 | +1.95 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | gnina alone | 2.70 | +1.20 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | fusion: rank mean | 1.95 | +0.45 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | fusion: redocking-weighted z | 1.95 | +0.45 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | fusion: reciprocal rank | 1.95 | +0.45 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | AutoDock-GPU alone | 1.50 | +0.00 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | gnina (Vina-type) per heavy atom^(1/3) | 1.35 | -0.15 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | AutoDock-GPU per heavy atom^(1/3) | 1.05 | -0.45 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | gnina (Vina-type) per heavy atom | 1.05 | -0.45 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | interaction recall | 0.94 | -0.56 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | AutoDock-GPU per heavy atom | 0.90 | -0.60 | 0 | 0 | 0 | 5 |
| DUD-E-protocol decoys | gnina (CNN affinity) | 0.60 | -0.90 | 0 | 0 | 0 | 5 |

**Table S23. Docking against each structure-free baseline separately**

| target | negative_set | baseline | score | baseline_auc | score_auc | diff | adj_low | adj_high | beats | below |
|---|---|---|---|---|---|---|---|---|---|---|
| ESR1_ant | random | n_heavy | autodock-gpu | 0.695 | 0.68 | -0.0147 | -0.0989 | 0.0588 | False | False |
| ESR1_ant | random | n_heavy | gnina | 0.695 | 0.606 | -0.0891 | -0.2679 | 0.0882 | False | False |
| ESR1_ant | random | n_heavy | consensus | 0.695 | 0.681 | -0.0136 | -0.143 | 0.0972 | False | False |
| ESR1_ant | random | mw | autodock-gpu | 0.696 | 0.68 | -0.0157 | -0.0971 | 0.0609 | False | False |
| ESR1_ant | random | mw | gnina | 0.696 | 0.606 | -0.0902 | -0.2621 | 0.0837 | False | False |
| ESR1_ant | random | mw | consensus | 0.696 | 0.681 | -0.0147 | -0.1381 | 0.0954 | False | False |
| ESR1_ant | random | clogp | autodock-gpu | 0.741 | 0.68 | -0.0606 | -0.184 | 0.0611 | False | False |
| ESR1_ant | random | clogp | gnina | 0.741 | 0.606 | -0.135 | -0.2991 | 0.0358 | False | False |
| ESR1_ant | random | clogp | consensus | 0.741 | 0.681 | -0.0595 | -0.1839 | 0.0758 | False | False |
| ESR1_ant | random | tpsa | autodock-gpu | 0.527 | 0.68 | 0.1534 | -0.0161 | 0.3174 | False | False |
| ESR1_ant | random | tpsa | gnina | 0.527 | 0.606 | 0.079 | -0.1425 | 0.2702 | False | False |
| ESR1_ant | random | tpsa | consensus | 0.527 | 0.681 | 0.1545 | -0.0361 | 0.3317 | False | False |
| ESR1_ant | random | hba | autodock-gpu | 0.497 | 0.68 | 0.1831 | 0.0363 | 0.3236 | True | False |
| ESR1_ant | random | hba | gnina | 0.497 | 0.606 | 0.1087 | -0.1194 | 0.3147 | False | False |
| ESR1_ant | random | hba | consensus | 0.497 | 0.681 | 0.1842 | 0.0026 | 0.3533 | True | False |
| ESR1_ant | random | hbd | autodock-gpu | 0.472 | 0.68 | 0.2086 | 0.014 | 0.3875 | True | False |
| ESR1_ant | random | hbd | gnina | 0.472 | 0.606 | 0.1342 | -0.0748 | 0.3184 | False | False |
| ESR1_ant | random | hbd | consensus | 0.472 | 0.681 | 0.2097 | 0.0112 | 0.3876 | True | False |
| ESR1_ant | random | rotb | autodock-gpu | 0.472 | 0.68 | 0.2084 | 0.0256 | 0.3814 | True | False |
| ESR1_ant | random | rotb | gnina | 0.472 | 0.606 | 0.134 | -0.068 | 0.3475 | False | False |
| ESR1_ant | random | rotb | consensus | 0.472 | 0.681 | 0.2095 | 0.0208 | 0.3802 | True | False |
| ESR1_ant | random | nn_sim | autodock-gpu | 0.685 | 0.68 | -0.0045 | -0.1923 | 0.1749 | False | False |
| ESR1_ant | random | nn_sim | gnina | 0.685 | 0.606 | -0.079 | -0.2348 | 0.0724 | False | False |
| ESR1_ant | random | nn_sim | consensus | 0.685 | 0.681 | -0.0035 | -0.1784 | 0.1595 | False | False |
| ESR1_ant | matched | n_heavy | autodock-gpu | 0.517 | 0.487 | -0.0302 | -0.122 | 0.0631 | False | False |
| ESR1_ant | matched | n_heavy | gnina | 0.517 | 0.491 | -0.0264 | -0.2287 | 0.1698 | False | False |
| ESR1_ant | matched | n_heavy | consensus | 0.517 | 0.495 | -0.0218 | -0.1643 | 0.1246 | False | False |
| ESR1_ant | matched | mw | autodock-gpu | 0.537 | 0.487 | -0.0506 | -0.1467 | 0.0453 | False | False |
| ESR1_ant | matched | mw | gnina | 0.537 | 0.491 | -0.0468 | -0.2464 | 0.151 | False | False |
| ESR1_ant | matched | mw | consensus | 0.537 | 0.495 | -0.0421 | -0.172 | 0.0917 | False | False |
| ESR1_ant | matched | clogp | autodock-gpu | 0.546 | 0.487 | -0.0587 | -0.1789 | 0.0565 | False | False |
| ESR1_ant | matched | clogp | gnina | 0.546 | 0.491 | -0.055 | -0.2302 | 0.1155 | False | False |
| ESR1_ant | matched | clogp | consensus | 0.546 | 0.495 | -0.0503 | -0.1778 | 0.0674 | False | False |
| ESR1_ant | matched | tpsa | autodock-gpu | 0.509 | 0.487 | -0.0225 | -0.2041 | 0.135 | False | False |
| ESR1_ant | matched | tpsa | gnina | 0.509 | 0.491 | -0.0187 | -0.2352 | 0.2041 | False | False |
| ESR1_ant | matched | tpsa | consensus | 0.509 | 0.495 | -0.014 | -0.2129 | 0.1774 | False | False |
| ESR1_ant | matched | hba | autodock-gpu | 0.481 | 0.487 | 0.0064 | -0.1351 | 0.1368 | False | False |
| ESR1_ant | matched | hba | gnina | 0.481 | 0.491 | 0.0101 | -0.2212 | 0.2269 | False | False |
| ESR1_ant | matched | hba | consensus | 0.481 | 0.495 | 0.0148 | -0.1688 | 0.1904 | False | False |
| ESR1_ant | matched | hbd | autodock-gpu | 0.51 | 0.487 | -0.0232 | -0.2181 | 0.1654 | False | False |
| ESR1_ant | matched | hbd | gnina | 0.51 | 0.491 | -0.0195 | -0.2019 | 0.1658 | False | False |
| ESR1_ant | matched | hbd | consensus | 0.51 | 0.495 | -0.0148 | -0.2221 | 0.1619 | False | False |
| ESR1_ant | matched | rotb | autodock-gpu | 0.506 | 0.487 | -0.019 | -0.1836 | 0.1445 | False | False |
| ESR1_ant | matched | rotb | gnina | 0.506 | 0.491 | -0.0152 | -0.2299 | 0.1917 | False | False |
| ESR1_ant | matched | rotb | consensus | 0.506 | 0.495 | -0.0105 | -0.1916 | 0.1814 | False | False |
| ESR1_ant | matched | nn_sim | autodock-gpu | 0.603 | 0.487 | -0.1166 | -0.3146 | 0.1142 | False | False |
| ESR1_ant | matched | nn_sim | gnina | 0.603 | 0.491 | -0.1128 | -0.2652 | 0.0423 | False | False |
| ESR1_ant | matched | nn_sim | consensus | 0.603 | 0.495 | -0.1081 | -0.2845 | 0.0799 | False | False |
| ESR1_ant | decoy | n_heavy | autodock-gpu | 0.452 | 0.427 | -0.0257 | -0.1256 | 0.0704 | False | False |
| ESR1_ant | decoy | n_heavy | gnina | 0.452 | 0.542 | 0.0898 | -0.1193 | 0.2827 | False | False |
| ESR1_ant | decoy | n_heavy | consensus | 0.452 | 0.489 | 0.0363 | -0.1164 | 0.1815 | False | False |
| ESR1_ant | decoy | mw | autodock-gpu | 0.498 | 0.427 | -0.0713 | -0.1687 | 0.0293 | False | False |
| ESR1_ant | decoy | mw | gnina | 0.498 | 0.542 | 0.0441 | -0.1546 | 0.2348 | False | False |
| ESR1_ant | decoy | mw | consensus | 0.498 | 0.489 | -0.0093 | -0.1558 | 0.1275 | False | False |
| ESR1_ant | decoy | clogp | autodock-gpu | 0.554 | 0.427 | -0.1268 | -0.237 | -0.0042 | False | True |
| ESR1_ant | decoy | clogp | gnina | 0.554 | 0.542 | -0.0113 | -0.184 | 0.1573 | False | False |
| ESR1_ant | decoy | clogp | consensus | 0.554 | 0.489 | -0.0648 | -0.2 | 0.067 | False | False |
| ESR1_ant | decoy | tpsa | autodock-gpu | 0.502 | 0.427 | -0.0748 | -0.2575 | 0.0974 | False | False |
| ESR1_ant | decoy | tpsa | gnina | 0.502 | 0.542 | 0.0406 | -0.1738 | 0.2615 | False | False |
| ESR1_ant | decoy | tpsa | consensus | 0.502 | 0.489 | -0.0129 | -0.2154 | 0.186 | False | False |
| ESR1_ant | decoy | hba | autodock-gpu | 0.451 | 0.427 | -0.024 | -0.1679 | 0.1137 | False | False |
| ESR1_ant | decoy | hba | gnina | 0.451 | 0.542 | 0.0914 | -0.1331 | 0.306 | False | False |
| ESR1_ant | decoy | hba | consensus | 0.451 | 0.489 | 0.038 | -0.1393 | 0.2297 | False | False |
| ESR1_ant | decoy | hbd | autodock-gpu | 0.475 | 0.427 | -0.0481 | -0.2427 | 0.1587 | False | False |
| ESR1_ant | decoy | hbd | gnina | 0.475 | 0.542 | 0.0673 | -0.1191 | 0.251 | False | False |
| ESR1_ant | decoy | hbd | consensus | 0.475 | 0.489 | 0.0139 | -0.1934 | 0.2187 | False | False |
| ESR1_ant | decoy | rotb | autodock-gpu | 0.467 | 0.427 | -0.0403 | -0.2165 | 0.1404 | False | False |
| ESR1_ant | decoy | rotb | gnina | 0.467 | 0.542 | 0.0752 | -0.1425 | 0.2802 | False | False |
| ESR1_ant | decoy | rotb | consensus | 0.467 | 0.489 | 0.0217 | -0.1723 | 0.2326 | False | False |
| ESR1_ant | decoy | nn_sim | autodock-gpu | 0.745 | 0.427 | -0.3187 | -0.5115 | -0.1043 | False | True |
| ESR1_ant | decoy | nn_sim | gnina | 0.745 | 0.542 | -0.2032 | -0.3615 | -0.0383 | False | True |
| ESR1_ant | decoy | nn_sim | consensus | 0.745 | 0.489 | -0.2567 | -0.4363 | -0.0646 | False | True |
| MAPK1 | random | n_heavy | autodock-gpu | 0.554 | 0.565 | 0.0111 | -0.0971 | 0.107 | False | False |
| MAPK1 | random | n_heavy | gnina | 0.554 | 0.58 | 0.0264 | -0.154 | 0.2045 | False | False |
| MAPK1 | random | n_heavy | consensus | 0.554 | 0.591 | 0.0372 | -0.1045 | 0.1687 | False | False |
| MAPK1 | random | mw | autodock-gpu | 0.541 | 0.565 | 0.024 | -0.0981 | 0.133 | False | False |
| MAPK1 | random | mw | gnina | 0.541 | 0.58 | 0.0393 | -0.1465 | 0.2236 | False | False |
| MAPK1 | random | mw | consensus | 0.541 | 0.591 | 0.0501 | -0.1011 | 0.1941 | False | False |
| MAPK1 | random | clogp | autodock-gpu | 0.697 | 0.565 | -0.1318 | -0.2676 | 0.0141 | False | False |
| MAPK1 | random | clogp | gnina | 0.697 | 0.58 | -0.1165 | -0.283 | 0.0477 | False | False |
| MAPK1 | random | clogp | consensus | 0.697 | 0.591 | -0.1057 | -0.2649 | 0.0493 | False | False |
| MAPK1 | random | tpsa | autodock-gpu | 0.419 | 0.565 | 0.1459 | -0.0132 | 0.2942 | False | False |
| MAPK1 | random | tpsa | gnina | 0.419 | 0.58 | 0.1612 | -0.0241 | 0.3337 | False | False |
| MAPK1 | random | tpsa | consensus | 0.419 | 0.591 | 0.1721 | 0.0041 | 0.3282 | True | False |
| MAPK1 | random | hba | autodock-gpu | 0.536 | 0.565 | 0.0284 | -0.1201 | 0.1768 | False | False |
| MAPK1 | random | hba | gnina | 0.536 | 0.58 | 0.0437 | -0.1405 | 0.2265 | False | False |
| MAPK1 | random | hba | consensus | 0.536 | 0.591 | 0.0546 | -0.1108 | 0.2098 | False | False |
| MAPK1 | random | hbd | autodock-gpu | 0.46 | 0.565 | 0.1049 | -0.1014 | 0.3008 | False | False |
| MAPK1 | random | hbd | gnina | 0.46 | 0.58 | 0.1202 | -0.074 | 0.2866 | False | False |
| MAPK1 | random | hbd | consensus | 0.46 | 0.591 | 0.1311 | -0.0805 | 0.3266 | False | False |
| MAPK1 | random | rotb | autodock-gpu | 0.469 | 0.565 | 0.096 | -0.0468 | 0.2439 | False | False |
| MAPK1 | random | rotb | gnina | 0.469 | 0.58 | 0.1114 | -0.0531 | 0.2764 | False | False |
| MAPK1 | random | rotb | consensus | 0.469 | 0.591 | 0.1222 | -0.0269 | 0.2791 | False | False |
| MAPK1 | random | nn_sim | autodock-gpu | 0.568 | 0.565 | -0.0034 | -0.1407 | 0.1334 | False | False |
| MAPK1 | random | nn_sim | gnina | 0.568 | 0.58 | 0.0119 | -0.1955 | 0.2041 | False | False |
| MAPK1 | random | nn_sim | consensus | 0.568 | 0.591 | 0.0227 | -0.1537 | 0.1837 | False | False |
| MAPK1 | matched | n_heavy | autodock-gpu | 0.503 | 0.506 | 0.0028 | -0.1168 | 0.1076 | False | False |
| MAPK1 | matched | n_heavy | gnina | 0.503 | 0.558 | 0.0545 | -0.1439 | 0.2346 | False | False |
| MAPK1 | matched | n_heavy | consensus | 0.503 | 0.556 | 0.0522 | -0.0891 | 0.1914 | False | False |
| MAPK1 | matched | mw | autodock-gpu | 0.502 | 0.506 | 0.0041 | -0.128 | 0.1282 | False | False |
| MAPK1 | matched | mw | gnina | 0.502 | 0.558 | 0.0558 | -0.1379 | 0.2468 | False | False |
| MAPK1 | matched | mw | consensus | 0.502 | 0.556 | 0.0535 | -0.1083 | 0.2094 | False | False |
| MAPK1 | matched | clogp | autodock-gpu | 0.516 | 0.506 | -0.0096 | -0.1614 | 0.1507 | False | False |
| MAPK1 | matched | clogp | gnina | 0.516 | 0.558 | 0.0421 | -0.1364 | 0.2048 | False | False |
| MAPK1 | matched | clogp | consensus | 0.516 | 0.556 | 0.0398 | -0.1326 | 0.2081 | False | False |
| MAPK1 | matched | tpsa | autodock-gpu | 0.483 | 0.506 | 0.0235 | -0.1404 | 0.183 | False | False |
| MAPK1 | matched | tpsa | gnina | 0.483 | 0.558 | 0.0752 | -0.1124 | 0.2658 | False | False |
| MAPK1 | matched | tpsa | consensus | 0.483 | 0.556 | 0.0729 | -0.1186 | 0.2462 | False | False |
| MAPK1 | matched | hba | autodock-gpu | 0.502 | 0.506 | 0.004 | -0.1502 | 0.1727 | False | False |
| MAPK1 | matched | hba | gnina | 0.502 | 0.558 | 0.0557 | -0.1171 | 0.233 | False | False |
| MAPK1 | matched | hba | consensus | 0.502 | 0.556 | 0.0534 | -0.1123 | 0.2263 | False | False |
| MAPK1 | matched | hbd | autodock-gpu | 0.502 | 0.506 | 0.0037 | -0.2046 | 0.2118 | False | False |
| MAPK1 | matched | hbd | gnina | 0.502 | 0.558 | 0.0554 | -0.1327 | 0.2416 | False | False |
| MAPK1 | matched | hbd | consensus | 0.502 | 0.556 | 0.0531 | -0.1504 | 0.2518 | False | False |
| MAPK1 | matched | rotb | autodock-gpu | 0.499 | 0.506 | 0.0076 | -0.1476 | 0.1609 | False | False |
| MAPK1 | matched | rotb | gnina | 0.499 | 0.558 | 0.0593 | -0.0926 | 0.222 | False | False |
| MAPK1 | matched | rotb | consensus | 0.499 | 0.556 | 0.057 | -0.0876 | 0.192 | False | False |
| MAPK1 | matched | nn_sim | autodock-gpu | 0.499 | 0.506 | 0.007 | -0.1346 | 0.1341 | False | False |
| MAPK1 | matched | nn_sim | gnina | 0.499 | 0.558 | 0.0587 | -0.1476 | 0.2657 | False | False |
| MAPK1 | matched | nn_sim | consensus | 0.499 | 0.556 | 0.0564 | -0.1216 | 0.2288 | False | False |
| MAPK1 | decoy | n_heavy | autodock-gpu | 0.537 | 0.505 | -0.0319 | -0.1521 | 0.0715 | False | False |
| MAPK1 | decoy | n_heavy | gnina | 0.537 | 0.568 | 0.0307 | -0.1702 | 0.2191 | False | False |
| MAPK1 | decoy | n_heavy | consensus | 0.537 | 0.558 | 0.0211 | -0.1287 | 0.1775 | False | False |
| MAPK1 | decoy | mw | autodock-gpu | 0.498 | 0.505 | 0.0066 | -0.1345 | 0.1284 | False | False |
| MAPK1 | decoy | mw | gnina | 0.498 | 0.568 | 0.0692 | -0.1326 | 0.2681 | False | False |
| MAPK1 | decoy | mw | consensus | 0.498 | 0.558 | 0.0595 | -0.1059 | 0.223 | False | False |
| MAPK1 | decoy | clogp | autodock-gpu | 0.511 | 0.505 | -0.0058 | -0.1628 | 0.1563 | False | False |
| MAPK1 | decoy | clogp | gnina | 0.511 | 0.568 | 0.0568 | -0.1275 | 0.2306 | False | False |
| MAPK1 | decoy | clogp | consensus | 0.511 | 0.558 | 0.0471 | -0.115 | 0.2168 | False | False |
| MAPK1 | decoy | tpsa | autodock-gpu | 0.5 | 0.505 | 0.0046 | -0.165 | 0.1632 | False | False |
| MAPK1 | decoy | tpsa | gnina | 0.5 | 0.568 | 0.0673 | -0.1255 | 0.2717 | False | False |
| MAPK1 | decoy | tpsa | consensus | 0.5 | 0.558 | 0.0576 | -0.1316 | 0.2414 | False | False |
| MAPK1 | decoy | hba | autodock-gpu | 0.54 | 0.505 | -0.0347 | -0.2011 | 0.1295 | False | False |
| MAPK1 | decoy | hba | gnina | 0.54 | 0.568 | 0.0279 | -0.1571 | 0.2158 | False | False |
| MAPK1 | decoy | hba | consensus | 0.54 | 0.558 | 0.0182 | -0.1601 | 0.1952 | False | False |
| MAPK1 | decoy | hbd | autodock-gpu | 0.501 | 0.505 | 0.0042 | -0.1968 | 0.2151 | False | False |
| MAPK1 | decoy | hbd | gnina | 0.501 | 0.568 | 0.0668 | -0.1156 | 0.2503 | False | False |
| MAPK1 | decoy | hbd | consensus | 0.501 | 0.558 | 0.0571 | -0.1562 | 0.2665 | False | False |
| MAPK1 | decoy | rotb | autodock-gpu | 0.49 | 0.505 | 0.0146 | -0.1458 | 0.1741 | False | False |
| MAPK1 | decoy | rotb | gnina | 0.49 | 0.568 | 0.0773 | -0.0782 | 0.2417 | False | False |
| MAPK1 | decoy | rotb | consensus | 0.49 | 0.558 | 0.0676 | -0.0823 | 0.2193 | False | False |
| MAPK1 | decoy | nn_sim | autodock-gpu | 0.621 | 0.505 | -0.1164 | -0.2493 | 0.0129 | False | False |
| MAPK1 | decoy | nn_sim | gnina | 0.621 | 0.568 | -0.0537 | -0.2515 | 0.1508 | False | False |
| MAPK1 | decoy | nn_sim | consensus | 0.621 | 0.558 | -0.0634 | -0.2338 | 0.1039 | False | False |
| MTORC1 | random | n_heavy | autodock-gpu | 0.494 | 0.533 | 0.0395 | -0.0683 | 0.1352 | False | False |
| MTORC1 | random | n_heavy | gnina | 0.494 | 0.473 | -0.0209 | -0.2225 | 0.1856 | False | False |
| MTORC1 | random | n_heavy | consensus | 0.494 | 0.507 | 0.0127 | -0.1406 | 0.1546 | False | False |
| MTORC1 | random | mw | autodock-gpu | 0.505 | 0.533 | 0.0288 | -0.0802 | 0.1291 | False | False |
| MTORC1 | random | mw | gnina | 0.505 | 0.473 | -0.0316 | -0.2409 | 0.1742 | False | False |
| MTORC1 | random | mw | consensus | 0.505 | 0.507 | 0.002 | -0.1635 | 0.1484 | False | False |
| MTORC1 | random | clogp | autodock-gpu | 0.507 | 0.533 | 0.0263 | -0.1426 | 0.205 | False | False |
| MTORC1 | random | clogp | gnina | 0.507 | 0.473 | -0.0341 | -0.2193 | 0.1209 | False | False |
| MTORC1 | random | clogp | consensus | 0.507 | 0.507 | -0.0004 | -0.1616 | 0.152 | False | False |
| MTORC1 | random | tpsa | autodock-gpu | 0.439 | 0.533 | 0.0943 | -0.0738 | 0.2599 | False | False |
| MTORC1 | random | tpsa | gnina | 0.439 | 0.473 | 0.034 | -0.1741 | 0.2284 | False | False |
| MTORC1 | random | tpsa | consensus | 0.439 | 0.507 | 0.0676 | -0.1049 | 0.2353 | False | False |
| MTORC1 | random | hba | autodock-gpu | 0.47 | 0.533 | 0.063 | -0.0968 | 0.2247 | False | False |
| MTORC1 | random | hba | gnina | 0.47 | 0.473 | 0.0026 | -0.178 | 0.1721 | False | False |
| MTORC1 | random | hba | consensus | 0.47 | 0.507 | 0.0362 | -0.1071 | 0.1926 | False | False |
| MTORC1 | random | hbd | autodock-gpu | 0.535 | 0.533 | -0.0018 | -0.1737 | 0.1607 | False | False |
| MTORC1 | random | hbd | gnina | 0.535 | 0.473 | -0.0622 | -0.2593 | 0.1361 | False | False |
| MTORC1 | random | hbd | consensus | 0.535 | 0.507 | -0.0285 | -0.2174 | 0.1541 | False | False |
| MTORC1 | random | rotb | autodock-gpu | 0.525 | 0.533 | 0.0084 | -0.159 | 0.1777 | False | False |
| MTORC1 | random | rotb | gnina | 0.525 | 0.473 | -0.052 | -0.2258 | 0.1337 | False | False |
| MTORC1 | random | rotb | consensus | 0.525 | 0.507 | -0.0183 | -0.1747 | 0.1418 | False | False |
| MTORC1 | random | nn_sim | autodock-gpu | 0.685 | 0.533 | -0.1513 | -0.3279 | 0.0204 | False | False |
| MTORC1 | random | nn_sim | gnina | 0.685 | 0.473 | -0.2117 | -0.3965 | -0.0204 | False | True |
| MTORC1 | random | nn_sim | consensus | 0.685 | 0.507 | -0.178 | -0.3427 | -0.0192 | False | True |
| MTORC1 | matched | n_heavy | autodock-gpu | 0.469 | 0.507 | 0.0378 | -0.074 | 0.1394 | False | False |
| MTORC1 | matched | n_heavy | gnina | 0.469 | 0.47 | 0.0013 | -0.2008 | 0.2077 | False | False |
| MTORC1 | matched | n_heavy | consensus | 0.469 | 0.474 | 0.0046 | -0.1474 | 0.1453 | False | False |
| MTORC1 | matched | mw | autodock-gpu | 0.499 | 0.507 | 0.0074 | -0.1021 | 0.1155 | False | False |
| MTORC1 | matched | mw | gnina | 0.499 | 0.47 | -0.0291 | -0.2383 | 0.1828 | False | False |
| MTORC1 | matched | mw | consensus | 0.499 | 0.474 | -0.0258 | -0.1838 | 0.1197 | False | False |
| MTORC1 | matched | clogp | autodock-gpu | 0.526 | 0.507 | -0.0193 | -0.1792 | 0.1604 | False | False |
| MTORC1 | matched | clogp | gnina | 0.526 | 0.47 | -0.0557 | -0.2172 | 0.1116 | False | False |
| MTORC1 | matched | clogp | consensus | 0.526 | 0.474 | -0.0525 | -0.1946 | 0.0984 | False | False |
| MTORC1 | matched | tpsa | autodock-gpu | 0.479 | 0.507 | 0.0279 | -0.1179 | 0.1871 | False | False |
| MTORC1 | matched | tpsa | gnina | 0.479 | 0.47 | -0.0086 | -0.2158 | 0.1891 | False | False |
| MTORC1 | matched | tpsa | consensus | 0.479 | 0.474 | -0.0053 | -0.1747 | 0.1609 | False | False |
| MTORC1 | matched | hba | autodock-gpu | 0.485 | 0.507 | 0.0215 | -0.1247 | 0.176 | False | False |
| MTORC1 | matched | hba | gnina | 0.485 | 0.47 | -0.015 | -0.189 | 0.165 | False | False |
| MTORC1 | matched | hba | consensus | 0.485 | 0.474 | -0.0117 | -0.1591 | 0.137 | False | False |
| MTORC1 | matched | hbd | autodock-gpu | 0.495 | 0.507 | 0.0121 | -0.1573 | 0.1835 | False | False |
| MTORC1 | matched | hbd | gnina | 0.495 | 0.47 | -0.0243 | -0.2255 | 0.1791 | False | False |
| MTORC1 | matched | hbd | consensus | 0.495 | 0.474 | -0.0211 | -0.2135 | 0.1472 | False | False |
| MTORC1 | matched | rotb | autodock-gpu | 0.504 | 0.507 | 0.0027 | -0.1631 | 0.1791 | False | False |
| MTORC1 | matched | rotb | gnina | 0.504 | 0.47 | -0.0338 | -0.2093 | 0.1466 | False | False |
| MTORC1 | matched | rotb | consensus | 0.504 | 0.474 | -0.0305 | -0.189 | 0.1359 | False | False |
| MTORC1 | matched | nn_sim | autodock-gpu | 0.633 | 0.507 | -0.1261 | -0.2878 | 0.0415 | False | False |
| MTORC1 | matched | nn_sim | gnina | 0.633 | 0.47 | -0.1626 | -0.3371 | 0.039 | False | False |
| MTORC1 | matched | nn_sim | consensus | 0.633 | 0.474 | -0.1593 | -0.3121 | -0.0123 | False | True |
| MTORC1 | decoy | n_heavy | autodock-gpu | 0.483 | 0.586 | 0.1035 | -0.0094 | 0.208 | False | False |
| MTORC1 | decoy | n_heavy | gnina | 0.483 | 0.534 | 0.0512 | -0.1504 | 0.2597 | False | False |
| MTORC1 | decoy | n_heavy | consensus | 0.483 | 0.587 | 0.1044 | -0.0406 | 0.244 | False | False |
| MTORC1 | decoy | mw | autodock-gpu | 0.502 | 0.586 | 0.084 | -0.0254 | 0.1921 | False | False |
| MTORC1 | decoy | mw | gnina | 0.502 | 0.534 | 0.0317 | -0.167 | 0.2499 | False | False |
| MTORC1 | decoy | mw | consensus | 0.502 | 0.587 | 0.0849 | -0.0645 | 0.2377 | False | False |
| MTORC1 | decoy | clogp | autodock-gpu | 0.563 | 0.586 | 0.0232 | -0.1281 | 0.1974 | False | False |
| MTORC1 | decoy | clogp | gnina | 0.563 | 0.534 | -0.0291 | -0.1948 | 0.147 | False | False |
| MTORC1 | decoy | clogp | consensus | 0.563 | 0.587 | 0.0241 | -0.1255 | 0.167 | False | False |
| MTORC1 | decoy | tpsa | autodock-gpu | 0.498 | 0.586 | 0.0883 | -0.061 | 0.2373 | False | False |
| MTORC1 | decoy | tpsa | gnina | 0.498 | 0.534 | 0.036 | -0.1756 | 0.2295 | False | False |
| MTORC1 | decoy | tpsa | consensus | 0.498 | 0.587 | 0.0892 | -0.0799 | 0.248 | False | False |
| MTORC1 | decoy | hba | autodock-gpu | 0.481 | 0.586 | 0.1049 | -0.0452 | 0.2555 | False | False |
| MTORC1 | decoy | hba | gnina | 0.481 | 0.534 | 0.0526 | -0.1328 | 0.2415 | False | False |
| MTORC1 | decoy | hba | consensus | 0.481 | 0.587 | 0.1058 | -0.0404 | 0.2577 | False | False |
| MTORC1 | decoy | hbd | autodock-gpu | 0.492 | 0.586 | 0.0936 | -0.0751 | 0.2615 | False | False |
| MTORC1 | decoy | hbd | gnina | 0.492 | 0.534 | 0.0413 | -0.17 | 0.2508 | False | False |
| MTORC1 | decoy | hbd | consensus | 0.492 | 0.587 | 0.0945 | -0.1018 | 0.2833 | False | False |
| MTORC1 | decoy | rotb | autodock-gpu | 0.518 | 0.586 | 0.0678 | -0.0951 | 0.2238 | False | False |
| MTORC1 | decoy | rotb | gnina | 0.518 | 0.534 | 0.0155 | -0.1649 | 0.201 | False | False |
| MTORC1 | decoy | rotb | consensus | 0.518 | 0.587 | 0.0687 | -0.0916 | 0.2417 | False | False |
| MTORC1 | decoy | nn_sim | autodock-gpu | 0.811 | 0.586 | -0.2246 | -0.3694 | -0.0738 | False | True |
| MTORC1 | decoy | nn_sim | gnina | 0.811 | 0.534 | -0.2769 | -0.4442 | -0.0924 | False | True |
| MTORC1 | decoy | nn_sim | consensus | 0.811 | 0.587 | -0.2237 | -0.3683 | -0.064 | False | True |
| PPARG | random | n_heavy | autodock-gpu | 0.767 | 0.773 | 0.0054 | -0.0854 | 0.0818 | False | False |
| PPARG | random | n_heavy | gnina | 0.767 | 0.783 | 0.0156 | -0.1372 | 0.1451 | False | False |
| PPARG | random | n_heavy | consensus | 0.767 | 0.809 | 0.0418 | -0.0622 | 0.132 | False | False |
| PPARG | random | mw | autodock-gpu | 0.744 | 0.773 | 0.0285 | -0.0733 | 0.117 | False | False |
| PPARG | random | mw | gnina | 0.744 | 0.783 | 0.0387 | -0.1166 | 0.1665 | False | False |
| PPARG | random | mw | consensus | 0.744 | 0.809 | 0.0649 | -0.0451 | 0.1548 | False | False |
| PPARG | random | clogp | autodock-gpu | 0.682 | 0.773 | 0.0909 | -0.0259 | 0.2058 | False | False |
| PPARG | random | clogp | gnina | 0.682 | 0.783 | 0.1011 | -0.0409 | 0.2358 | False | False |
| PPARG | random | clogp | consensus | 0.682 | 0.809 | 0.1273 | 0.0177 | 0.2402 | True | False |
| PPARG | random | tpsa | autodock-gpu | 0.592 | 0.773 | 0.1809 | 0.0253 | 0.3093 | True | False |
| PPARG | random | tpsa | gnina | 0.592 | 0.783 | 0.191 | -0.0015 | 0.3717 | False | False |
| PPARG | random | tpsa | consensus | 0.592 | 0.809 | 0.2172 | 0.0581 | 0.3675 | True | False |
| PPARG | random | hba | autodock-gpu | 0.658 | 0.773 | 0.1154 | -0.0228 | 0.256 | False | False |
| PPARG | random | hba | gnina | 0.658 | 0.783 | 0.1256 | -0.0855 | 0.3083 | False | False |
| PPARG | random | hba | consensus | 0.658 | 0.809 | 0.1517 | -0.0122 | 0.3021 | False | False |
| PPARG | random | hbd | autodock-gpu | 0.275 | 0.773 | 0.4977 | 0.3239 | 0.6342 | True | False |
| PPARG | random | hbd | gnina | 0.275 | 0.783 | 0.5079 | 0.3075 | 0.676 | True | False |
| PPARG | random | hbd | consensus | 0.275 | 0.809 | 0.5341 | 0.3631 | 0.6904 | True | False |
| PPARG | random | rotb | autodock-gpu | 0.688 | 0.773 | 0.0847 | -0.097 | 0.2508 | False | False |
| PPARG | random | rotb | gnina | 0.688 | 0.783 | 0.0949 | -0.0808 | 0.2937 | False | False |
| PPARG | random | rotb | consensus | 0.688 | 0.809 | 0.121 | -0.0423 | 0.2878 | False | False |
| PPARG | random | nn_sim | autodock-gpu | 0.75 | 0.773 | 0.0231 | -0.1055 | 0.1793 | False | False |
| PPARG | random | nn_sim | gnina | 0.75 | 0.783 | 0.0333 | -0.1494 | 0.219 | False | False |
| PPARG | random | nn_sim | consensus | 0.75 | 0.809 | 0.0594 | -0.0821 | 0.2093 | False | False |
| PPARG | matched | n_heavy | autodock-gpu | 0.566 | 0.581 | 0.0154 | -0.0927 | 0.1111 | False | False |
| PPARG | matched | n_heavy | gnina | 0.566 | 0.624 | 0.0584 | -0.1289 | 0.2304 | False | False |
| PPARG | matched | n_heavy | consensus | 0.566 | 0.599 | 0.0334 | -0.1196 | 0.1651 | False | False |
| PPARG | matched | mw | autodock-gpu | 0.543 | 0.581 | 0.0376 | -0.0845 | 0.1469 | False | False |
| PPARG | matched | mw | gnina | 0.543 | 0.624 | 0.0806 | -0.1126 | 0.2507 | False | False |
| PPARG | matched | mw | consensus | 0.543 | 0.599 | 0.0556 | -0.103 | 0.1922 | False | False |
| PPARG | matched | clogp | autodock-gpu | 0.556 | 0.581 | 0.0245 | -0.1342 | 0.1832 | False | False |
| PPARG | matched | clogp | gnina | 0.556 | 0.624 | 0.0675 | -0.0996 | 0.2402 | False | False |
| PPARG | matched | clogp | consensus | 0.556 | 0.599 | 0.0424 | -0.1081 | 0.1995 | False | False |
| PPARG | matched | tpsa | autodock-gpu | 0.549 | 0.581 | 0.0323 | -0.1284 | 0.1854 | False | False |
| PPARG | matched | tpsa | gnina | 0.549 | 0.624 | 0.0753 | -0.1501 | 0.2906 | False | False |
| PPARG | matched | tpsa | consensus | 0.549 | 0.599 | 0.0502 | -0.1503 | 0.248 | False | False |
| PPARG | matched | hba | autodock-gpu | 0.532 | 0.581 | 0.049 | -0.1243 | 0.219 | False | False |
| PPARG | matched | hba | gnina | 0.532 | 0.624 | 0.092 | -0.1343 | 0.3133 | False | False |
| PPARG | matched | hba | consensus | 0.532 | 0.599 | 0.067 | -0.1427 | 0.2679 | False | False |
| PPARG | matched | hbd | autodock-gpu | 0.463 | 0.581 | 0.1177 | -0.0902 | 0.3127 | False | False |
| PPARG | matched | hbd | gnina | 0.463 | 0.624 | 0.1607 | -0.0857 | 0.388 | False | False |
| PPARG | matched | hbd | consensus | 0.463 | 0.599 | 0.1356 | -0.0955 | 0.3556 | False | False |
| PPARG | matched | rotb | autodock-gpu | 0.517 | 0.581 | 0.064 | -0.1481 | 0.2552 | False | False |
| PPARG | matched | rotb | gnina | 0.517 | 0.624 | 0.107 | -0.1017 | 0.327 | False | False |
| PPARG | matched | rotb | consensus | 0.517 | 0.599 | 0.082 | -0.1279 | 0.2856 | False | False |
| PPARG | matched | nn_sim | autodock-gpu | 0.656 | 0.581 | -0.0748 | -0.2225 | 0.089 | False | False |
| PPARG | matched | nn_sim | gnina | 0.656 | 0.624 | -0.0318 | -0.256 | 0.2048 | False | False |
| PPARG | matched | nn_sim | consensus | 0.656 | 0.599 | -0.0568 | -0.2341 | 0.1456 | False | False |
| PPARG | decoy | n_heavy | autodock-gpu | 0.503 | 0.482 | -0.0213 | -0.1402 | 0.074 | False | False |
| PPARG | decoy | n_heavy | gnina | 0.503 | 0.546 | 0.0435 | -0.1366 | 0.2278 | False | False |
| PPARG | decoy | n_heavy | consensus | 0.503 | 0.516 | 0.0131 | -0.1419 | 0.1432 | False | False |
| PPARG | decoy | mw | autodock-gpu | 0.498 | 0.482 | -0.0165 | -0.1507 | 0.0991 | False | False |
| PPARG | decoy | mw | gnina | 0.498 | 0.546 | 0.0483 | -0.1342 | 0.2295 | False | False |
| PPARG | decoy | mw | consensus | 0.498 | 0.516 | 0.0179 | -0.1433 | 0.1597 | False | False |
| PPARG | decoy | clogp | autodock-gpu | 0.579 | 0.482 | -0.0978 | -0.2658 | 0.0558 | False | False |
| PPARG | decoy | clogp | gnina | 0.579 | 0.546 | -0.0331 | -0.1975 | 0.1492 | False | False |
| PPARG | decoy | clogp | consensus | 0.579 | 0.516 | -0.0634 | -0.2108 | 0.0878 | False | False |
| PPARG | decoy | tpsa | autodock-gpu | 0.499 | 0.482 | -0.0171 | -0.1708 | 0.1276 | False | False |
| PPARG | decoy | tpsa | gnina | 0.499 | 0.546 | 0.0477 | -0.1666 | 0.2731 | False | False |
| PPARG | decoy | tpsa | consensus | 0.499 | 0.516 | 0.0173 | -0.1851 | 0.2143 | False | False |
| PPARG | decoy | hba | autodock-gpu | 0.502 | 0.482 | -0.02 | -0.2065 | 0.1615 | False | False |
| PPARG | decoy | hba | gnina | 0.502 | 0.546 | 0.0448 | -0.191 | 0.2698 | False | False |
| PPARG | decoy | hba | consensus | 0.502 | 0.516 | 0.0144 | -0.2042 | 0.2248 | False | False |
| PPARG | decoy | hbd | autodock-gpu | 0.451 | 0.482 | 0.0306 | -0.172 | 0.2328 | False | False |
| PPARG | decoy | hbd | gnina | 0.451 | 0.546 | 0.0953 | -0.1558 | 0.3324 | False | False |
| PPARG | decoy | hbd | consensus | 0.451 | 0.516 | 0.065 | -0.1723 | 0.2933 | False | False |
| PPARG | decoy | rotb | autodock-gpu | 0.557 | 0.482 | -0.0752 | -0.2796 | 0.1066 | False | False |
| PPARG | decoy | rotb | gnina | 0.557 | 0.546 | -0.0105 | -0.2547 | 0.2118 | False | False |
| PPARG | decoy | rotb | consensus | 0.557 | 0.516 | -0.0408 | -0.2679 | 0.1662 | False | False |
| PPARG | decoy | nn_sim | autodock-gpu | 0.698 | 0.482 | -0.2161 | -0.3843 | -0.0486 | False | True |
| PPARG | decoy | nn_sim | gnina | 0.698 | 0.546 | -0.1514 | -0.3871 | 0.0964 | False | False |
| PPARG | decoy | nn_sim | consensus | 0.698 | 0.516 | -0.1818 | -0.3844 | 0.0148 | False | False |
| TP53 | random | n_heavy | autodock-gpu | 0.593 | 0.597 | 0.0042 | -0.1062 | 0.0941 | False | False |
| TP53 | random | n_heavy | gnina | 0.593 | 0.56 | -0.0323 | -0.213 | 0.143 | False | False |
| TP53 | random | n_heavy | consensus | 0.593 | 0.621 | 0.0286 | -0.107 | 0.1498 | False | False |
| TP53 | random | mw | autodock-gpu | 0.603 | 0.597 | -0.0065 | -0.1493 | 0.1148 | False | False |
| TP53 | random | mw | gnina | 0.603 | 0.56 | -0.043 | -0.2304 | 0.1422 | False | False |
| TP53 | random | mw | consensus | 0.603 | 0.621 | 0.0179 | -0.1448 | 0.1558 | False | False |
| TP53 | random | clogp | autodock-gpu | 0.475 | 0.597 | 0.1214 | -0.0263 | 0.2698 | False | False |
| TP53 | random | clogp | gnina | 0.475 | 0.56 | 0.0849 | -0.1015 | 0.2549 | False | False |
| TP53 | random | clogp | consensus | 0.475 | 0.621 | 0.1458 | -0.0188 | 0.2982 | False | False |
| TP53 | random | tpsa | autodock-gpu | 0.584 | 0.597 | 0.0127 | -0.162 | 0.1937 | False | False |
| TP53 | random | tpsa | gnina | 0.584 | 0.56 | -0.0237 | -0.2205 | 0.1659 | False | False |
| TP53 | random | tpsa | consensus | 0.584 | 0.621 | 0.0372 | -0.1462 | 0.222 | False | False |
| TP53 | random | hba | autodock-gpu | 0.597 | 0.597 | -0.0002 | -0.1712 | 0.1658 | False | False |
| TP53 | random | hba | gnina | 0.597 | 0.56 | -0.0367 | -0.2322 | 0.1521 | False | False |
| TP53 | random | hba | consensus | 0.597 | 0.621 | 0.0242 | -0.1541 | 0.2053 | False | False |
| TP53 | random | hbd | autodock-gpu | 0.583 | 0.597 | 0.0139 | -0.1638 | 0.2006 | False | False |
| TP53 | random | hbd | gnina | 0.583 | 0.56 | -0.0226 | -0.2089 | 0.1581 | False | False |
| TP53 | random | hbd | consensus | 0.583 | 0.621 | 0.0383 | -0.1423 | 0.2204 | False | False |
| TP53 | random | rotb | autodock-gpu | 0.435 | 0.597 | 0.162 | -0.0274 | 0.3363 | False | False |
| TP53 | random | rotb | gnina | 0.435 | 0.56 | 0.1255 | -0.0902 | 0.3299 | False | False |
| TP53 | random | rotb | consensus | 0.435 | 0.621 | 0.1864 | -0.0181 | 0.3819 | False | False |
| TP53 | random | nn_sim | autodock-gpu | 0.578 | 0.597 | 0.0191 | -0.1693 | 0.214 | False | False |
| TP53 | random | nn_sim | gnina | 0.578 | 0.56 | -0.0174 | -0.1656 | 0.1437 | False | False |
| TP53 | random | nn_sim | consensus | 0.578 | 0.621 | 0.0435 | -0.1218 | 0.2006 | False | False |
| TP53 | matched | n_heavy | autodock-gpu | 0.54 | 0.562 | 0.0219 | -0.1025 | 0.13 | False | False |
| TP53 | matched | n_heavy | gnina | 0.54 | 0.555 | 0.0143 | -0.18 | 0.1862 | False | False |
| TP53 | matched | n_heavy | consensus | 0.54 | 0.581 | 0.0407 | -0.1031 | 0.1699 | False | False |
| TP53 | matched | mw | autodock-gpu | 0.533 | 0.562 | 0.0291 | -0.1164 | 0.1588 | False | False |
| TP53 | matched | mw | gnina | 0.533 | 0.555 | 0.0214 | -0.1933 | 0.2062 | False | False |
| TP53 | matched | mw | consensus | 0.533 | 0.581 | 0.0479 | -0.1192 | 0.1936 | False | False |
| TP53 | matched | clogp | autodock-gpu | 0.489 | 0.562 | 0.0732 | -0.0651 | 0.2168 | False | False |
| TP53 | matched | clogp | gnina | 0.489 | 0.555 | 0.0656 | -0.1024 | 0.2472 | False | False |
| TP53 | matched | clogp | consensus | 0.489 | 0.581 | 0.0921 | -0.055 | 0.2343 | False | False |
| TP53 | matched | tpsa | autodock-gpu | 0.52 | 0.562 | 0.0426 | -0.1452 | 0.2275 | False | False |
| TP53 | matched | tpsa | gnina | 0.52 | 0.555 | 0.035 | -0.1492 | 0.2236 | False | False |
| TP53 | matched | tpsa | consensus | 0.52 | 0.581 | 0.0615 | -0.1332 | 0.2359 | False | False |
| TP53 | matched | hba | autodock-gpu | 0.53 | 0.562 | 0.0327 | -0.1477 | 0.2106 | False | False |
| TP53 | matched | hba | gnina | 0.53 | 0.555 | 0.0251 | -0.1687 | 0.2261 | False | False |
| TP53 | matched | hba | consensus | 0.53 | 0.581 | 0.0516 | -0.1349 | 0.2464 | False | False |
| TP53 | matched | hbd | autodock-gpu | 0.525 | 0.562 | 0.0376 | -0.1633 | 0.2276 | False | False |
| TP53 | matched | hbd | gnina | 0.525 | 0.555 | 0.03 | -0.156 | 0.1918 | False | False |
| TP53 | matched | hbd | consensus | 0.525 | 0.581 | 0.0564 | -0.1329 | 0.2261 | False | False |
| TP53 | matched | rotb | autodock-gpu | 0.478 | 0.562 | 0.0841 | -0.1197 | 0.2738 | False | False |
| TP53 | matched | rotb | gnina | 0.478 | 0.555 | 0.0765 | -0.1531 | 0.2916 | False | False |
| TP53 | matched | rotb | consensus | 0.478 | 0.581 | 0.103 | -0.0961 | 0.3091 | False | False |
| TP53 | matched | nn_sim | autodock-gpu | 0.563 | 0.562 | -0.0006 | -0.1739 | 0.1948 | False | False |
| TP53 | matched | nn_sim | gnina | 0.563 | 0.555 | -0.0082 | -0.1477 | 0.153 | False | False |
| TP53 | matched | nn_sim | consensus | 0.563 | 0.581 | 0.0182 | -0.1306 | 0.1816 | False | False |
| TP53 | decoy | n_heavy | autodock-gpu | 0.468 | 0.466 | -0.0016 | -0.1269 | 0.1156 | False | False |
| TP53 | decoy | n_heavy | gnina | 0.468 | 0.509 | 0.0417 | -0.1611 | 0.2263 | False | False |
| TP53 | decoy | n_heavy | consensus | 0.468 | 0.497 | 0.029 | -0.1404 | 0.1769 | False | False |
| TP53 | decoy | mw | autodock-gpu | 0.493 | 0.466 | -0.027 | -0.1858 | 0.1033 | False | False |
| TP53 | decoy | mw | gnina | 0.493 | 0.509 | 0.0163 | -0.1846 | 0.1982 | False | False |
| TP53 | decoy | mw | consensus | 0.493 | 0.497 | 0.0036 | -0.1694 | 0.1644 | False | False |
| TP53 | decoy | clogp | autodock-gpu | 0.471 | 0.466 | -0.0048 | -0.1496 | 0.1473 | False | False |
| TP53 | decoy | clogp | gnina | 0.471 | 0.509 | 0.0385 | -0.127 | 0.2151 | False | False |
| TP53 | decoy | clogp | consensus | 0.471 | 0.497 | 0.0259 | -0.1272 | 0.1671 | False | False |
| TP53 | decoy | tpsa | autodock-gpu | 0.507 | 0.466 | -0.041 | -0.2245 | 0.1442 | False | False |
| TP53 | decoy | tpsa | gnina | 0.507 | 0.509 | 0.0023 | -0.1882 | 0.1935 | False | False |
| TP53 | decoy | tpsa | consensus | 0.507 | 0.497 | -0.0104 | -0.2014 | 0.1645 | False | False |
| TP53 | decoy | hba | autodock-gpu | 0.515 | 0.466 | -0.0494 | -0.2344 | 0.1343 | False | False |
| TP53 | decoy | hba | gnina | 0.515 | 0.509 | -0.0061 | -0.2037 | 0.194 | False | False |
| TP53 | decoy | hba | consensus | 0.515 | 0.497 | -0.0187 | -0.2031 | 0.1737 | False | False |
| TP53 | decoy | hbd | autodock-gpu | 0.54 | 0.466 | -0.0738 | -0.2663 | 0.1154 | False | False |
| TP53 | decoy | hbd | gnina | 0.54 | 0.509 | -0.0305 | -0.2084 | 0.1417 | False | False |
| TP53 | decoy | hbd | consensus | 0.54 | 0.497 | -0.0431 | -0.2297 | 0.1372 | False | False |
| TP53 | decoy | rotb | autodock-gpu | 0.45 | 0.466 | 0.0163 | -0.179 | 0.2071 | False | False |
| TP53 | decoy | rotb | gnina | 0.45 | 0.509 | 0.0596 | -0.177 | 0.2896 | False | False |
| TP53 | decoy | rotb | consensus | 0.45 | 0.497 | 0.047 | -0.1654 | 0.2726 | False | False |
| TP53 | decoy | nn_sim | autodock-gpu | 0.665 | 0.466 | -0.199 | -0.381 | -0.002 | False | True |
| TP53 | decoy | nn_sim | gnina | 0.665 | 0.509 | -0.1557 | -0.305 | 0.0163 | False | False |
| TP53 | decoy | nn_sim | consensus | 0.665 | 0.497 | -0.1683 | -0.3186 | 0.002 | False | False |
