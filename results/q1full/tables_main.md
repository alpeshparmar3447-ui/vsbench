**Table 1. Panel, template validation and compounds docked**

| target | actives | random inactives | matched inactives | decoys | templates validated | pocket ions | deposition checked | selected | RMSD AutoDock-GPU (A) | RMSD gnina (A) | validated |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ESR1_ant | 102 | 300 | 816 | 816 | 10/15 | none | yes | 6b0f | 0.48 | 0.40 | yes |
| MAPK1 | 308 | 300 | 2464 | 2464 | 4/15 | none | yes | 4qte | 0.60 | 0.44 | yes |
| MTORC1 | 97 | 300 | 776 | 776 | 0/11 | none | yes | 4jsx | 6.13 | 6.21 | NO |
| PPARG | 27 | 300 | 216 | 216 | 2/15 | none | yes | 2p4y | 0.96 | 0.45 | yes |
| TP53 | 79 | 300 | 632 | 632 | 1/6 | none | yes | 4agq | 1.78 | 1.37 | yes |

**Table 2. ROC-AUC against each negative set**

| negative set | target | autodock-gpu | gnina | consensus |
|---|---|---|---|---|
| random inactives | ESR1_ant | 0.65 [0.58, 0.71] | 0.64 [0.57, 0.70] | 0.69 [0.62, 0.75] |
| random inactives | MAPK1 | 0.57 [0.53, 0.62] | 0.57 [0.53, 0.61] | 0.59 [0.55, 0.64] |
| random inactives | MTORC1 | 0.49 [0.42, 0.56] | 0.49 [0.43, 0.56] | 0.50 [0.43, 0.56] |
| random inactives | PPARG | 0.77 [0.69, 0.84] | 0.78 [0.68, 0.88] | 0.81 [0.72, 0.88] |
| random inactives | TP53 | 0.59 [0.52, 0.65] | 0.57 [0.50, 0.64] | 0.63 [0.56, 0.70] |
| property-matched inactives | ESR1_ant | 0.52 [0.46, 0.58] | 0.56 [0.49, 0.62] | 0.55 [0.49, 0.61] |
| property-matched inactives | MAPK1 | 0.51 [0.47, 0.54] | 0.55 [0.52, 0.59] | 0.54 [0.51, 0.57] |
| property-matched inactives | MTORC1 | 0.50 [0.43, 0.56] | 0.48 [0.42, 0.54] | 0.48 [0.43, 0.53] |
| property-matched inactives | PPARG | 0.58 [0.47, 0.69] | 0.62 [0.51, 0.73] | 0.60 [0.48, 0.71] |
| property-matched inactives | TP53 | 0.56 [0.49, 0.62] | 0.57 [0.50, 0.63] | 0.59 [0.52, 0.65] |
| DUD-E-protocol decoys | ESR1_ant | 0.43 [0.37, 0.49] | 0.57 [0.50, 0.63] | 0.50 [0.43, 0.56] |
| DUD-E-protocol decoys | MAPK1 | 0.50 [0.47, 0.54] | 0.56 [0.53, 0.59] | 0.54 [0.51, 0.58] |
| DUD-E-protocol decoys | MTORC1 | 0.59 [0.53, 0.64] | 0.55 [0.49, 0.60] | 0.60 [0.55, 0.65] |
| DUD-E-protocol decoys | PPARG | 0.49 [0.38, 0.61] | 0.56 [0.43, 0.69] | 0.53 [0.40, 0.65] |
| DUD-E-protocol decoys | TP53 | 0.42 [0.35, 0.49] | 0.50 [0.44, 0.56] | 0.46 [0.40, 0.53] |

**Table 3. ROC-AUC of single properties used as scores**

| target | decoy: mw | decoy: n_heavy | decoy: tpsa | matched: mw | matched: n_heavy | matched: tpsa | random: mw | random: n_heavy | random: tpsa |
|---|---|---|---|---|---|---|---|---|---|
| ESR1_ant | 0.49 | 0.45 | 0.5 | 0.55 | 0.54 | 0.52 | 0.68 | 0.68 | 0.53 |
| MAPK1 | 0.5 | 0.52 | 0.5 | 0.51 | 0.51 | 0.49 | 0.57 | 0.57 | 0.42 |
| MTORC1 | 0.5 | 0.5 | 0.5 | 0.5 | 0.48 | 0.48 | 0.5 | 0.49 | 0.4 |
| PPARG | 0.51 | 0.51 | 0.51 | 0.54 | 0.57 | 0.55 | 0.74 | 0.77 | 0.59 |
| TP53 | 0.49 | 0.46 | 0.5 | 0.57 | 0.57 | 0.52 | 0.61 | 0.61 | 0.6 |

**Table 4. Size-controlled ROC-AUC against each negative set**

| negative set | score | ROC-AUC | size-controlled AUC | change | SD of change | targets |
|---|---|---|---|---|---|---|
| random inactives | gnina | 0.611 | 0.603 | -0.008 | 0.034 | 5 |
| random inactives | fusion: rank mean | 0.644 | 0.593 | -0.051 | 0.045 | 5 |
| random inactives | gnina (CNN score) | 0.549 | 0.593 | +0.044 | 0.034 | 5 |
| random inactives | fusion: reciprocal rank | 0.644 | 0.593 | -0.051 | 0.062 | 5 |
| random inactives | fusion: redocking-weighted z | 0.643 | 0.590 | -0.053 | 0.047 | 5 |
| random inactives | consensus | 0.643 | 0.588 | -0.055 | 0.049 | 5 |
| random inactives | gnina (CNN affinity) | 0.625 | 0.567 | -0.058 | 0.062 | 5 |
| random inactives | gnina (Vina-type) per heavy atom^(1/3) | 0.574 | 0.565 | -0.008 | 0.022 | 5 |
| random inactives | gnina (Vina-type) | 0.607 | 0.565 | -0.042 | 0.056 | 5 |
| random inactives | gnina (Vina-type) per heavy atom | 0.444 | 0.553 | +0.109 | 0.106 | 5 |
| random inactives | autodock-gpu | 0.616 | 0.525 | -0.091 | 0.098 | 5 |
| random inactives | AutoDock-GPU per heavy atom^(1/3) | 0.587 | 0.522 | -0.065 | 0.074 | 5 |
| random inactives | autodock-gpu (largest cluster) | 0.613 | 0.522 | -0.091 | 0.095 | 5 |
| random inactives | interaction recall | 0.538 | 0.517 | -0.022 | 0.026 | 5 |
| random inactives | AutoDock-GPU per heavy atom | 0.430 | 0.512 | +0.083 | 0.085 | 5 |
| property-matched inactives | fusion: reciprocal rank | 0.559 | 0.565 | +0.006 | 0.020 | 5 |
| property-matched inactives | gnina | 0.556 | 0.560 | +0.004 | 0.017 | 5 |
| property-matched inactives | gnina (CNN score) | 0.540 | 0.554 | +0.014 | 0.011 | 5 |
| property-matched inactives | fusion: redocking-weighted z | 0.553 | 0.554 | +0.001 | 0.024 | 5 |
| property-matched inactives | fusion: rank mean | 0.552 | 0.553 | +0.001 | 0.023 | 5 |
| property-matched inactives | consensus | 0.552 | 0.551 | -0.001 | 0.024 | 5 |
| property-matched inactives | gnina (Vina-type) per heavy atom^(1/3) | 0.533 | 0.534 | +0.001 | 0.010 | 5 |
| property-matched inactives | gnina (CNN affinity) | 0.535 | 0.533 | -0.003 | 0.026 | 5 |
| property-matched inactives | gnina (Vina-type) | 0.533 | 0.532 | -0.001 | 0.020 | 5 |
| property-matched inactives | gnina (Vina-type) per heavy atom | 0.500 | 0.523 | +0.023 | 0.029 | 5 |
| property-matched inactives | autodock-gpu | 0.533 | 0.517 | -0.015 | 0.021 | 5 |
| property-matched inactives | AutoDock-GPU per heavy atom^(1/3) | 0.528 | 0.515 | -0.013 | 0.014 | 5 |
| property-matched inactives | autodock-gpu (largest cluster) | 0.527 | 0.513 | -0.014 | 0.021 | 5 |
| property-matched inactives | interaction recall | 0.506 | 0.507 | +0.001 | 0.010 | 5 |
| property-matched inactives | AutoDock-GPU per heavy atom | 0.493 | 0.504 | +0.011 | 0.020 | 5 |
| DUD-E-protocol decoys | gnina (CNN score) | 0.562 | 0.568 | +0.006 | 0.010 | 5 |
| DUD-E-protocol decoys | gnina | 0.546 | 0.557 | +0.011 | 0.010 | 5 |
| DUD-E-protocol decoys | fusion: reciprocal rank | 0.535 | 0.556 | +0.021 | 0.025 | 5 |
| DUD-E-protocol decoys | fusion: redocking-weighted z | 0.529 | 0.547 | +0.018 | 0.022 | 5 |
| DUD-E-protocol decoys | consensus | 0.526 | 0.546 | +0.020 | 0.024 | 5 |
| DUD-E-protocol decoys | fusion: rank mean | 0.523 | 0.544 | +0.021 | 0.026 | 5 |
| DUD-E-protocol decoys | gnina (CNN affinity) | 0.514 | 0.527 | +0.013 | 0.015 | 5 |
| DUD-E-protocol decoys | gnina (Vina-type) per heavy atom | 0.526 | 0.527 | +0.001 | 0.024 | 5 |
| DUD-E-protocol decoys | gnina (Vina-type) per heavy atom^(1/3) | 0.516 | 0.519 | +0.003 | 0.011 | 5 |
| DUD-E-protocol decoys | gnina (Vina-type) | 0.507 | 0.518 | +0.011 | 0.015 | 5 |
| DUD-E-protocol decoys | interaction recall | 0.493 | 0.497 | +0.003 | 0.007 | 5 |
| DUD-E-protocol decoys | autodock-gpu | 0.486 | 0.495 | +0.009 | 0.029 | 5 |
| DUD-E-protocol decoys | autodock-gpu (largest cluster) | 0.483 | 0.494 | +0.011 | 0.025 | 5 |
| DUD-E-protocol decoys | AutoDock-GPU per heavy atom | 0.504 | 0.494 | -0.010 | 0.024 | 5 |
| DUD-E-protocol decoys | AutoDock-GPU per heavy atom^(1/3) | 0.486 | 0.494 | +0.007 | 0.020 | 5 |

**Table 5. Docking against the best structure-free baseline**

| target | negatives | score | best baseline | baseline AUC | score AUC | difference | 95% CI | adjusted CI | verdict |
|---|---|---|---|---|---|---|---|---|---|
| ESR1_ant | DUD-E-protocol decoys | autodock-gpu | nn_sim | 0.769 | 0.426 | -0.343 | [-0.438, -0.245] | [-0.471, -0.209] | below |
| ESR1_ant | DUD-E-protocol decoys | consensus | nn_sim | 0.769 | 0.498 | -0.271 | [-0.360, -0.178] | [-0.395, -0.150] | below |
| ESR1_ant | DUD-E-protocol decoys | gnina | nn_sim | 0.769 | 0.567 | -0.202 | [-0.283, -0.118] | [-0.319, -0.093] | below |
| ESR1_ant | property-matched inactives | autodock-gpu | nn_sim | 0.666 | 0.519 | -0.148 | [-0.245, -0.052] | [-0.286, -0.014] | below |
| ESR1_ant | property-matched inactives | consensus | nn_sim | 0.666 | 0.549 | -0.118 | [-0.204, -0.032] | [-0.244, -0.002] | below |
| ESR1_ant | property-matched inactives | gnina | nn_sim | 0.666 | 0.555 | -0.111 | [-0.192, -0.032] | [-0.225, -0.002] | below |
| ESR1_ant | random inactives | autodock-gpu | nn_sim | 0.715 | 0.651 | -0.064 | [-0.160, +0.035] | [-0.199, +0.072] | not resolved |
| ESR1_ant | random inactives | consensus | nn_sim | 0.715 | 0.685 | -0.03 | [-0.114, +0.053] | [-0.148, +0.086] | not resolved |
| ESR1_ant | random inactives | gnina | nn_sim | 0.715 | 0.637 | -0.077 | [-0.160, +0.004] | [-0.197, +0.036] | not resolved |
| MAPK1 | DUD-E-protocol decoys | autodock-gpu | nn_sim | 0.77 | 0.505 | -0.266 | [-0.311, -0.220] | [-0.328, -0.203] | below |
| MAPK1 | DUD-E-protocol decoys | consensus | nn_sim | 0.77 | 0.544 | -0.227 | [-0.271, -0.180] | [-0.289, -0.163] | below |
| MAPK1 | DUD-E-protocol decoys | gnina | nn_sim | 0.77 | 0.56 | -0.211 | [-0.260, -0.162] | [-0.279, -0.140] | below |
| MAPK1 | property-matched inactives | autodock-gpu | nn_sim | 0.637 | 0.507 | -0.129 | [-0.173, -0.083] | [-0.191, -0.067] | below |
| MAPK1 | property-matched inactives | consensus | nn_sim | 0.637 | 0.541 | -0.096 | [-0.140, -0.051] | [-0.158, -0.031] | below |
| MAPK1 | property-matched inactives | gnina | nn_sim | 0.637 | 0.552 | -0.085 | [-0.132, -0.038] | [-0.152, -0.018] | below |
| MAPK1 | random inactives | autodock-gpu | clogp | 0.703 | 0.572 | -0.131 | [-0.183, -0.078] | [-0.206, -0.062] | below |
| MAPK1 | random inactives | consensus | clogp | 0.703 | 0.592 | -0.111 | [-0.166, -0.057] | [-0.191, -0.039] | below |
| MAPK1 | random inactives | gnina | clogp | 0.703 | 0.571 | -0.132 | [-0.192, -0.071] | [-0.218, -0.048] | below |
| MTORC1 | DUD-E-protocol decoys | autodock-gpu | nn_sim | 0.854 | 0.586 | -0.268 | [-0.341, -0.193] | [-0.368, -0.166] | below |
| MTORC1 | DUD-E-protocol decoys | consensus | nn_sim | 0.854 | 0.601 | -0.253 | [-0.326, -0.179] | [-0.354, -0.149] | below |
| MTORC1 | DUD-E-protocol decoys | gnina | nn_sim | 0.854 | 0.548 | -0.305 | [-0.385, -0.225] | [-0.416, -0.190] | below |
| MTORC1 | property-matched inactives | autodock-gpu | nn_sim | 0.677 | 0.496 | -0.181 | [-0.264, -0.097] | [-0.298, -0.066] | below |
| MTORC1 | property-matched inactives | consensus | nn_sim | 0.677 | 0.481 | -0.196 | [-0.278, -0.113] | [-0.314, -0.085] | below |
| MTORC1 | property-matched inactives | gnina | nn_sim | 0.677 | 0.481 | -0.196 | [-0.285, -0.108] | [-0.322, -0.070] | below |
| MTORC1 | random inactives | autodock-gpu | nn_sim | 0.736 | 0.495 | -0.242 | [-0.332, -0.153] | [-0.369, -0.114] | below |
| MTORC1 | random inactives | consensus | nn_sim | 0.736 | 0.495 | -0.241 | [-0.333, -0.151] | [-0.367, -0.117] | below |
| MTORC1 | random inactives | gnina | nn_sim | 0.736 | 0.491 | -0.245 | [-0.342, -0.148] | [-0.378, -0.108] | below |
| PPARG | DUD-E-protocol decoys | autodock-gpu | nn_sim | 0.682 | 0.491 | -0.191 | [-0.316, -0.063] | [-0.363, -0.012] | below |
| PPARG | DUD-E-protocol decoys | consensus | nn_sim | 0.682 | 0.527 | -0.155 | [-0.305, -0.003] | [-0.368, +0.048] | not resolved |
| PPARG | DUD-E-protocol decoys | gnina | nn_sim | 0.682 | 0.558 | -0.124 | [-0.303, +0.059] | [-0.379, +0.113] | not resolved |
| PPARG | property-matched inactives | autodock-gpu | nn_sim | 0.656 | 0.581 | -0.075 | [-0.194, +0.048] | [-0.236, +0.087] | not resolved |
| PPARG | property-matched inactives | consensus | nn_sim | 0.656 | 0.599 | -0.057 | [-0.195, +0.083] | [-0.243, +0.133] | not resolved |
| PPARG | property-matched inactives | gnina | nn_sim | 0.656 | 0.624 | -0.032 | [-0.195, +0.139] | [-0.258, +0.196] | not resolved |
| PPARG | random inactives | autodock-gpu | n_heavy | 0.767 | 0.773 | 0.005 | [-0.055, +0.062] | [-0.085, +0.083] | not resolved |
| PPARG | random inactives | consensus | n_heavy | 0.767 | 0.81 | 0.042 | [-0.032, +0.110] | [-0.067, +0.134] | not resolved |
| PPARG | random inactives | gnina | n_heavy | 0.767 | 0.783 | 0.016 | [-0.094, +0.118] | [-0.143, +0.152] | not resolved |
| TP53 | DUD-E-protocol decoys | autodock-gpu | nn_sim | 0.769 | 0.421 | -0.348 | [-0.450, -0.244] | [-0.488, -0.210] | below |
| TP53 | DUD-E-protocol decoys | consensus | nn_sim | 0.769 | 0.461 | -0.309 | [-0.397, -0.218] | [-0.434, -0.185] | below |
| TP53 | DUD-E-protocol decoys | gnina | nn_sim | 0.769 | 0.498 | -0.271 | [-0.356, -0.185] | [-0.393, -0.151] | below |
| TP53 | property-matched inactives | autodock-gpu | nn_sim | 0.723 | 0.56 | -0.163 | [-0.260, -0.067] | [-0.300, -0.030] | below |
| TP53 | property-matched inactives | consensus | nn_sim | 0.723 | 0.588 | -0.135 | [-0.224, -0.049] | [-0.260, -0.010] | below |
| TP53 | property-matched inactives | gnina | nn_sim | 0.723 | 0.566 | -0.157 | [-0.241, -0.071] | [-0.275, -0.036] | below |
| TP53 | random inactives | autodock-gpu | nn_sim | 0.734 | 0.589 | -0.145 | [-0.246, -0.048] | [-0.286, -0.002] | below |
| TP53 | random inactives | consensus | nn_sim | 0.734 | 0.632 | -0.102 | [-0.195, -0.010] | [-0.233, +0.031] | not resolved |
| TP53 | random inactives | gnina | nn_sim | 0.734 | 0.573 | -0.161 | [-0.252, -0.067] | [-0.292, -0.031] | below |
