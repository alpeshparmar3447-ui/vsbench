# Vendored: AutoDock4Zn

Copied unmodified from this project's docking module (itself vendored from
https://github.com/ccsb-scripps/AutoDock-Vina, `develop` branch).
Paper: Santos-Martins et al., J. Chem. Inf. Model. 2014, 54, 2371-2379.

* `AD4Zn.dat` -- AutoDock4 force-field parameters extended with `Zn` and `TZ`
  (tetrahedral zinc pseudo-atom). GNU GPLv2 per its own header.
* `zinc_pseudo.py` -- adds `TZ` pseudo-atoms around each Zn in a receptor PDBQT and
  zeroes the Zn partial charge. Apache-2.0 per the parent repository.

Known upstream behaviour handled in `vsbench.engines`: `zinc_pseudo.py` writes every
AutoDock atom type in upper case, which turns `Zn`, `Mg`, `Ca`, ... into `ZN`, `MG`,
`CA` and makes autogrid4 reject the receptor; the original casing is restored after
the script runs.
