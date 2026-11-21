"""Command line entry point.

    vsbench doctor                                  # is this machine set up?
    vsbench check-data --data-root DIR              # is this data usable, before docking?
    vsbench run        --config cfg.json [--stage ...] [--targets ...]
    vsbench report     --config cfg.json

``doctor`` and ``check-data`` need no configuration and no GPU; everything else reads a run
configuration. See README for the full set.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .pipeline import Config, run_all

STUDIES = (("ablation", "quality-gate ablation"),
           ("grid", "docking-box strategies"),
           ("perf", "CPU against GPU performance"),
           ("rescore", "template difficulty, size-controlled AUC and score variants"),
           ("interactions", "interaction recovery for every screened compound"),
           ("qc", "structure quality of every candidate template"),
           ("sites", "pocket detection and druggability of every candidate template"),
           ("metals", "metal-ion audit of every template against its deposition"))

#: Studies that describe the receptors rather than the screen, so they may be asked to cover
#: every target in the data root and not only the ones this config screens.
RECEPTOR_STUDIES = ("qc", "sites", "metals")


def _print_network_report() -> None:
    """State exactly what would leave this machine, and whether it has to.

    Everything that does scientific work runs locally. Two things can reach the network, both
    optional and both replaceable, so a reader can see at a glance whether a run is offline.
    """
    from . import tools
    prank, java = tools.find_p2rank(), tools.find_java()
    local_ranker = "available" if (prank and java) else "NOT available (needs P2Rank and Java)"
    print("\nNetwork")
    print("  Everything that docks, scores or analyses runs locally. Two optional steps can "
          "reach out:")
    print("    DoGSiteScorer   https://proteins.plus   pocket ranking; sends a PDB accession, "
          "not your structure")
    print(f"                    local alternative: P2Rank, {local_ranker}")
    print("    RCSB            https://files.rcsb.org  deposited structures for the metal audit "
          "and QC, cached on disk")
    print('  Set "offline": true in the config\'s receptors block to forbid both; cached '
          "depositions still work.")


def _doctor() -> int:
    """Report external executables, Python packages and the GPU. Exit 1 if a required piece is missing."""
    from .tools import check_python_packages, check_tools, gpu_report

    print(f"vsbench {__version__}  (python {'.'.join(map(str, sys.version_info[:3]))})\n")
    print("External executables")
    missing = []
    for s in check_tools():
        mark = "ok " if s.ok else ("MISSING" if s.required else "absent ")
        print(f"  [{mark:>7}] {s.name:<24} {s.path if s.ok else s.purpose}")
        if s.required and not s.ok:
            missing.append(s.name)
    print("\nPython packages")
    for name, version, required, purpose in check_python_packages():
        mark = "ok " if version else ("MISSING" if required else "absent ")
        print(f"  [{mark:>7}] {name:<24} {version or purpose}")
        if required and not version:
            missing.append(name)
    print(f"\nGPU\n  {gpu_report()}")
    _print_network_report()
    if missing:
        print(f"\n{len(missing)} required component(s) missing: {', '.join(missing)}")
        print("Set VSBENCH_BIN to the directory holding the docking binaries, or install the "
              "packages with:  pip install 'vsbench[all]'")
        return 1
    print("\nAll required components present.")
    return 0


def _check_data(data_root: Path, targets, n_actives: int, per_active_matched: int,
                n_random: int) -> int:
    """Validate target directories before any docking. Exit 1 if any target is unusable."""
    from .data import validate_data_root

    checks = validate_data_root(data_root, targets, n_actives=n_actives,
                                per_active_matched=per_active_matched, n_random=n_random)
    if not checks:
        print(f"no target directory found under {data_root}")
        return 1
    bad = 0
    for c in checks:
        print(f"{'ok  ' if c.ok else 'FAIL'}  {c.name}: {c.n_actives} actives, "
              f"{c.n_inactives} inactives, {len(c.templates)} template(s)"
              f"{' ' + ', '.join(c.templates) if c.templates else ''}")
        for p in c.problems:
            print(f"        - {p}")
        bad += not c.ok
    if len(checks) < 2:
        print("\nnote: the DUD-E-protocol decoy arm draws from the other targets' inactives, so it "
              "needs at least two target directories. Set decoy_per_active to 0 to disable it.")
    print(f"\n{len(checks) - bad} of {len(checks)} target(s) usable.")
    return 1 if bad else 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="vsbench", description=__doc__.splitlines()[0])
    ap.add_argument("--version", action="version", version=f"vsbench {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)

    sub.add_parser("doctor", help="check external executables, Python packages and the GPU")

    chk = sub.add_parser("check-data", help="validate target directories before docking")
    chk.add_argument("--data-root", type=Path, help="directory of target subdirectories")
    chk.add_argument("--config", type=Path, help="take the data root and targets from this config")
    chk.add_argument("--targets", nargs="*")
    chk.add_argument("--n-actives", type=int, default=40)
    chk.add_argument("--per-active-matched", type=int, default=8)
    chk.add_argument("--n-random", type=int, default=300)

    run = sub.add_parser("run", help="run benchmark stages")
    run.add_argument("--config", required=True, type=Path)
    run.add_argument("--stage", choices=["templates", "screen", "analyze", "all"], default="all")
    run.add_argument("--targets", nargs="*", help="restrict to these targets")
    rep = sub.add_parser("report", help="aggregate, build the tables and render the figures")
    rep.add_argument("--config", required=True, type=Path)
    for name, helptext in STUDIES:
        sp = sub.add_parser(name, help=helptext)
        sp.add_argument("--config", required=True, type=Path)
        sp.add_argument("--targets", nargs="*")
        if name == "rescore":
            sp.add_argument("--n-boot", type=int, default=None,
                            help="bootstrap replicates (default: the config's n_boot)")
            # A multiplicity-adjusted interval reads a far quantile, so it needs more draws than
            # the pointwise one to be anything but a couple of extreme replicates. Only the
            # variant comparison carries one, so only it has to pay for them.
            sp.add_argument("--variant-n-boot", type=int, default=None,
                            help="replicates for the score-variant comparison alone "
                                 "(default: --n-boot; about 5000 resolves a 12-variant interval)")
        if name in RECEPTOR_STUDIES:
            sp.add_argument("--all-targets", action="store_true",
                            help="cover every target directory under the data root, "
                                 "not only the ones this config screens")
    box = sub.add_parser("boxdiag", help="redock one template across box sizes")
    box.add_argument("--config", required=True, type=Path)
    box.add_argument("--target", required=True)
    box.add_argument("--template", required=True)
    box.add_argument("--edges", nargs="*", type=float, default=[16.0, 20.0, 24.0, 30.0])
    args = ap.parse_args(argv)

    if args.cmd == "doctor":
        return _doctor()
    if args.cmd == "check-data":
        if args.config:
            cfg = Config.from_json(args.config)
            root, targets = Path(cfg.data_root), args.targets or cfg.targets
            defaults = (cfg.n_actives, cfg.per_active_matched, cfg.n_random)
        elif args.data_root:
            root, targets = args.data_root, args.targets
            defaults = (args.n_actives, args.per_active_matched, args.n_random)
        else:
            ap.error("check-data needs --data-root or --config")
        return _check_data(root, targets, *defaults)

    cfg = Config.from_json(args.config)
    if args.cmd == "report":
        from . import pipeline, report, studies
        pipeline.aggregate(cfg)
        report.build_report(Path(cfg.results_dir), cfg.targets)
        studies.make_figures(cfg)
        return 0
    if args.cmd == "boxdiag":
        from . import studies
        if args.target not in cfg.targets:
            raise SystemExit(f"unknown target {args.target!r}; configured: {', '.join(cfg.targets)}")
        studies.box_diagnostic(cfg, args.target, args.template, tuple(args.edges))
        return 0
    if args.cmd in RECEPTOR_STUDIES:
        from . import receptors
        targets = args.targets or (receptors.all_targets(cfg) if args.all_targets else cfg.targets)
        {"qc": receptors.qc_study,
         "sites": receptors.sites_study,
         "metals": receptors.metal_study}[args.cmd](cfg, targets)
        return 0
    if args.cmd in {name for name, _ in STUDIES}:
        from . import studies
        targets = args.targets or cfg.targets
        if args.cmd == "ablation":
            studies.ablation_study(cfg, targets)
        elif args.cmd == "grid":
            studies.grid_study(cfg, targets)
        elif args.cmd == "rescore":
            studies.rescore_study(cfg, targets, n_boot=args.n_boot,
                                  variant_n_boot=args.variant_n_boot)
        elif args.cmd == "interactions":
            studies.interaction_full_study(cfg, targets)
        else:
            studies.perf_study(cfg, targets[0])
        return 0
    if args.targets:
        unknown = set(args.targets) - set(cfg.targets)
        if unknown:
            ap.error(f"targets not in config: {sorted(unknown)}")
        cfg.targets = args.targets
    stages = ("templates", "screen", "analyze") if args.stage == "all" else (args.stage,)
    failed = run_all(cfg, stages)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
