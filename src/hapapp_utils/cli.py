"""The hapapp-utils command."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .config import TEMPLATE, ConfigError, load_config

EXIT_OK, EXIT_ERROR, EXIT_DUPLICATES = 0, 1, 3
ASSIGN_STAGES = ("update-ids", "filter", "trim", "blast", "assign-ids", "db-check", "hap-pa")


def _build_parser():
    parser = argparse.ArgumentParser(
        prog="hapapp-utils",
        description="Build and maintain microhaplotype databases from targeted genotyping reports.")
    parser.add_argument("--version", action="version", version=f"hapapp-utils {__version__}")
    commands = parser.add_subparsers(dest="command", required=True)

    init = commands.add_parser("init", help="write a commented settings file to fill in")
    init.add_argument("path", nargs="?", default="run.yaml", help="where to write it (default: run.yaml)")
    init.add_argument("--force", action="store_true", help="replace the file if it already exists")

    def add_run_options(sub):
        sub.add_argument("-c", "--config", required=True, help="settings file (YAML)")
        sub.add_argument("--set", dest="overrides", action="append", default=[], metavar="KEY=VALUE",
                         help="override a setting, e.g. --set filter.min_reads=30 (can be repeated)")
        sub.add_argument("-o", "--outdir", help="output folder, replacing output.dir")

    check = commands.add_parser("check", help="check tools, inputs and settings without running anything")
    add_run_options(check)

    assign = commands.add_parser("assign", help="turn a MADC into fixed allele IDs and a new database version")
    add_run_options(assign)
    assign.add_argument("--until", choices=ASSIGN_STAGES, help="stop after this stage")
    assign.add_argument("-v", "--verbose", action="store_true", help="also show each step's output, not only in the log")

    pa = commands.add_parser("hap-pa", help="build presence/absence tables from a report with fixed allele IDs")
    pa.add_argument("--madc", required=True, help="report with fixed allele IDs")
    pa.add_argument("--descriptor", required=True, help="column name for this project in the matrix")
    pa.add_argument("--read-threshold", type=int, default=10, help="reads needed to count as present (default: 10)")

    merge = commands.add_parser("merge-pa", help="merge presence/absence matrices from several projects")
    merge.add_argument("reports", nargs="+", help="matrices to merge")
    merge.add_argument("-o", "--output", required=True, help="merged matrix to write")
    return parser


def _cmd_init(args):
    path = Path(args.path)
    if path.exists() and not args.force:
        print(f"[ERROR] {path} already exists. Use --force to replace it.", file=sys.stderr)
        return EXIT_ERROR
    path.write_text(TEMPLATE)
    print(f"Wrote {path}. Fill it in, then run: hapapp-utils check -c {path}")
    return EXIT_OK


def _cmd_check(args):
    from .checks import OK, PROBLEM, WARNING, preflight

    cfg = load_config(args.config, args.overrides, args.outdir)
    findings, _ = preflight(cfg)
    labels = {OK: "  ok     ", WARNING: "  WARNING", PROBLEM: "  PROBLEM"}
    for finding in findings:
        print(f"{labels[finding.level]}  {finding.message}")
    problems = sum(finding.level == PROBLEM for finding in findings)
    warnings = sum(finding.level == WARNING for finding in findings)
    print(f"\n{problems} problem(s), {warnings} warning(s).")
    return EXIT_ERROR if problems else EXIT_OK


def _cmd_assign(args):
    from .workflows.assign import run_assign

    cfg = load_config(args.config, args.overrides, args.outdir)
    summary = run_assign(cfg, until=args.until, verbose=args.verbose, command=" ".join(sys.argv))
    print(f"\n{summary['status'].capitalize()}. Results are in {cfg.output.dir}")
    counts = summary["counts"]
    outputs = summary["outputs"]
    if "report" in outputs:
        print(f"  Report with fixed IDs  {outputs['report']}")
        if "report_rows" in counts:
            print(f"                         {counts['report_rows']} alleles, {counts['report_samples']} samples")
    if "hap_pa" in outputs:
        print(f"  Allele database        {outputs['allele_db']}")
        print(f"                         {counts['database_after']} sequences, was {counts['database_before']}")
        print(f"  Count table            {outputs['match_counts']}")
        print(f"  Presence/absence       {outputs['hap_pa']}")
        print(f"  Sample status          {outputs['hap_status']}")
    print("  Full log               pipeline.log")
    print("  Settings and versions  resolved_config.yaml, run_manifest.json")
    return EXIT_OK


def _cmd_hap_pa(args):
    from .steps.hap_pa import util_haplotype_presence_absence_matrix as hap_pa

    hap_pa.main(["--madc", args.madc, "--descriptor", args.descriptor, "--read_threshold", str(args.read_threshold)])
    return EXIT_OK


def _cmd_merge_pa(args):
    from .steps.hap_pa import util_hapDB_merge_presence_absence_PIhaps as merge_pa

    merge_pa.concat_madc(args.reports, args.output)
    return EXIT_OK


def main(argv=None):
    args = _build_parser().parse_args(argv)
    handlers = {"init": _cmd_init, "check": _cmd_check, "assign": _cmd_assign,
                "hap-pa": _cmd_hap_pa, "merge-pa": _cmd_merge_pa}
    try:
        return handlers[args.command](args)
    except ConfigError as err:
        print(f"[ERROR] {err}", file=sys.stderr)
        return EXIT_ERROR
    except Exception as err:
        from .tools import ToolError
        from .workflows.assign import DuplicateAllelesFound, PipelineError

        if isinstance(err, DuplicateAllelesFound):
            print(f"[STOPPED] {err}", file=sys.stderr)
            return EXIT_DUPLICATES
        if isinstance(err, (PipelineError, ToolError)):
            print(f"[ERROR] {err}\nSee pipeline.log in the output folder for details.", file=sys.stderr)
            return EXIT_ERROR
        raise
