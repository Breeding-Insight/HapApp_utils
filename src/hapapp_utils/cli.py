"""The hapapp-utils command."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .config import BUILD_DB_SCHEMA, BUILD_DB_TEMPLATE, TEMPLATE, ConfigError, load_config, workflow_of

EXIT_OK, EXIT_ERROR, EXIT_DUPLICATES, EXIT_INTERRUPTED = 0, 1, 3, 130
ASSIGN_STAGES = ("update-ids", "filter", "trim", "blast", "assign-ids", "db-check", "hap-pa")
BUILD_DB_STAGES = ("lookup", "flanks", "orient", "alleles", "duplicates", "finish")


def _build_parser():
    parser = argparse.ArgumentParser(
        prog="hapapp-utils",
        description="Build and maintain microhaplotype databases from targeted genotyping reports.")
    parser.add_argument("--version", action="version", version=f"hapapp-utils {__version__}")
    parser.add_argument("--debug", action="store_true", help="show the full Python error if something unexpected fails")
    commands = parser.add_subparsers(dest="command", required=True)

    init = commands.add_parser("init", help="write a commented settings file to fill in")
    init.add_argument("path", nargs="?", default="run.yaml", help="where to write it (default: run.yaml)")
    init.add_argument("--workflow", choices=("assign", "build-db"), default="assign",
                      help="which command the settings are for (default: assign)")
    init.add_argument("--force", action="store_true", help="replace the file if it already exists")

    def add_run_options(sub):
        sub.add_argument("-c", "--config", required=True, help="settings file (YAML)")
        sub.add_argument("--set", dest="overrides", action="append", default=[], metavar="KEY=VALUE",
                         help="override a setting, e.g. --set filter.min_reads=30 (can be repeated)")
        sub.add_argument("-o", "--outdir", help="output folder, replacing output.dir")

    check = commands.add_parser("check", help="check tools, inputs and settings without running anything")
    add_run_options(check)

    build_db = commands.add_parser("build-db", help="build a panel's first allele database from the reference genome")
    add_run_options(build_db)
    build_db.add_argument("--until", choices=BUILD_DB_STAGES, help="stop after this stage")
    build_db.add_argument("-v", "--verbose", action="store_true", help="also show each step's output, not only in the log")

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
    path.write_text(BUILD_DB_TEMPLATE if args.workflow == "build-db" else TEMPLATE)
    print(f"Wrote {path}. Fill it in, then run: hapapp-utils check -c {path}")
    return EXIT_OK


def _cmd_check(args):
    from .checks import OK, PROBLEM, WARNING, preflight, preflight_build_db

    if not Path(args.config).is_file():
        raise ConfigError(f"Settings file not found: {args.config}")
    if workflow_of(args.config) == "build-db":
        findings = preflight_build_db(load_config(args.config, args.overrides, args.outdir, BUILD_DB_SCHEMA))
    else:
        findings, _ = preflight(load_config(args.config, args.overrides, args.outdir))
    labels = {OK: "  ok     ", WARNING: "  WARNING", PROBLEM: "  PROBLEM"}
    for finding in findings:
        print(f"{labels[finding.level]}  {finding.message}")
    problems = sum(finding.level == PROBLEM for finding in findings)
    warnings = sum(finding.level == WARNING for finding in findings)
    print(f"\n{problems} problem(s), {warnings} warning(s).")
    return EXIT_ERROR if problems else EXIT_OK


def _require_workflow(path, wanted):
    """Stop early with a clear message if a settings file is for the other command."""
    if Path(path).is_file() and workflow_of(path) != wanted:
        other = "assign" if wanted == "build-db" else "build-db"
        raise ConfigError(f"{path} holds settings for {other}, not {wanted}. Run hapapp-utils {other} with it, "
                          f"or start a settings file with: hapapp-utils init --workflow {wanted}")


def _print_run_record(summary):
    print(f"  Full log               {summary['process_readme']}")
    print("  Settings and versions  resolved_config.yaml, run_manifest.json")
    if summary["warnings"]:
        print(f"\n{len(summary['warnings'])} warning(s), shown above and saved in run_manifest.json.")


def _cmd_build_db(args):
    from .workflows.build_db import run_build_db

    _require_workflow(args.config, "build-db")
    cfg = load_config(args.config, args.overrides, args.outdir, BUILD_DB_SCHEMA)
    summary = run_build_db(cfg, until=args.until, verbose=args.verbose, command=" ".join(sys.argv))
    print(f"\n{summary['status'].capitalize()}. Results are in {cfg.output.dir}")
    counts = summary["counts"]
    outputs = summary["outputs"]
    if "alleles" in counts:
        print(f"  Allele database        {outputs['allele_db']}")
        print(f"                         {counts['alleles']} alleles from {counts['markers_in_database']} markers")
        print(f"  Count table            {outputs['match_counts']}")
        print(f"  Lookup table           {outputs['marker_lut']}")
        print(f"  Minus-strand markers   {outputs['bottom_strand_markers']}")
        removed = counts["markers_removed_as_duplicates"] + counts["markers_removed_without_alleles"]
        if removed:
            print(f"  Markers removed        {outputs['markers_removed']}")
            print(f"                         {', '.join(removed[:5])}{' ...' if len(removed) > 5 else ''}")
    _print_run_record(summary)
    if "alleles" in counts:
        print("\nTo use it, put the allele database, count table and lookup table in your assign settings.")
    return EXIT_OK


def _cmd_assign(args):
    from .workflows.assign import run_assign

    _require_workflow(args.config, "assign")
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
    _print_run_record(summary)
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
    handlers = {"init": _cmd_init, "check": _cmd_check, "build-db": _cmd_build_db, "assign": _cmd_assign,
                "hap-pa": _cmd_hap_pa, "merge-pa": _cmd_merge_pa}
    try:
        return handlers[args.command](args)
    except ConfigError as err:
        print(f"[ERROR] {err}", file=sys.stderr)
        return EXIT_ERROR
    except Exception as err:
        from .tools import ToolError
        from .workflows.assign import DuplicateAllelesFound
        from .workflows.common import CannotStart, PipelineError

        if isinstance(err, DuplicateAllelesFound):
            print(f"[STOPPED] {err}", file=sys.stderr)
            return EXIT_DUPLICATES
        if isinstance(err, CannotStart):
            print(f"[ERROR] {err}", file=sys.stderr)
            return EXIT_ERROR
        if isinstance(err, (PipelineError, ToolError)):
            print(f"[ERROR] {err}\nSee the .readme file in the output folder's versions folder for details.",
                  file=sys.stderr)
            return EXIT_ERROR
        if args.debug:
            raise
        print(f"[ERROR] Something unexpected went wrong: {type(err).__name__}: {err}\n"
              "Run again with hapapp-utils --debug to see the full Python error.", file=sys.stderr)
        return EXIT_ERROR
    except KeyboardInterrupt:
        print("\n[STOPPED] Interrupted. The output folder holds a partial run.", file=sys.stderr)
        return EXIT_INTERRUPTED
