"""The hapapp-utils command."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .config import BUILD_DB_SCHEMA, BUILD_DB_TEMPLATE, TEMPLATE, ConfigError, load_config, workflow_of

EXIT_OK, EXIT_ERROR, EXIT_USAGE, EXIT_DUPLICATES, EXIT_INTERRUPTED = 0, 1, 2, 3, 130

ASSIGN_STAGE_HELP = {
    "update-ids": "replace marker names with your standard IDs",
    "filter": "keep alleles with enough reads in enough samples",
    "trim": "remove adapters if reads run past the design length",
    "blast": "align the kept alleles to the database",
    "assign-ids": "reuse or issue IDs, then write the report and new database",
    "db-check": "confirm each new ID has a sequence",
    "hap-pa": "build the presence/absence tables",
}
BUILD_DB_STAGE_HELP = {
    "lookup": "build the marker lookup table from the probe file",
    "flanks": "fetch the sequence around each marker from the genome",
    "orient": "find each marker's strand by aligning the MADC's amplicons",
    "alleles": "cut allele_length bases for each Ref and Alt allele",
    "duplicates": "remove markers with near-identical Ref alleles",
    "finish": "write the database and count table as _v001",
}
ASSIGN_STAGES = tuple(ASSIGN_STAGE_HELP)
BUILD_DB_STAGES = tuple(BUILD_DB_STAGE_HELP)

DEBUG_HELP = "show the full Python error if something unexpected fails"

# argparse keeps the line breaks in these, so tests/test_cli_help.py holds them to 80 columns
MAIN_DESCRIPTION = """\
hapapp-utils gives the alleles in MADC reports fixed IDs that stay the same
from project to project. It keeps a versioned allele database for each panel
and builds presence/absence tables from the results."""
MAIN_EPILOG = """\
typical workflow:
  1. first database for a panel   init --workflow build-db, check, build-db
  2. each project                 init, check, assign
  3. redo tables (optional)       hap-pa, merge-pa

Run 'hapapp-utils <command> --help' for a command's options.
Documentation: https://hapapp.readthedocs.io"""
INIT_DESCRIPTION = """\
Write a settings file for assign or build-db. The file lists every setting.
Optional settings start commented out at their default values."""
CHECK_DESCRIPTION = """\
Check the tools, inputs and settings for a run without running anything. The
settings file's sections show whether it is for assign or build-db. A warning
lets a run go ahead. A problem stops it, and check exits with code 1."""
BUILD_DB_DESCRIPTION = """\
Build a panel's first allele database from the reference genome. The tool takes
the Ref and Alt alleles from the genome, because a MADC's sequences can contain
IUPAC codes. It also writes the lookup table and count table that assign needs."""
ASSIGN_DESCRIPTION = """\
Turn a MADC into a report with fixed allele IDs, a new database version and
presence/absence tables. The run copies its inputs into a new folder and works
there, so your input database stays as it was."""
HAP_PA_DESCRIPTION = """\
Build presence/absence tables from a report with fixed allele IDs. The assign
command does this as its last stage, so use hap-pa to redo the tables with
another read threshold. It writes <report>_hapStatus.csv and <report>_hapPA.csv
next to the report."""
MERGE_PA_DESCRIPTION = """\
Merge presence/absence matrices from several projects into one table. It keeps
every allele, fills gaps with 0 and recounts Presence_count. The assign command
can do this for you if you set inputs.previous_hap_pa."""


def _build_parser():
    parser = argparse.ArgumentParser(prog="hapapp-utils", description=MAIN_DESCRIPTION, epilog=MAIN_EPILOG,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", action="version", version=f"hapapp-utils {__version__}")
    parser.add_argument("--debug", action="store_true", help=DEBUG_HELP)
    # main() prints this help if no command is given, so the command is not required here
    commands = parser.add_subparsers(dest="command", title="commands", metavar="<command>")

    def add_command(name, summary, description, examples, stages=None):
        sections = ["examples:\n" + "\n".join(f"  hapapp-utils {name} {example}" for example in examples)]
        if stages:
            sections.insert(0, "stages:\n" + "\n".join(f"  {stage:<12} {text}" for stage, text in stages.items()))
        return commands.add_parser(name, help=summary, description=description,
                                   epilog="\n\n".join(sections), formatter_class=argparse.RawDescriptionHelpFormatter)

    def add_run_options(sub, config_help):
        sub.add_argument("-c", "--config", required=True, metavar="FILE", help=config_help)
        sub.add_argument("--set", dest="overrides", action="append", default=[], metavar="KEY=VALUE",
                         help="change a setting for this run; use once per setting")
        sub.add_argument("-o", "--outdir", metavar="DIR", help="output folder, in place of output.dir")

    def add_stage_options(sub, stages):
        sub.add_argument("--until", choices=stages, metavar="STAGE", help="stop after this stage")
        sub.add_argument("-v", "--verbose", action="store_true",
                         help="show each step's output on screen as well as in the process readme")

    init = add_command("init", "write a commented settings file to fill in", INIT_DESCRIPTION,
                       ["run.yaml", "build.yaml --workflow build-db"])
    init.add_argument("path", nargs="?", default="run.yaml", help="file to write (default: run.yaml)")
    init.add_argument("--workflow", choices=("assign", "build-db"), default="assign",
                      help="the command these settings are for (default: assign)")
    init.add_argument("--force", action="store_true", help="replace the file if it already exists")

    check = add_command("check", "check tools, inputs and settings without running anything", CHECK_DESCRIPTION,
                        ["-c run.yaml", "-c run.yaml --set filter.min_reads=30",
                         "-c build.yaml --set database.allele_length=120"])
    add_run_options(check, "settings file (YAML) for assign or build-db")

    build_db = add_command("build-db", "build a panel's first allele database from the reference genome",
                           BUILD_DB_DESCRIPTION,
                           ["-c build.yaml", "-c build.yaml --set database.allele_length=120",
                            "-c build.yaml --until orient -o runs/orient_only"], BUILD_DB_STAGE_HELP)
    add_run_options(build_db, "settings file (YAML); 'hapapp-utils init --workflow build-db' writes one")
    add_stage_options(build_db, BUILD_DB_STAGES)

    assign = add_command("assign", "turn a MADC into fixed allele IDs and a new database version", ASSIGN_DESCRIPTION,
                         ["-c run.yaml", "-c run.yaml --set filter.min_reads=30 -o runs/trial",
                          "-c run.yaml --until filter -o runs/filter_only"], ASSIGN_STAGE_HELP)
    add_run_options(assign, "settings file (YAML); 'hapapp-utils init' writes one")
    add_stage_options(assign, ASSIGN_STAGES)

    pa = add_command("hap-pa", "build presence/absence tables from a report with fixed IDs", HAP_PA_DESCRIPTION,
                     ["--madc report.csv --descriptor PROJECT01 --min-reads 20"])
    pa.add_argument("--madc", required=True, metavar="FILE",
                    help="report with fixed allele IDs, without the Code_version column")
    pa.add_argument("--descriptor", required=True, metavar="NAME", help="column name for this project in the matrix")
    pa.add_argument("--min-reads", type=int, default=10, metavar="N", help="reads needed to count as present (default: 10)")

    merge = add_command("merge-pa", "merge presence/absence matrices from several projects", MERGE_PA_DESCRIPTION,
                        ["earlier_hapPA.csv this_hapPA.csv -o all_hapPA.csv"])
    merge.add_argument("matrices", nargs="+", metavar="MATRIX", help="_hapPA.csv files or earlier merged matrices")
    merge.add_argument("-o", "--output", required=True, metavar="FILE", help="merged matrix to write")

    # The same option after the command. SUPPRESS stops the command's default from undoing a --debug given before it.
    for sub in commands.choices.values():
        sub.add_argument("--debug", action="store_true", default=argparse.SUPPRESS, help=DEBUG_HELP)
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

    hap_pa.main(["--madc", args.madc, "--descriptor", args.descriptor, "--read_threshold", str(args.min_reads)])
    return EXIT_OK


def _cmd_merge_pa(args):
    from .steps.hap_pa import util_hapDB_merge_presence_absence_PIhaps as merge_pa

    merge_pa.concat_madc(args.matrices, args.output)
    return EXIT_OK


def main(argv=None):
    parser = _build_parser()
    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help(sys.stderr)
        return EXIT_USAGE
    handlers ={"init": _cmd_init, "check": _cmd_check, "build-db": _cmd_build_db, "assign": _cmd_assign,
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
