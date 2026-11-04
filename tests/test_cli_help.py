"""Keep the command-line help complete and its examples working."""

import argparse
import shlex

import pandas as pd
import pytest

from hapapp_utils.cli import ASSIGN_STAGES, BUILD_DB_STAGES, EXIT_OK, EXIT_USAGE, _build_parser, main
from hapapp_utils.config import BUILD_DB_SCHEMA, BUILD_DB_TEMPLATE, SCHEMA, TEMPLATE, load_config, workflow_of

COMMANDS = ["init", "check", "build-db", "assign", "hap-pa", "merge-pa"]


def subcommands():
    parser = _build_parser()
    action = next(item for item in parser._actions if isinstance(item, argparse._SubParsersAction))
    return action.choices


def examples(command):
    prefix = f"hapapp-utils {command} "
    lines = [line.strip() for line in subcommands()[command].format_help().splitlines()]
    return [shlex.split(line)[2:] for line in lines if line.startswith(prefix)]


def test_main_help_lists_the_commands_and_the_way_forward():
    text = _build_parser().format_help()
    for command in COMMANDS:
        assert f"\n    {command} " in text, f"{command} is not in the command list"
    assert "typical workflow:" in text
    assert "hapapp-utils <command> --help" in text
    assert "https://hapapp.readthedocs.io" in text


def test_no_command_prints_the_help_and_exits_as_a_usage_error(capsys):
    assert main([]) == EXIT_USAGE
    assert "typical workflow:" in capsys.readouterr().err


@pytest.mark.parametrize("command", COMMANDS)
def test_each_command_says_what_it_does(command):
    parser = subcommands()[command]
    assert parser.description and parser.description in parser.format_help()


@pytest.mark.parametrize("command", COMMANDS)
def test_each_option_has_help(command):
    for action in subcommands()[command]._actions:
        assert action.help, f"{command} {action.option_strings or action.dest} has no help"


@pytest.mark.parametrize("command,stages", [("assign", ASSIGN_STAGES), ("build-db", BUILD_DB_STAGES)])
def test_stages_are_explained_and_accepted(command, stages):
    text = subcommands()[command].format_help()
    listing = text.split("stages:\n")[1].split("\n\n")[0].splitlines()
    assert [line.split()[0] for line in listing] == list(stages)
    assert all(len(line.split()) > 1 for line in listing)
    for stage in stages:
        args = _build_parser().parse_args([command, "-c", "settings.yaml", "--until", stage])
        assert args.until == stage


def test_hap_pa_help_names_the_files_it_writes():
    text = subcommands()["hap-pa"].format_help()
    assert "_hapStatus.csv" in text and "_hapPA.csv" in text


@pytest.mark.parametrize("command", [""] + COMMANDS)
def test_help_fits_an_80_column_terminal(command, monkeypatch):
    monkeypatch.setenv("COLUMNS", "80")
    parser = subcommands()[command] if command else _build_parser()
    too_long = [line for line in parser.format_help().splitlines() if len(line) > 80]
    assert not too_long


@pytest.mark.parametrize("command", COMMANDS)
def test_each_command_shows_an_example(command):
    assert examples(command)


@pytest.mark.parametrize("command,example", [(command, example) for command in COMMANDS[1:4]
                                             for example in examples(command)])
def test_examples_for_settings_files_load(command, example, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "run.yaml").write_text(TEMPLATE)
    (tmp_path / "build.yaml").write_text(BUILD_DB_TEMPLATE)
    args = _build_parser().parse_args([command] + example)
    workflow = workflow_of(args.config)
    if command != "check":
        assert workflow == command, f"this example gives {command} a settings file for {workflow}"
    schema = BUILD_DB_SCHEMA if workflow == "build-db" else SCHEMA
    load_config(args.config, args.overrides, args.outdir, schema)


@pytest.mark.parametrize("example", examples("init"))
def test_init_examples_run(example, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert main(["init"] + example) == EXIT_OK
    assert (tmp_path / example[0]).is_file()


def test_table_examples_run(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    example = examples("hap-pa")[0]
    # The example asks for 20 reads, so only the second allele counts as present
    pd.DataFrame({"AlleleID": ["Chr01_000000100|Ref_0001", "Chr01_000000100|Alt_0002"],
                  "CloneID": ["Chr01_000000100"] * 2, "AlleleSequence": ["ACGT", "ACGA"],
                  "S1": [15, 30], "S2": [1, 1]}).to_csv(example[1], index=False)
    assert main(["hap-pa"] + example) == EXIT_OK
    written = pd.read_csv(tmp_path / "report_hapPA.csv")
    assert written.columns.tolist() == ["AlleleID", "PROJECT01"]
    assert written["AlleleID"].tolist() == ["Chr01_000000100|Alt_0002"]

    example = examples("merge-pa")[0]
    for name, column in zip(example[:2], ("A", "B")):
        (tmp_path / name).write_text(f"AlleleID,{column}\nChr01_000000100|Ref_0001,1\n")
    assert main(["merge-pa"] + example) == EXIT_OK
    assert (tmp_path / example[3]).is_file()
