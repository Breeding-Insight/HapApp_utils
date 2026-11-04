import inspect

import pytest
import yaml

from hapapp_utils import __version__
from hapapp_utils.checks import (
    OK,
    PROBLEM,
    WARNING,
    counted_ids_missing_from_fasta,
    minted_ids_missing_from_fasta,
    preflight,
)
from hapapp_utils.cli import ASSIGN_STAGES, main
from hapapp_utils.config import (
    BUILD_DB_SCHEMA,
    BUILD_DB_TEMPLATE,
    SCHEMA,
    TEMPLATE,
    ConfigError,
    config_to_dict,
    load_config,
)
from hapapp_utils.madc import read_layout, read_marker_ids
from hapapp_utils.steps.assign import (
    step01_filter_missing_AND_ext_matchAlleles_from_madc_v1 as filter_alleles,
    step05_parse_madc_matchAlleles_blastn as assign_ids,
)
from hapapp_utils.steps.hap_pa import util_haplotype_presence_absence_matrix as hap_pa
from hapapp_utils.workflows.assign import STAGES

from conftest import MARKERS, SAMPLES, make_dataset


def test_template_is_valid_and_covers_every_setting(tmp_path):
    path = tmp_path / "run.yaml"
    path.write_text(TEMPLATE)
    cfg = load_config(path)
    assert cfg.filter.min_samples == 10
    assert cfg.filter.min_reads == 2
    assert cfg.blast.min_coverage == 90
    assert cfg.hap_pa.min_reads == 10
    for key in SCHEMA:
        assert key.split(".")[-1] in TEMPLATE, key


def test_relative_paths_resolve_from_the_settings_file(dataset):
    cfg = load_config(dataset.config)
    assert cfg.inputs.madc == dataset.madc.resolve()


def test_override_beats_the_file_and_outdir_replaces_output_dir(dataset, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    cfg = load_config(dataset.config, ["filter.min_reads=30", "blast.min_coverage=90"], "elsewhere")
    assert cfg.filter.min_reads == 30
    assert cfg.blast.min_coverage == 90
    assert cfg.output.dir == (tmp_path / "elsewhere").resolve()


def edited(dataset, section, key, value):
    settings = dataset.settings()
    settings[section][key] = value
    dataset.write_settings(settings)
    return dataset.config


@pytest.mark.parametrize("section,key,value,message", [
    ("filter", "min_reads", "ten", "whole number"),
    ("filter", "min_samples", 0, "at least 1"),
    ("blast", "min_coverage", 180, "at most 100"),
    ("blast", "allele_type_truth", "maybe", "report or db"),
    ("panel", "colour", "red", "Unknown setting"),
])
def test_bad_settings_are_rejected(dataset, section, key, value, message):
    with pytest.raises(ConfigError, match=message):
        load_config(edited(dataset, section, key, value))


def test_filter_and_blast_use_the_step_script_defaults(dataset):
    settings = dataset.settings()
    del settings["filter"]
    del settings["blast"]
    dataset.write_settings(settings)
    cfg = load_config(dataset.config)
    assert (cfg.filter.min_samples, cfg.filter.min_reads) == (10, 2)
    assert (cfg.blast.min_coverage, cfg.blast.min_identity) == (90, 90)
    assert cfg.blast.allele_type_truth == "report"
    assert cfg.hap_pa.min_reads == 10


def test_a_setting_can_replace_one_default(dataset):
    settings = dataset.settings()
    settings["filter"] = {"min_reads": 30}
    dataset.write_settings(settings)
    cfg = load_config(dataset.config)
    assert (cfg.filter.min_samples, cfg.filter.min_reads) == (10, 30)


def test_defaults_match_the_step_scripts_themselves():
    filter_source = inspect.getsource(filter_alleles.main)
    assert "'--min_samples', type=int, default=10" in filter_source
    assert "'--min_reads', type=int, default=2" in filter_source
    blast_source = inspect.getsource(assign_ids.main)
    assert "'--cov_threshold', type=float, default=90.0" in blast_source
    assert "'--iden_threshold', type=float, default=90.0" in blast_source
    assert "default='report'" in blast_source
    assert "'--read_threshold', type=int, default=10" in inspect.getsource(hap_pa.main)


def test_required_settings_are_still_required(dataset):
    settings = dataset.settings()
    del settings["panel"]
    dataset.write_settings(settings)
    with pytest.raises(ConfigError, match="panel.design_length, panel.sequencing_length"):
        load_config(dataset.config)


def test_words_yaml_reads_as_booleans_are_rejected(dataset):
    settings = dataset.settings()
    settings["project_id"] = False  # what YAML makes of an unquoted NO
    dataset.write_settings(settings)
    with pytest.raises(ConfigError, match="project_id was read as False. Put its value in quotes."):
        load_config(dataset.config)


def test_unknown_override_is_rejected(dataset):
    with pytest.raises(ConfigError, match="Unknown setting"):
        load_config(dataset.config, ["filter.min_read=30"])


def test_settings_round_trip_through_yaml(dataset):
    saved = yaml.safe_load(yaml.safe_dump(config_to_dict(load_config(dataset.config))))
    assert saved["filter"]["min_samples"] == 2
    assert saved["inputs"]["madc"] == str(dataset.madc.resolve())


def test_code_version_follows_the_package_version(dataset):
    assert __version__ == "2.0.0"
    assert load_config(dataset.config).output.code_version == "v2"
    assert load_config(dataset.config, ["output.code_version=v1"]).output.code_version == "v1"


def test_cli_stage_names_match_the_workflow():
    assert ASSIGN_STAGES == tuple(STAGES)


def test_layout_of_a_raw_madc(dataset):
    layout = read_layout(dataset.madc)
    assert layout.first_sample_column == 17
    assert layout.samples == SAMPLES
    assert layout.filler_columns == 16


def test_layout_of_a_three_column_madc(tmp_path):
    data = make_dataset(tmp_path / "data", raw_format=False)
    layout = read_layout(data.madc)
    assert layout.first_sample_column == 4
    assert layout.samples == SAMPLES


def test_layout_needs_a_header(tmp_path):
    path = tmp_path / "bad.csv"
    path.write_text("not,a,madc\n")
    with pytest.raises(ValueError, match="no header row"):
        read_layout(path)


def test_marker_ids_come_from_allele_ids(dataset):
    assert read_marker_ids(dataset.madc) == set(MARKERS)


def write_counts_and_db(folder, old, new, fasta_ids):
    old_path, new_path, fasta_path = folder / "old.txt", folder / "new.txt", folder / "db.fa"
    old_path.write_text("".join(f"{key}\t{value}\t0\n" for key, value in old.items()))
    new_path.write_text("".join(f"{key}\t{value}\t0\n" for key, value in new.items()))
    fasta_path.write_text("".join(f">{allele}\nACGT\n" for allele in fasta_ids))
    return old_path, new_path, fasta_path


def test_new_id_without_a_sequence_is_caught(tmp_path):
    # The state Breeding-Insight/HapApp_utils#1 reported: counter advanced, no sequence written
    paths = write_counts_and_db(tmp_path, {"M|Other": 0}, {"M|Other": 2}, ["M|Other_0001"])
    assert minted_ids_missing_from_fasta(*paths) == ["M|Other_0002"]


def test_older_gaps_are_not_blamed_on_this_run(tmp_path):
    paths = write_counts_and_db(tmp_path, {"M|Other": 3}, {"M|Other": 4}, ["M|Other_0004"])
    assert minted_ids_missing_from_fasta(*paths) == []
    assert counted_ids_missing_from_fasta(paths[1], paths[2]) == ["M|Other_0001", "M|Other_0002", "M|Other_0003"]


def levels_and_messages(cfg):
    findings, _ = preflight(cfg)
    return {finding.level for finding in findings}, " ".join(finding.message for finding in findings)


def test_preflight_passes_on_good_inputs(dataset):
    levels, _ = levels_and_messages(load_config(dataset.config))
    assert PROBLEM not in levels and WARNING not in levels and OK in levels


def test_preflight_stops_a_first_sample_column_inside_the_summary_columns(dataset):
    levels, messages = levels_and_messages(load_config(dataset.config, ["panel.first_sample_column=5"]))
    assert PROBLEM in levels and "are summary columns, not samples" in messages


def test_preflight_warns_when_the_first_sample_column_skips_samples(dataset):
    levels, messages = levels_and_messages(load_config(dataset.config, ["panel.first_sample_column=18"]))
    assert PROBLEM not in levels and WARNING in levels
    assert "the first 1 sample column(s) will be left out" in messages


def test_preflight_warns_about_the_wrong_sequencing_length(dataset):
    levels, messages = levels_and_messages(load_config(dataset.config, ["panel.sequencing_length=150", "panel.design_length=150"]))
    assert WARNING in levels and "the MADC's Ref alleles are 109 bp" in messages


def test_unknown_setting_suggests_the_closest_one(dataset):
    with pytest.raises(ConfigError, match=r"filter.min_read \(did you mean filter.min_reads\?\)"):
        load_config(dataset.config, ["filter.min_read=30"])


@pytest.mark.parametrize("template,schema,bad_setting,listed,not_listed", [
    (TEMPLATE, SCHEMA, "database.name=x", "min_samples", "allele_length"),
    (BUILD_DB_TEMPLATE, BUILD_DB_SCHEMA, "filter.min_reads=30", "allele_length", "min_samples"),
])
def test_unknown_setting_error_gives_a_command_for_the_right_template(
        template, schema, bad_setting, listed, not_listed, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    settings = tmp_path / "settings.yaml"
    settings.write_text(template)
    with pytest.raises(ConfigError) as caught:
        load_config(settings, [bad_setting], schema=schema)
    command = str(caught.value).split("run: hapapp-utils ")[1].split()
    assert main(command) == 0
    written = (tmp_path / "template.yaml").read_text()
    assert listed in written and not_listed not in written


def test_preflight_reports_a_folder_given_as_an_input(dataset):
    levels, messages = levels_and_messages(load_config(dataset.config, [f"inputs.madc={dataset.folder}"]))
    assert PROBLEM in levels and "is a folder, not a file" in messages


def test_preflight_reports_a_missing_input(dataset):
    dataset.madc.unlink()
    levels, messages = levels_and_messages(load_config(dataset.config))
    assert PROBLEM in levels and "inputs.madc not found" in messages


def test_preflight_warns_about_markers_missing_from_lookup(dataset):
    lines = dataset.lookup.read_text().splitlines(keepends=True)
    dataset.lookup.write_text("".join(lines[1:]))
    levels, messages = levels_and_messages(load_config(dataset.config))
    assert WARNING in levels and MARKERS[0] in messages


def test_preflight_warns_about_dup_tags_that_match_nothing(dataset):
    tags = dataset.folder / "tags.txt"
    tags.write_text("old_style_marker_name\n")
    settings = dataset.settings()
    settings["inputs"]["dup_tags"] = tags.name
    dataset.write_settings(settings)
    levels, messages = levels_and_messages(load_config(dataset.config))
    assert WARNING in levels and "old_style_marker_name" in messages


def test_preflight_notes_when_the_five_percent_rule_lowers_min_samples(dataset):
    settings = dataset.settings()
    settings["filter"]["min_samples"] = 10
    dataset.write_settings(settings)
    _, messages = levels_and_messages(load_config(dataset.config))
    assert "lowers it to 2" in messages


def test_preflight_refuses_a_used_output_folder(dataset):
    dataset.outdir.mkdir()
    (dataset.outdir / "old_result.csv").write_text("x")
    levels, messages = levels_and_messages(load_config(dataset.config))
    assert PROBLEM in levels and "not empty" in messages
