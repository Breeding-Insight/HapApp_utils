"""The whole assign workflow through the command-line entry point. These need BLAST installed."""

import hashlib
import json
import os
import shutil
from pathlib import Path

import pandas as pd
import pytest

from hapapp_utils.checks import minted_ids_missing_from_fasta, read_fasta_ids
from hapapp_utils.cli import EXIT_DUPLICATES, EXIT_ERROR, EXIT_OK, main

from conftest import MARKERS, SAMPLES, make_dataset, requires_blast, requires_cutadapt

pytestmark = requires_blast


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fasta_records(path):
    records = {}
    name = None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].strip()
        else:
            records[name] = line.strip()
    return records


def run(data, *extra):
    return main(["assign", "-c", str(data.config), *extra])


def manifest(data):
    return json.loads((data.outdir / "run_manifest.json").read_text())


def final_report(data):
    return pd.read_csv(data.outdir / manifest(data)["outputs"]["report"])


def assert_expected_alleles(frame):
    """Six rows per marker: Ref, Alt, the known RefMatch and one new allele of each class."""
    for marker in MARKERS:
        classes = sorted(frame.loc[frame["CloneID"] == marker, "AlleleID"].str.split("|").str[1])
        assert classes == ["AltMatch_0001", "Alt_0002", "Other_0001", "RefMatch_0001", "RefMatch_0002", "Ref_0001"]


def test_assign_end_to_end(dataset, capsys):
    before = {path: sha256(path) for path in (dataset.madc, dataset.allele_db, dataset.match_counts, dataset.lookup)}
    assert run(dataset) == EXIT_OK

    frame = final_report(dataset)
    assert_expected_alleles(frame)
    assert list(frame.columns) == ["Code_version", "AlleleID", "CloneID", "AlleleSequence"] + SAMPLES
    assert set(frame["Code_version"]) == {"v2"}
    assert manifest(dataset)["outputs"]["report"].endswith("_rename_v2.csv")

    versions = dataset.outdir / "versions"
    database = fasta_records(versions / "test_allele_db_v002.fa")
    assert len(database) == 18
    for marker in MARKERS:
        assert dataset.noise[marker] not in database.values()
        assert database[f"{marker}|Other_0001"] == dataset.new_alleles[marker]["Other"]
    assert minted_ids_missing_from_fasta(dataset.match_counts, versions / "test_allele_db_v002_matchCnt_lut.txt",
                                         versions / "test_allele_db_v002.fa") == []

    # Read counts stay with the right sample
    row = frame[frame["AlleleID"] == f"{MARKERS[0]}|Other_0001"].iloc[0]
    assert [row[sample] for sample in SAMPLES[:5]] == [14, 14, 14, 14, 0]

    # Inputs are only ever read
    for path, digest in before.items():
        assert sha256(path) == digest
    assert not list(dataset.folder.glob("*_snpID*"))
    assert not (dataset.folder / "test_allele_db_v002.fa").exists()

    # A record of the run is saved with the results
    record = manifest(dataset)
    assert record["status"] == "finished"
    assert record["first_sample_column"] == 17
    assert record["counts"] == {"database_before": 9, "database_after": 18, "report_rows": 18,
                                "report_samples": len(SAMPLES)}
    assert record["inputs"]["madc"]["sha256"] == before[dataset.madc]
    assert "blastn" in record["tools"]
    assert (dataset.outdir / "resolved_config.yaml").exists()
    assert record["process_readme"] == "versions/test_allele_db_v002_process.readme"
    assert (dataset.outdir / record["process_readme"]).read_text()

    matrix = pd.read_csv(dataset.outdir / record["outputs"]["hap_pa"])
    assert set(frame["AlleleID"]) == set(matrix["AlleleID"])
    assert list(matrix.columns) == ["AlleleID", "PROJECT"]
    assert "Finished" in capsys.readouterr().out


def test_step_output_goes_to_the_log_unless_verbose(tmp_path, capsys):
    quiet = make_dataset(tmp_path / "quiet")
    assert run(quiet) == EXIT_OK
    assert "[INFO] Analyze BLAST results" not in capsys.readouterr().out
    assert "[INFO] Analyze BLAST results" in (quiet.outdir / manifest(quiet)["process_readme"]).read_text()

    loud = make_dataset(tmp_path / "loud")
    assert run(loud, "--verbose") == EXIT_OK
    assert "[INFO] Analyze BLAST results" in capsys.readouterr().out


def test_three_column_madc_keeps_every_sample(tmp_path):
    # step05 used to skip 16 columns regardless, losing the first 13 samples of a file like this
    data = make_dataset(tmp_path / "data", raw_format=False)
    assert run(data) == EXIT_OK
    frame = final_report(data)
    assert list(frame.columns) == ["Code_version", "AlleleID", "CloneID", "AlleleSequence"] + SAMPLES
    assert manifest(data)["first_sample_column"] == 4
    assert_expected_alleles(frame)


def test_a_used_output_folder_is_refused(dataset, capsys):
    assert run(dataset) == EXIT_OK
    assert run(dataset) == EXIT_ERROR
    assert "not empty" in capsys.readouterr().err


def test_a_trial_run_leaves_earlier_results_alone(dataset, tmp_path):
    assert run(dataset) == EXIT_OK
    first = sha256(dataset.outdir / manifest(dataset)["outputs"]["allele_db"])
    trial = tmp_path / "trial"
    assert run(dataset, "--set", "filter.min_samples=1", "--set", "filter.min_reads=30", "-o", str(trial)) == EXIT_OK
    assert sha256(dataset.outdir / manifest(dataset)["outputs"]["allele_db"]) == first
    # Only the noise allele has 30 or more reads in a sample, besides Ref and Alt
    database = fasta_records(trial / "versions" / "test_allele_db_v002.fa")
    assert dataset.noise[MARKERS[0]] in database.values()
    assert dataset.new_alleles[MARKERS[0]]["Other"] not in database.values()


def test_new_ids_without_sequences_stop_the_run(dataset, monkeypatch, capsys):
    from hapapp_utils.workflows import assign as assign_workflow

    monkeypatch.setattr(assign_workflow, "minted_ids_missing_from_fasta", lambda *paths: ["M|Other_0001"])
    assert run(dataset) == EXIT_ERROR
    assert "no sequence" in capsys.readouterr().err
    assert manifest(dataset)["status"] == "failed"


def test_until_stops_early(dataset, capsys):
    assert run(dataset, "--until", "filter") == EXIT_OK
    record = manifest(dataset)
    assert record["status"] == "stopped after filter"
    assert not (dataset.outdir / "hapSearch").exists()
    assert (dataset.outdir / "report" / "PROJECT_MADC_snpID_match.fa").exists()


def test_duplicate_alleles_stop_with_their_own_exit_code(tmp_path, capsys):
    data = make_dataset(tmp_path / "data", duplicate_allele=True)
    assert run(data) == EXIT_DUPLICATES
    assert "same sequence" in capsys.readouterr().err
    assert manifest(data)["status"] == "failed"
    assert (data.outdir / "report" / "PROJECT_MADC_snpID_dup.csv").exists()


@pytest.fixture
def tagged(dataset):
    tags = dataset.folder / "tags.txt"
    tags.write_text(MARKERS[1] + "\n")
    settings = dataset.settings()
    settings["inputs"]["dup_tags"] = tags.name
    dataset.write_settings(settings)
    return dataset


def test_dup_tags_remove_markers(tagged):
    assert run(tagged) == EXIT_OK
    frame = final_report(tagged)
    assert MARKERS[1] not in set(frame["CloneID"])
    assert {MARKERS[0], MARKERS[2]} == set(frame["CloneID"])
    assert frame.shape[1] == 4 + len(SAMPLES)


def test_duplicate_check_still_works_with_dup_tags(tmp_path):
    # The bash workflow looked for the wrong file here, so a duplicate did not stop it
    data = make_dataset(tmp_path / "data", duplicate_allele=True)
    tags = data.folder / "tags.txt"
    tags.write_text(MARKERS[1] + "\n")
    settings = data.settings()
    settings["inputs"]["dup_tags"] = tags.name
    data.write_settings(settings)
    assert run(data) == EXIT_DUPLICATES
    assert (data.outdir / "report" / "PROJECT_MADC_snpID_rmDup_dup.csv").exists()


def test_merging_into_an_earlier_matrix(dataset, tmp_path):
    previous = tmp_path / "earlier_hapPA.csv"
    previous.write_text(f"AlleleID,OLD\n{MARKERS[0]}|Ref_0001,1\n{MARKERS[0]}|RefMatch_0099,1\n")
    settings = dataset.settings()
    settings["inputs"]["previous_hap_pa"] = str(previous)
    dataset.write_settings(settings)
    assert run(dataset) == EXIT_OK
    matrix = pd.read_csv(dataset.outdir / manifest(dataset)["outputs"]["hap_pa"]).set_index("AlleleID")
    assert set(matrix.columns) == {"CloneID", "OLD", "PROJECT", "Presence_count"}
    assert matrix.loc[f"{MARKERS[0]}|RefMatch_0099", "OLD"] == 1
    assert matrix.loc[f"{MARKERS[0]}|Other_0001", "PROJECT"] == 1
    assert matrix.loc[f"{MARKERS[0]}|Ref_0001", "Presence_count"] == 2


def test_merging_a_matrix_again_keeps_presence_counts_correct(tmp_path):
    # An old Presence_count used to be summed like a sample, so 3 projects reported 5
    (tmp_path / "a.csv").write_text("AlleleID,A\nM1|Ref_0001,1\nM1|Alt_0002,1\n")
    (tmp_path / "b.csv").write_text("AlleleID,B\nM1|Ref_0001,1\nM1|Alt_0002,0\n")
    (tmp_path / "c.csv").write_text("AlleleID,C\nM1|Ref_0001,1\n")
    first = tmp_path / "ab.csv"
    second = tmp_path / "abc.csv"
    assert main(["merge-pa", str(tmp_path / "a.csv"), str(tmp_path / "b.csv"), "-o", str(first)]) == EXIT_OK
    assert main(["merge-pa", str(first), str(tmp_path / "c.csv"), "-o", str(second)]) == EXIT_OK
    matrix = pd.read_csv(second).set_index("AlleleID")
    assert list(matrix.columns) == ["CloneID", "A", "B", "C", "Presence_count"]
    assert matrix.loc["M1|Ref_0001", "Presence_count"] == 3
    assert matrix.loc["M1|Alt_0002", "Presence_count"] == 1


def test_process_readme_lists_inputs_and_settings(dataset, capsys):
    assert run(dataset, "--set", "filter.min_reads=11") == EXIT_OK
    readme = (dataset.outdir / manifest(dataset)["process_readme"]).read_text()
    header = readme[:readme.index("[1/7]")]
    assert "--- Scripts, files, and parameters used ---" in header
    for path in (dataset.madc, dataset.lookup, dataset.allele_db, dataset.match_counts):
        assert str(path.resolve()) in header
    assert "# First sample column: 17" in header
    assert "# Panel design length: 109" in header
    assert "# Sequencing length: 109" in header
    assert "# BLAST filtering thresholds: 80.0% coverage, 85.0% identity" in header
    assert "# Read filter: at least 2 samples with at least 11 reads" in header
    assert "# Process README: versions/test_allele_db_v002_process.readme" in header
    assert "# Tool version: blastn " in header
    # The step output follows in the same file
    assert "[INFO] Analyze BLAST results" in readme
    assert "Full log               versions/test_allele_db_v002_process.readme" in capsys.readouterr().out


def test_failed_run_keeps_its_process_readme(tmp_path):
    data = make_dataset(tmp_path / "data", duplicate_allele=True)
    assert run(data) == EXIT_DUPLICATES
    readme = data.outdir / "versions" / "test_allele_db_v002_process.readme"
    assert "[ERROR] The MADC has the same sequence" in readme.read_text()
    assert manifest(data)["process_readme"] == "versions/test_allele_db_v002_process.readme"


def test_no_new_alleles_writes_no_new_database(tmp_path):
    data = make_dataset(tmp_path / "data")
    # Keep only Ref, Alt and alleles already in the database
    lines = data.madc.read_text().splitlines(keepends=True)
    known = set(fasta_records(data.allele_db).values())
    data.madc.write_text("".join(line for line in lines
                                 if line.startswith(("*", "AlleleID")) or line.split(",")[2] in known))
    assert run(data) == EXIT_OK
    record = manifest(data)
    assert record["outputs"]["allele_db"] == "versions/test_allele_db_v001.fa"
    assert not (data.outdir / "versions" / "test_allele_db_v002.fa").exists()
    assert record["counts"]["report_rows"] == 9
    assert (data.outdir / record["outputs"]["hap_pa"]).name == "test_allele_db_v001_PROJECT_hapPA.csv"
    # As in the bash workflow, the readme is named after the database that was kept
    assert record["process_readme"] == "versions/test_allele_db_v001_noNewAllele.readme"
    assert (data.outdir / record["process_readme"]).exists()
    assert not (data.outdir / "versions" / "test_allele_db_v002_process.readme").exists()


def test_two_alleles_with_one_sequence_in_a_new_database_are_collapsed(dataset):
    # A database that already holds the same sequence twice gets one more version with one copy
    existing = dataset.allele_db.read_text()
    dataset.allele_db.write_text(existing + f">{MARKERS[0]}|RefMatch_0002\n"
                                 + fasta_records(dataset.allele_db)[f"{MARKERS[0]}|RefMatch_0001"] + "\n")
    counts = dataset.match_counts.read_text().replace(f"{MARKERS[0]}|RefMatch\t1", f"{MARKERS[0]}|RefMatch\t2")
    dataset.match_counts.write_text(counts)
    assert run(dataset) == EXIT_OK
    record = manifest(dataset)
    assert record["outputs"]["allele_db"] == "versions/test_allele_db_v003.fa"
    assert record["outputs"]["match_counts"] == "versions/test_allele_db_v003_matchCnt_lut.txt"
    final = fasta_records(dataset.outdir / record["outputs"]["allele_db"])
    assert len(final) == len(set(final.values()))


@requires_cutadapt
@pytest.mark.parametrize("duplicate_after_trim", [False, True], ids=["no_duplicates", "duplicates"])
def test_adapters_are_removed(tmp_path, duplicate_after_trim):
    data = make_dataset(tmp_path / "data", trim_design=81, duplicate_after_trim=duplicate_after_trim)
    assert run(data) == EXIT_OK
    frame = final_report(data)
    assert_expected_alleles(frame)
    assert set(frame["AlleleSequence"].str.len()) == {81}

    versions = data.outdir / "versions"
    database = fasta_records(versions / "test_allele_db_v002.fa")
    assert {len(sequence) for sequence in database.values()} == {81}
    for marker in MARKERS:
        assert database[f"{marker}|Other_0001"] == data.new_alleles[marker]["Other"]
    assert minted_ids_missing_from_fasta(data.match_counts, versions / "test_allele_db_v002_matchCnt_lut.txt",
                                         versions / "test_allele_db_v002.fa") == []
    expected_reads = 26 if duplicate_after_trim else 15
    new_ref = frame[frame["AlleleID"] == f"{MARKERS[0]}|RefMatch_0002"].iloc[0]
    assert new_ref[SAMPLES[0]] == expected_reads


@requires_cutadapt
def test_mixed_insert_sizes_are_trimmed_per_marker(tmp_path):
    inserts = {MARKERS[0]: 81, MARKERS[1]: 109, MARKERS[2]: 81}
    data = make_dataset(tmp_path / "data", trim_design=81, insert_sizes=inserts)
    settings = data.settings()
    settings["panel"]["insert_size_column"] = "insert_size"
    data.write_settings(settings)
    assert run(data) == EXIT_OK
    frame = final_report(data)
    assert_expected_alleles(frame)
    lengths = frame.groupby("CloneID")["AlleleSequence"].apply(lambda column: set(column.str.len()))
    assert lengths[MARKERS[0]] == {81}
    assert lengths[MARKERS[1]] == {109}
    assert lengths[MARKERS[2]] == {81}


def test_check_command(dataset, capsys):
    assert main(["check", "-c", str(dataset.config)]) == EXIT_OK
    out = capsys.readouterr().out
    assert "40 samples, first sample in column 17" in out
    assert "0 problem(s)" in out
    dataset.lookup.unlink()
    assert main(["check", "-c", str(dataset.config)]) == EXIT_ERROR


def test_init_command(tmp_path, capsys):
    target = tmp_path / "new.yaml"
    assert main(["init", str(target)]) == EXIT_OK
    assert "filter:" in target.read_text()
    assert main(["init", str(target)]) == EXIT_ERROR
    assert main(["init", str(target), "--force"]) == EXIT_OK


def test_hap_pa_and_merge_pa_commands(dataset, tmp_path):
    assert run(dataset) == EXIT_OK
    report = tmp_path / "report.csv"
    frame = final_report(dataset).drop(columns="Code_version")
    frame.to_csv(report, index=False)
    assert main(["hap-pa", "--madc", str(report), "--descriptor", "A", "--read-threshold", "20"]) == EXIT_OK
    first = tmp_path / "report_hapPA.csv"
    assert pd.read_csv(first).columns.tolist() == ["AlleleID", "A"]
    second = tmp_path / "other.csv"
    second.write_text(f"AlleleID,B\n{MARKERS[0]}|Ref_0001,1\n")
    merged = tmp_path / "merged.csv"
    assert main(["merge-pa", str(first), str(second), "-o", str(merged)]) == EXIT_OK
    assert set(pd.read_csv(merged).columns) == {"AlleleID", "CloneID", "A", "B", "Presence_count"}


def test_bad_settings_give_a_clean_error(dataset, capsys):
    assert main(["assign", "-c", str(dataset.config), "--set", "filter.min_reads=many"]) == EXIT_ERROR
    assert "whole number" in capsys.readouterr().err


@pytest.mark.skipif(not os.environ.get("HAPAPP_BASELINE_DIR"), reason="HAPAPP_BASELINE_DIR is not set")
def test_real_data_matches_saved_results(tmp_path):
    """Compare a real run with earlier results: HAPAPP_BASELINE_DIR and HAPAPP_BASELINE_CONFIG."""
    baseline = Path(os.environ["HAPAPP_BASELINE_DIR"])
    outdir = tmp_path / "run"
    # Saved results from the bash workflow carry its code version
    assert main(["assign", "-c", os.environ["HAPAPP_BASELINE_CONFIG"], "-o", str(outdir),
                 "--set", "output.code_version=v1"]) == EXIT_OK
    record = json.loads((outdir / "run_manifest.json").read_text())
    for key in ("report", "allele_db", "match_counts", "hap_pa", "hap_status"):
        name = record["outputs"][key]
        assert sha256(outdir / name) == sha256(baseline / name), name
    shutil.rmtree(outdir)
