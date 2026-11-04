"""The step modules on their own, including regressions for bugs found on real data."""

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

import hapapp_utils
from hapapp_utils.steps.assign import (
    step00_rm_dupTag_closeTag_rawMADC_v1 as remove_tags,
    step01_filter_missing_AND_ext_matchAlleles_from_madc_v1 as filter_alleles,
    step03_check_cutadapt_allele_uniqueness_AND_update_tmp_rename_report_v1_1 as update_trimmed,
)

from conftest import MARKERS, SAMPLES, make_dataset

STEP_FILES = sorted(Path(hapapp_utils.__file__).parent.glob("steps/*/[!_]*.py"))


@pytest.mark.parametrize("path", STEP_FILES, ids=lambda path: path.stem)
def test_each_step_still_runs_as_a_plain_script(path):
    result = subprocess.run([sys.executable, str(path), "--help"], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


def test_step_modules_are_all_found():
    assert len(STEP_FILES) == 18


def fasta_sequences(path):
    return [line.strip() for line in open(path) if not line.startswith(">")]


def test_filter_ignores_madc_summary_columns(dataset):
    # The noise allele has reads in one sample. AvgCountRef and AvgCountSnp must not count as two more.
    filter_alleles.main([str(dataset.madc), "--min_samples", "2", "--min_reads", "10", "--split_by_colID", ""])
    kept = fasta_sequences(str(dataset.madc).replace(".csv", "_match.fa"))
    for marker in MARKERS:
        assert dataset.noise[marker] not in kept
        assert dataset.new_alleles[marker]["Other"] in kept


@pytest.mark.parametrize("min_samples,n_samples,needed", [
    (None, 94, 5), (None, 400, 10), (None, 10, 1), (2, 94, 2), (10, 94, 10), (20, 94, 20),
])
def test_five_percent_rule_only_applies_without_a_chosen_min_samples(min_samples, n_samples, needed):
    assert filter_alleles.required_samples(min_samples, n_samples) == needed


def test_filter_uses_a_chosen_min_samples_as_written(dataset):
    # 5% of the 40 samples is 2. Lowering this setting to 2 would keep the AltMatch alleles, which have 3 samples.
    filter_alleles.main([str(dataset.madc), "--min_samples", "4", "--min_reads", "10", "--split_by_colID", ""])
    kept = fasta_sequences(str(dataset.madc).replace(".csv", "_match.fa"))
    for marker in MARKERS:
        assert dataset.new_alleles[marker]["Other"] in kept
        assert dataset.new_alleles[marker]["AltMatch"] not in kept
    log = open(str(dataset.madc).replace(".csv", "_match.log")).read()
    assert "Retain an allele with at least 4 samples, each having >= 10 reads" in log
    assert "whichever is smaller" not in log


def test_filter_uses_the_default_rule_without_min_samples(dataset):
    filter_alleles.main([str(dataset.madc), "--min_reads", "10", "--split_by_colID", ""])
    kept = fasta_sequences(str(dataset.madc).replace(".csv", "_match.fa"))
    for marker in MARKERS:
        assert dataset.new_alleles[marker]["AltMatch"] in kept
        assert dataset.noise[marker] not in kept
    log = open(str(dataset.madc).replace(".csv", "_match.log")).read()
    assert "whichever is smaller" in log and "At least 5% samples (2)" in log


def test_filter_keeps_summary_columns_in_temporary_report(dataset):
    # step05 removes these columns itself, using its read_col setting
    filter_alleles.main([str(dataset.madc), "--min_samples", "2", "--min_reads", "10", "--split_by_colID", ""])
    header = open(str(dataset.madc).replace(".csv", "_tmp_rename.csv")).readline().strip().split(",")
    assert len(header) == 16 + len(SAMPLES)


def test_filter_uses_given_first_sample_column(tmp_path):
    data = make_dataset(tmp_path / "data")
    # Without leading * rows the column has to be supplied
    lines = data.madc.read_text().splitlines(keepends=True)
    data.madc.write_text("".join(line for line in lines if not line.startswith("*")))
    filter_alleles.main([str(data.madc), "--min_samples", "2", "--min_reads", "10",
                         "--first_sample_col", "17", "--split_by_colID", ""])
    kept = fasta_sequences(str(data.madc).replace(".csv", "_match.fa"))
    assert data.noise[MARKERS[0]] not in kept


def write_report(path, rows):
    frame = pd.DataFrame(rows, columns=["AlleleID", "CloneID", "AlleleSequence", "sample1"])
    frame.to_csv(path, index=False)
    return frame


MARKER = "Chr01_000000001"
FULL_MARKER_ROWS = [
    (f"{MARKER}|Ref_0001", MARKER, "AAAA", 10),
    (f"{MARKER}|Alt_0002", MARKER, "AAAC", 8),
    (f"{MARKER}|RefMatch_tmp_0003", MARKER, "AATA", 5),
    (f"{MARKER}|RefMatch_tmp_0004", MARKER, "AATG", 4),
    (f"{MARKER}|AltMatch_tmp_0005", MARKER, "AACC", 3),
    (f"{MARKER}|Other_tmp_0006", MARKER, "GGGG", 2),
]


def test_reorder_keeps_other_alleles(tmp_path):
    # From Breeding-Insight/HapApp_utils#1: Other alleles were dropped, leaving the count table ahead of the FASTA
    report = tmp_path / "report.csv"
    frame = write_report(report, FULL_MARKER_ROWS)
    empty = pd.DataFrame(columns=frame.columns).set_index("AlleleID")
    update_trimmed.update_tmp_rename_report(str(report), [], empty, {}, {})
    out = pd.read_csv(tmp_path / "report_updatedSeq.csv")
    assert [allele.split("|")[1] for allele in out["AlleleID"]] == [
        "Ref_0001", "Alt_0002", "RefMatch_tmp_0003", "RefMatch_tmp_0004", "AltMatch_tmp_0005", "Other_tmp_0006"]


def test_reorder_stops_when_a_marker_has_no_ref_row(tmp_path):
    report = tmp_path / "report.csv"
    write_report(report, FULL_MARKER_ROWS[1:])
    with pytest.raises(ValueError, match="changed its rows"):
        update_trimmed.update_tmp_rename_report(str(report), [], None, {}, {})


@pytest.fixture
def trimmed_inputs(tmp_path):
    folder = tmp_path / "folder with space"  # the shell grep used to fail on this
    folder.mkdir()
    report = folder / "report.csv"
    write_report(report, [
        (f"{MARKER}|Ref_0001", MARKER, "AAAATTTT", 10),
        (f"{MARKER}|Alt_0002", MARKER, "AAACTTTT", 8),
        (f"{MARKER}|RefMatch_tmp_0001", MARKER, "AATATTTT", 5),
        (f"{MARKER}|Other_tmp_0001", MARKER, "GGGGTTTT", 2),
    ])
    db = folder / "db.fa"
    db.write_text(f">{MARKER}|Ref_0001\nAAAA\n>{MARKER}|Alt_0002\nAAAC\n")
    return folder, report, db


def test_trimmed_sequences_are_used_when_there_are_no_duplicates(trimmed_inputs):
    folder, report, db = trimmed_inputs
    fasta = folder / "alleles_cutadapt.fa"
    fasta.write_text(f">{MARKER}|RefMatch_tmp_0001\nAATA\n>{MARKER}|Other_tmp_0001\nGGGG\n")
    update_trimmed.main([str(fasta), str(db), str(report)])
    assert (folder / "alleles_cutadapt_unique.fa").exists()
    out = pd.read_csv(folder / "report_updatedSeq.csv").set_index("AlleleID")
    assert out.loc[f"{MARKER}|RefMatch_tmp_0001", "AlleleSequence"] == "AATA"
    assert out.loc[f"{MARKER}|Other_tmp_0001", "AlleleSequence"] == "GGGG"
    assert out.loc[f"{MARKER}|Ref_0001", "AlleleSequence"] == "AAAA"


def test_duplicates_after_trimming_are_merged(trimmed_inputs):
    folder, report, db = trimmed_inputs
    fasta = folder / "alleles_cutadapt.fa"
    fasta.write_text(f">{MARKER}|RefMatch_tmp_0001\nAATA\n>{MARKER}|Other_tmp_0001\nAATA\n")
    update_trimmed.main([str(fasta), str(db), str(report)])
    out = pd.read_csv(folder / "report_updatedSeq.csv").set_index("AlleleID")
    assert out.loc[f"{MARKER}|RefMatch_tmp_0001", "sample1"] == 7
    assert f"{MARKER}|Other_tmp_0001" not in out.index


def test_remove_tags_drops_only_the_listed_markers(dataset, tmp_path):
    tags = tmp_path / "tags.txt"
    tags.write_text(MARKERS[0] + "\n")
    out = tmp_path / "without.csv"
    remove_tags.main([str(tags), str(dataset.madc), str(out)])
    kept = [line.split(",")[0] for line in out.read_text().splitlines()]
    assert "AlleleID" in kept
    assert not any(allele.startswith(MARKERS[0]) for allele in kept)
    assert any(allele.startswith(MARKERS[1]) for allele in kept)
