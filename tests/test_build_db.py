"""The core database workflow on a small synthetic genome, plus its FASTA and probe file readers."""

import json
import random
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import openpyxl
import pandas as pd
import pytest
import yaml

from hapapp_utils.cli import BUILD_DB_STAGES, EXIT_ERROR, EXIT_OK, main
from hapapp_utils.config import BUILD_DB_SCHEMA, BUILD_DB_TEMPLATE, ConfigError, load_config, workflow_of
from hapapp_utils.fasta import fetch_subsequences, iter_fasta, read_lengths
from hapapp_utils.probe_design import _cell_text, read_probe_rows, write_probe_csv
from hapapp_utils.steps.core_db import (
    db00_prep_lut_from_probeDesign as prepare_lookup,
    db02_update_snpID_in_madc_v1 as update_ids,
)
from hapapp_utils.workflows.build_db import STAGES

from conftest import SAMPLES, SUMMARY_COLUMNS, reads_in

requires_tools = pytest.mark.skipif(
    not all(shutil.which(tool) for tool in ("blastn", "mmseqs")), reason="BLAST and MMseqs2 are not installed")
requires_esl_sfetch = pytest.mark.skipif(shutil.which("esl-sfetch") is None, reason="esl-sfetch is not installed")

PROBE_HEADER = ["Marker Name", "Target Sequence", "Reference Genome", "Chrom ", "Chrom Pos", "VariantAllelesDef ",
                "Marker_Type", "Essential Marker ", "Comments"]
AMPLICON = 109
OFFSET = 40  # where the variant sits in the read
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


@dataclass
class BuildDataset:
    folder: Path
    config: Path
    genome: Path
    madc: Path
    outdir: Path
    expected: dict = field(default_factory=dict)
    duplicate_pair: tuple = ()
    placed_badly: tuple = ()

    def settings(self):
        return yaml.safe_load(self.config.read_text())

    def write_settings(self, settings):
        self.config.write_text(yaml.safe_dump(settings))


def make_build_dataset(folder, probe_format="xlsx", use_genome_names=False):
    """A genome, probe file and MADC, with markers of every kind the workflow has to handle."""
    folder.mkdir(parents=True, exist_ok=True)
    rng = random.Random(11)
    chromosomes = {name: [rng.choice("ACGT") for _ in range(6000)] for name in ("Chr01", "Chr02")}
    # The same stretch on both chromosomes makes two markers whose alleles are identical
    chromosomes["Chr02"][3500:4500] = chromosomes["Chr01"][3500:4500]
    chromosomes = {name: "".join(bases) for name, bases in chromosomes.items()}
    genome_names = {"Chr01": "NC_1.1", "Chr02": "NC_2.1"}

    probe_rows = []
    madc_rows = []
    expected = {}

    def add(marker, chrom, pos, variant, marker_type, ref_read, alt_read):
        probe_rows.append([marker, "SEQ", "test_v1", chrom, str(pos), variant, marker_type, "2", "Marker"])
        madc_rows.append((f"{marker}|Ref", marker, ref_read, [50] * len(SAMPLES)))
        madc_rows.append((f"{marker}|Alt", marker, alt_read, [30] * len(SAMPLES)))
        # A new allele for assign to find later, one base from Ref
        bases = list(ref_read)
        bases[80] = "A" if bases[80] != "A" else "C"
        madc_rows.append((f"{marker}|RefMatch", marker, "".join(bases), reads_in(5, 15)))

    def snp(marker, chrom, pos, strand="+", variant=None):
        sequence = chromosomes[chrom]
        ref = sequence[pos - 1]
        alt = {"A": "G", "G": "A", "C": "T", "T": "C"}[ref]
        if strand == "+":
            start = pos - 1 - OFFSET
            ref_read = sequence[start:start + AMPLICON]
            alt_read = ref_read[:OFFSET] + alt + ref_read[OFFSET + 1:]
        else:
            start = pos - AMPLICON + OFFSET
            ref_read = reverse_complement(sequence[start:start + AMPLICON])
            alt_read = reverse_complement(sequence[start:pos - 1] + alt + sequence[pos:start + AMPLICON])
        expected[marker] = (ref_read, alt_read)
        add(marker, chrom, pos, variant or f"[{ref}/{alt}]", "SNP", ref_read, alt_read)
        return ref_read, alt_read

    snp("Chr01_000000500", "Chr01", 500)
    snp("Chr01_000001000", "Chr01", 1000, strand="-")
    snp("Chr02_000000500", "Chr02", 500)

    # IUPAC codes in the MADC are why the database is built from the genome instead
    marker = "Chr01_000001500"
    ref_read, alt_read = snp(marker, "Chr01", 1500)
    madc_rows[-3] = (f"{marker}|Ref", marker, ref_read[:5] + "Y" + ref_read[6:], madc_rows[-3][3])

    # An insertion and a deletion
    marker, pos = "Chr01_000002000", 2000
    sequence = chromosomes["Chr01"]
    anchor = sequence[pos - 1]
    start = pos - 1 - OFFSET
    ref_read = sequence[start:start + AMPLICON]
    alt_read = (sequence[start:pos] + "GC" + sequence[pos:start + AMPLICON])[:AMPLICON]
    expected[marker] = (ref_read, alt_read)
    add(marker, "Chr01", pos, f"[{anchor}/{anchor}GC]", "Indel", ref_read, alt_read)

    marker, pos = "Chr01_000002500", 2500
    deleted = sequence[pos - 1:pos + 2]
    start = pos - 1 - OFFSET
    ref_read = sequence[start:start + AMPLICON]
    alt_read = (sequence[start:pos] + sequence[pos + 2:start + AMPLICON + 2])[:AMPLICON]
    expected[marker] = (ref_read, alt_read)
    add(marker, "Chr01", pos, f"[{deleted}/{deleted[0]}]", "Indel", ref_read, alt_read)

    # Two markers with identical alleles
    pair = ("Chr01_000004000", "Chr02_000004000")
    snp(pair[0], "Chr01", 4000)
    snp(pair[1], "Chr02", 4000)

    # A read that is nowhere near its flanks, and a probe whose bases do not match the genome
    no_hit = "Chr02_000001500"
    junk = "".join(random.Random(99).choice("ACGT") for _ in range(AMPLICON))
    add(no_hit, "Chr02", 1500, f"[{chromosomes['Chr02'][1499]}/N]", "SNP", junk, junk[:OFFSET] + "N" + junk[OFFSET + 1:])
    wrong_base = "Chr02_000002500"
    real = chromosomes["Chr02"][2499]
    others = [base for base in "ACGT" if base != real]
    start = 2500 - 1 - OFFSET
    read = chromosomes["Chr02"][start:start + AMPLICON]
    add(wrong_base, "Chr02", 2500, f"[{others[0]}/{others[1]}]", "SNP", read, read[:OFFSET] + others[1] + read[OFFSET + 1:])

    # In the probe file only
    probe_rows.append(["Chr02_000003000", "SEQ", "test_v1", "Chr02", "3000", "[A/T]", "SNP", "2", "Marker"])

    genome = folder / "genome.fna"
    with open(genome, "w") as out:
        for name, sequence in chromosomes.items():
            title = genome_names[name] if not use_genome_names else name
            out.write(f">{title} some description\n")
            lowered = sequence[:300] + sequence[300:360].lower() + sequence[360:]  # soft-masked stretch
            for start in range(0, len(lowered), 70):
                out.write(lowered[start:start + 70] + "\n")
    names_file = folder / "chromosome_names.csv"
    names_file.write_text("".join(f"{probe},{genome_name}\n" for probe, genome_name in genome_names.items()))

    if probe_format == "xlsx":
        probe = folder / "probe.xlsx"
        book = openpyxl.Workbook()
        sheet = book.active
        sheet.append([None, "Requirements for submitting sequences"])
        sheet.append([None, None])
        sheet.append(["Column 1", "Column 2"])
        sheet.append(PROBE_HEADER)
        for row in probe_rows:
            sheet.append(row[:4] + [int(row[4])] + row[5:])
        book.save(probe)
    else:
        delimiter = "," if probe_format == "csv" else "\t"
        probe = folder / f"probe.{probe_format}"
        probe.write_text("".join(delimiter.join(row) + "\n" for row in [PROBE_HEADER] + probe_rows))

    madc = folder / "PROJECT_MADC.csv"
    filler = ["*"] * (3 + len(SUMMARY_COLUMNS))
    lines = [",".join(filler + ["x"] * len(SAMPLES)) + "\n" for _ in range(7)]
    lines.append(",".join(["AlleleID", "CloneID", "AlleleSequence"] + SUMMARY_COLUMNS + SAMPLES) + "\n")
    for allele_id, marker, sequence, counts in madc_rows:
        summary = [sequence, "0.5", "0.5", "0.5", "0", "0", "0", "0", "0", "0", "25.0", "25.0", "1.0"]
        lines.append(",".join([allele_id, marker, sequence] + summary + [str(value) for value in counts]) + "\n")
    madc.write_text("".join(lines))

    config = folder / "build.yaml"
    settings = {
        "inputs": {"probe_design": probe.name, "madc": madc.name, "genome": genome.name},
        "database": {"name": "core_db", "allele_length": AMPLICON, "flank_length": 150, "threads": 2},
        "output": {"dir": "out"},
    }
    if not use_genome_names:
        settings["inputs"]["chromosome_map"] = names_file.name
    config.write_text(yaml.safe_dump(settings))
    return BuildDataset(folder, config, genome, madc, folder / "out", expected, pair, (no_hit, wrong_base))


@pytest.fixture
def build_data(tmp_path):
    return make_build_dataset(tmp_path / "data")


def build(data, *extra):
    return main(["build-db", "-c", str(data.config), *extra])


def manifest(data):
    return json.loads((data.outdir / "run_manifest.json").read_text())


def database(data):
    return dict(iter_fasta(data.outdir / manifest(data)["outputs"]["allele_db"]))


@requires_tools
def test_build_db_end_to_end(build_data, capsys):
    genome_folder_before = sorted(path.name for path in build_data.folder.iterdir())
    assert build(build_data) == EXIT_OK
    alleles = database(build_data)

    # SNPs on either strand match the reads, in uppercase even where the genome is soft-masked
    for marker in ("Chr01_000000500", "Chr01_000001000", "Chr02_000000500"):
        ref_read, alt_read = build_data.expected[marker]
        assert alleles[f"{marker}|Ref_0001"] == ref_read
        assert alleles[f"{marker}|Alt_0002"] == alt_read
    assert all(sequence == sequence.upper() for sequence in alleles.values())

    # The genome's base replaces an IUPAC code in the MADC
    assert alleles["Chr01_000001500|Ref_0001"] == build_data.expected["Chr01_000001500"][0]
    assert "Y" not in alleles["Chr01_000001500|Ref_0001"]

    # Indels keep the read's alignment: longer for an insertion, shorter for a deletion
    marker = "Chr01_000002000"
    assert alleles[f"{marker}|Ref_0001"] == build_data.expected[marker][0]
    assert len(alleles[f"{marker}|Alt_0002"]) == AMPLICON + 2
    assert alleles[f"{marker}|Alt_0002"].startswith(build_data.expected[marker][1])
    marker = "Chr01_000002500"
    assert len(alleles[f"{marker}|Alt_0002"]) == AMPLICON - 2
    assert build_data.expected[marker][1].startswith(alleles[f"{marker}|Alt_0002"])

    record = manifest(build_data)
    counts = record["counts"]
    assert record["status"] == "finished"
    assert counts["markers_in_lookup"] == 10
    assert counts["markers_in_database"] == 7
    assert counts["alleles"] == 14
    assert counts["short_snp_alleles"] == []
    assert counts["markers_on_minus_strand"] == 1

    # Inputs are only read, and no index files land next to the genome
    assert sorted(path.name for path in build_data.folder.iterdir() if path.name != "out") == genome_folder_before
    assert "core_db_v001.fa" in capsys.readouterr().out


@requires_tools
def test_process_readme_lists_inputs_and_settings(build_data, capsys):
    assert build(build_data) == EXIT_OK
    assert manifest(build_data)["process_readme"] == "versions/core_db_v001_process.readme"
    readme = (build_data.outdir / "versions" / "core_db_v001_process.readme").read_text()
    header = readme[:readme.index("[1/6]")]
    assert "--- Scripts, files, and parameters used ---" in header
    for name in ("probe.xlsx", "PROJECT_MADC.csv", "genome.fna", "chromosome_names.csv"):
        assert str((build_data.folder / name).resolve()) in header
    assert "# Allele length: 109" in header
    assert "# Flank length: 150" in header
    assert "# Tool version: mmseqs " in header
    assert "Number of markers in probe design file" in readme  # the step output


@requires_tools
def test_problem_markers_are_removed_from_all_outputs(build_data):
    assert build(build_data) == EXIT_OK
    counts = manifest(build_data)["counts"]
    outputs = manifest(build_data)["outputs"]

    assert counts["markers_removed_as_duplicates"] == [build_data.duplicate_pair[1]]
    assert counts["markers_removed_without_alleles"] == list(build_data.placed_badly)
    removed = [build_data.duplicate_pair[1], *build_data.placed_badly]

    alleles = database(build_data)
    assert f"{build_data.duplicate_pair[0]}|Ref_0001" in alleles
    lookup = pd.read_csv(build_data.outdir / outputs["marker_lut"])
    count_table = (build_data.outdir / outputs["match_counts"]).read_text()
    for marker in removed:
        assert f"{marker}|Ref_0001" not in alleles
        assert marker not in set(lookup["Marker_ID"])
        assert marker not in count_table
    assert "Chr02_000003000" not in set(lookup["Marker_ID"])  # in the probe file but not the MADC


@requires_tools
def test_two_runs_are_identical(tmp_path):
    data = make_build_dataset(tmp_path / "data")
    assert build(data) == EXIT_OK
    assert build(data, "-o", str(tmp_path / "again")) == EXIT_OK
    for name in ("core_db_v001.fa", "core_db_v001_matchCnt_lut.txt"):
        assert (data.outdir / "versions" / name).read_bytes() == (tmp_path / "again" / "versions" / name).read_bytes()


@requires_tools
def test_core_database_feeds_assign(build_data, tmp_path, capsys):
    assert build(build_data) == EXIT_OK
    outputs = manifest(build_data)["outputs"]
    settings = {
        "project_id": "PROJECT",
        "inputs": {
            "madc": str(build_data.madc),
            "marker_lut": str(build_data.outdir / outputs["marker_lut"]),
            "allele_db": str(build_data.outdir / outputs["allele_db"]),
            "match_counts": str(build_data.outdir / outputs["match_counts"]),
        },
        "panel": {"design_length": AMPLICON, "sequencing_length": AMPLICON},
        "filter": {"min_samples": 2, "min_reads": 10},
        "blast": {"min_coverage": 80, "min_identity": 85, "threads": 2},
        "output": {"dir": str(tmp_path / "assigned")},
    }
    config = tmp_path / "assign.yaml"
    config.write_text(yaml.safe_dump(settings))
    assert main(["check", "-c", str(config)]) == EXIT_OK
    assert main(["assign", "-c", str(config)]) == EXIT_OK

    record = json.loads((tmp_path / "assigned" / "run_manifest.json").read_text())
    frame = pd.read_csv(tmp_path / "assigned" / record["outputs"]["report"])
    assert len(SAMPLES) == record["counts"]["report_samples"]
    assert set(frame["CloneID"]) == {marker.split("|")[0] for marker in database(build_data)}
    assert record["outputs"]["allele_db"] == "versions/core_db_v002.fa"
    assert frame["AlleleID"].str.contains(r"\|RefMatch_0001$").sum() == 7


@requires_tools
def test_renamed_genome_needs_no_chromosome_map(tmp_path):
    with_map = make_build_dataset(tmp_path / "mapped")
    renamed = make_build_dataset(tmp_path / "renamed", use_genome_names=True)
    assert build(with_map) == EXIT_OK
    assert build(renamed) == EXIT_OK
    assert database(with_map) == database(renamed)


@requires_tools
@pytest.mark.parametrize("probe_format", ["csv", "txt"])
def test_every_probe_format_gives_the_same_database(tmp_path, probe_format):
    from_excel = make_build_dataset(tmp_path / "excel")
    other = make_build_dataset(tmp_path / probe_format, probe_format=probe_format)
    assert build(from_excel) == EXIT_OK
    assert build(other) == EXIT_OK
    assert database(from_excel) == database(other)


@requires_tools
def test_mismatched_chromosome_names_are_caught_before_starting(build_data, capsys):
    settings = build_data.settings()
    del settings["inputs"]["chromosome_map"]
    build_data.write_settings(settings)
    assert main(["check", "-c", str(build_data.config)]) == EXIT_ERROR
    assert "Set inputs.chromosome_map" in capsys.readouterr().out
    assert build(build_data) == EXIT_ERROR
    assert "Chromosome(s) not in the genome" in capsys.readouterr().err
    assert not build_data.outdir.exists()


@requires_tools
def test_a_genome_that_matches_nothing_stops_the_run(build_data, capsys):
    # Right names, wrong sequence
    other = random.Random(5)
    with open(build_data.genome, "w") as out:
        for name in ("NC_1.1", "NC_2.1"):
            out.write(f">{name}\n{''.join(other.choice('ACGT') for _ in range(6000))}\n")
    assert build(build_data) == EXIT_ERROR
    err = capsys.readouterr().err
    assert "the assembly the probes were designed on" in err
    assert manifest(build_data)["status"] == "failed"


@requires_tools
def test_until_stops_early_and_used_folder_is_refused(build_data, capsys):
    assert build(build_data, "--until", "orient") == EXIT_OK
    record = manifest(build_data)
    assert record["status"] == "stopped after orient"
    assert not (build_data.outdir / "versions" / "core_db_v001.fa").exists()
    assert build(build_data) == EXIT_ERROR
    assert "not empty" in capsys.readouterr().err


@requires_tools
def test_short_flanks_warn_about_cut_alleles(build_data, capsys):
    assert build(build_data, "--set", "database.flank_length=60") == EXIT_OK
    assert manifest(build_data)["counts"]["short_snp_alleles"]
    assert "shorter than requested" in capsys.readouterr().out


def test_stage_names_match_the_workflow():
    assert BUILD_DB_STAGES == tuple(STAGES)


def test_build_db_template_is_valid(tmp_path):
    path = tmp_path / "build.yaml"
    assert main(["init", str(path), "--workflow", "build-db"]) == EXIT_OK
    assert path.read_text() == BUILD_DB_TEMPLATE
    cfg = load_config(path, schema=BUILD_DB_SCHEMA)
    assert cfg.database.allele_length == 109
    assert cfg.database.flank_length == 150
    assert cfg.inputs.chromosome_map is None
    assert workflow_of(path) == "build-db"


def test_assign_settings_are_rejected_by_build_db(dataset, capsys):
    assert workflow_of(dataset.config) == "assign"
    assert main(["build-db", "-c", str(dataset.config)]) == EXIT_ERROR
    assert "Unknown setting" in capsys.readouterr().err


def test_duplicate_thresholds_are_validated(build_data):
    with pytest.raises(ConfigError, match="at most 1"):
        load_config(build_data.config, ["database.dup_min_identity=95"], schema=BUILD_DB_SCHEMA)


# --- FASTA reading


@pytest.fixture
def small_fasta(tmp_path):
    path = tmp_path / "t.fa"
    path.write_text(">a first sequence\nACGTACGTAC\nGTACGT\n>b\nTTTTTGGGGG\n")
    return path


def fetch(tmp_path, fasta, key_lines, rename=None):
    keys = tmp_path / "keys.txt"
    keys.write_text(key_lines)
    out = tmp_path / "out.fa"
    fetch_subsequences(fasta, keys, out, rename)
    return list(iter_fasta(out))


def test_read_lengths_and_names(small_fasta):
    assert read_lengths(small_fasta) == {"a": 16, "b": 10}


def test_fetch_keeps_key_order_and_uses_inclusive_coordinates(tmp_path, small_fasta):
    assert fetch(tmp_path, small_fasta, "n1\t6\t10\tb\nn2\t3\t7\ta\n") == [("n1", "GGGGG"), ("n2", "GTACG")]


def test_fetch_can_translate_names(tmp_path, small_fasta):
    assert fetch(tmp_path, small_fasta, "n\t1\t4\tChr1\n", {"Chr1": "a"}) == [("n", "ACGT")]


def test_fetch_cuts_ranges_past_the_end(tmp_path, small_fasta):
    assert fetch(tmp_path, small_fasta, "n\t8\t100\tb\n") == [("n", "GGG")]


def test_fetch_reports_missing_names(tmp_path, small_fasta):
    with pytest.raises(ValueError, match="missing"):
        fetch(tmp_path, small_fasta, "n\t1\t4\tmissing\n")


@requires_esl_sfetch
def test_fetch_matches_esl_sfetch(tmp_path, small_fasta):
    keys = "n1\t3\t7\ta\nn2\t9\t16\ta\nn3\t1\t10\tb\n"
    ours = fetch(tmp_path, small_fasta, keys)
    subprocess.run(["esl-sfetch", "--index", str(small_fasta)], check=True, capture_output=True)
    result = subprocess.run(["esl-sfetch", "-Cf", str(small_fasta), str(tmp_path / "keys.txt")],
                            check=True, capture_output=True, text=True)
    theirs_path = tmp_path / "theirs.fa"
    theirs_path.write_text(result.stdout)
    assert ours == list(iter_fasta(theirs_path))


# --- Probe files


def write_sheet(path, rows):
    book = openpyxl.Workbook()
    for row in rows:
        book.active.append(row)
    book.save(path)


def test_probe_reader_skips_notes_above_the_header(tmp_path):
    path = tmp_path / "probe.xlsx"
    row = ["m1", "SEQ", "g", "Chr01", 120735, "[T/A]", "SNP", 2, "Marker"]
    write_sheet(path, [[None, "notes"], ["Column 1", "Column 2"], PROBE_HEADER, row, [None] * 9])
    header, rows = read_probe_rows(path)
    assert header[0] == "Marker Name"
    assert rows == [["m1", "SEQ", "g", "Chr01", "120735", "[T/A]", "SNP", "2", "Marker"]]


def test_whole_number_cells_are_written_without_a_decimal():
    # Spreadsheet tools can store a position as 120735.0
    assert _cell_text(120735.0) == "120735"
    assert _cell_text(1.5) == "1.5"
    assert _cell_text(None) == ""


def test_probe_reader_needs_a_header(tmp_path):
    path = tmp_path / "probe.csv"
    path.write_text("m1,SEQ,g,Chr01,1,[T/A],SNP\n")
    with pytest.raises(ValueError, match="Marker Name"):
        read_probe_rows(path)


def test_probe_reader_needs_enough_columns(tmp_path):
    path = tmp_path / "probe.csv"
    path.write_text("Marker Name,Seq\nm1,SEQ\n")
    with pytest.raises(ValueError, match="fewer than"):
        read_probe_rows(path)


def test_commas_inside_cells_do_not_shift_columns(tmp_path):
    source = tmp_path / "probe.csv"
    source.write_text(",".join(PROBE_HEADER) + '\nm1,SEQ,g,Chr01,100,[T/A],SNP,2,"Marker, important"\n')
    target = tmp_path / "clean.csv"
    write_probe_csv(*read_probe_rows(source), target)
    prepare_lookup.main([str(target)])
    row = (tmp_path / "clean_snpID_lut.csv").read_text().splitlines()[1].split(",")
    assert row[:7] == ["m1", "Chr01_000000100", "Chr01", "100", "T", "A", "SNP"]
    assert row[-1] == "Marker; important"


# --- Step fixes


def test_db00_writes_the_same_lookup_table_from_a_tab_delimited_probe(tmp_path):
    rows = [PROBE_HEADER, ["m1", "SEQ", "g", "Chr01", "100", "[T/A]", "SNP", "2", "Marker"]]
    (tmp_path / "tabs.txt").write_text("".join("\t".join(row) + "\n" for row in rows))
    (tmp_path / "commas.csv").write_text("".join(",".join(row) + "\n" for row in rows))
    prepare_lookup.main([str(tmp_path / "tabs.txt")])
    prepare_lookup.main([str(tmp_path / "commas.csv")])
    assert (tmp_path / "tabs_snpID_lut.csv").read_text() == (tmp_path / "commas_snpID_lut.csv").read_text()


def test_db02_count_table_has_a_fixed_order(tmp_path):
    markers = [f"Chr{number:02d}_1" for number in range(20, 0, -1)]
    lookup = tmp_path / "lookup.csv"
    lookup.write_text("Panel_markerID,Marker_ID\n" + "".join(f"m{name},{name}\n" for name in markers))
    report = tmp_path / "report.csv"
    report.write_text("AlleleID,CloneID,AlleleSequence,S1\n" + "".join(f"m{name}|Ref,m{name},AAA,1\n" for name in markers))
    update_ids.main([str(lookup), str(report)])
    rows = (tmp_path / "report_snpID_matchCnt_lut.txt").read_text().splitlines()
    written = [row.split("|")[0] for row in rows]
    assert written == sorted(written)
    assert len(set(written)) == len(markers)
