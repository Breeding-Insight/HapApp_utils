"""A small synthetic MADC, allele database and lookup table, built fresh for each test."""

import random
import shutil
from dataclasses import dataclass
from pathlib import Path

import pytest
import yaml

SAMPLES = [f"S{number:02d}" for number in range(1, 41)]
SUMMARY_COLUMNS = ["ClusterConsensusSequence", "CallRate", "OneRatioRef", "OneRatioSnp", "FreqHomRef", "FreqHomSnp",
                   "FreqHets", "PICRef", "PICSnp", "AvgPIC", "AvgCountRef", "AvgCountSnp", "RatioAvgCountRefAvgCountSnp"]
MARKERS = ["Chr01_000000100", "Chr01_000000200", "Chr02_000000300"]
ADAPTER = "ACCGATCTG"

requires_blast = pytest.mark.skipif(shutil.which("blastn") is None, reason="BLAST is not installed")
requires_cutadapt = pytest.mark.skipif(shutil.which("cutadapt") is None, reason="cutadapt is not installed")


def mutate(sequence, *positions):
    bases = list(sequence)
    for position in positions:
        bases[position] = {"A": "C", "C": "G", "G": "T", "T": "A"}[bases[position]]
    return "".join(bases)


def reads_in(n_samples, reads):
    return [reads] * n_samples + [0] * (len(SAMPLES) - n_samples)


@dataclass
class Dataset:
    folder: Path
    config: Path
    madc: Path
    allele_db: Path
    match_counts: Path
    lookup: Path
    outdir: Path
    noise: dict
    new_alleles: dict

    def settings(self):
        return yaml.safe_load(self.config.read_text())

    def write_settings(self, settings):
        self.config.write_text(yaml.safe_dump(settings))


def make_dataset(folder, *, length=109, trim_design=None, duplicate_after_trim=False,
                 raw_format=True, duplicate_allele=False, insert_sizes=None):
    """Build the input files.

    Each marker has Ref and Alt, a RefMatch already in the database, one new allele of each class,
    and a noise allele with reads in only one sample. trim_design makes reads carry an adapter
    after that many bases. insert_sizes gives each marker its own design length instead.
    """
    rng = random.Random(7)
    folder.mkdir(parents=True, exist_ok=True)
    rows = []
    database = []
    counts = []
    noise = {}
    new_alleles = {}

    def to_read(genomic, tail_seed, read_length):
        """A sequence as reported: genomic bases, then adapter and tail if the read runs past the insert."""
        if len(genomic) >= read_length:
            return genomic
        tail_rng = random.Random(tail_seed)
        tail = "".join(tail_rng.choice("ACGT") for _ in range(read_length))
        return (genomic + ADAPTER + tail)[:read_length]

    for marker in MARKERS:
        design = (insert_sizes or {}).get(marker, trim_design or length)
        ref = "".join(rng.choice("ACGT") for _ in range(design))
        alt = mutate(ref, 54)
        known = mutate(ref, 20)
        new_ref = mutate(ref, 30)
        new_alt = mutate(alt, 70)
        new_other = mutate(ref, 10, 75)
        noisy = mutate(ref, 40)
        noise[marker] = noisy
        new_alleles[marker] = {"RefMatch": new_ref, "AltMatch": new_alt, "Other": new_other}

        def row(allele_id, genomic, tail_seed, counts_by_sample):
            rows.append((f"{marker}|{allele_id}", marker, to_read(genomic, tail_seed, length), counts_by_sample))

        row("Ref", ref, 1, [50] * len(SAMPLES))
        row("Alt", alt, 2, [30] * len(SAMPLES))
        row("RefMatch", known, 3, reads_in(10, 20))
        row("RefMatch", new_ref, 4, reads_in(5, 15))
        if duplicate_after_trim:
            row("RefMatch", new_ref, 44, reads_in(3, 11))
        row("AltMatch", new_alt, 5, reads_in(3, 12))
        row("Other", new_other, 6, reads_in(4, 14))
        row("RefMatch", noisy, 7, reads_in(1, 40))
        database += [(f"{marker}|Ref_0001", ref), (f"{marker}|Alt_0002", alt), (f"{marker}|RefMatch_0001", known)]
        counts += [f"{marker}|RefMatch\t1\t0\n", f"{marker}|AltMatch\t0\t0\n"]

    if duplicate_allele:
        original = rows[2]
        rows.insert(3, (f"{MARKERS[0]}|Other", MARKERS[0], original[2], original[3]))

    madc = folder / "PROJECT_MADC.csv"
    lines = []
    if raw_format:
        filler = ["*"] * (3 + len(SUMMARY_COLUMNS))
        lines += [",".join(filler + ["x"] * len(SAMPLES)) + "\n" for _ in range(7)]
        lines.append(",".join(["AlleleID", "CloneID", "AlleleSequence"] + SUMMARY_COLUMNS + SAMPLES) + "\n")
    else:
        lines.append(",".join(["AlleleID", "CloneID", "AlleleSequence"] + SAMPLES) + "\n")
    for allele_id, marker, sequence, counts_by_sample in rows:
        cells = [allele_id, marker, sequence]
        if raw_format:
            # AvgCountRef and AvgCountSnp pass any read threshold, which is what the filter must ignore
            cells += [sequence, "0.5", "0.5", "0.5", "0", "0", "0", "0", "0", "0", "25.0", "25.0", "1.0"]
        lines.append(",".join(cells + [str(value) for value in counts_by_sample]) + "\n")
    madc.write_text("".join(lines))

    lookup = folder / "marker_lookup.csv"
    if insert_sizes:
        body = "".join(f"{marker},{marker},{insert_sizes[marker]}\n" for marker in MARKERS)
        lookup.write_text("Panel_markerID,markerID,insert_size\n" + body)
    else:
        lookup.write_text("".join(f"{marker},{marker},chr,1,A,T,SNP\n" for marker in MARKERS))
    allele_db = folder / "test_allele_db_v001.fa"
    allele_db.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in database))
    match_counts = folder / "test_allele_db_v001_matchCnt_lut.txt"
    match_counts.write_text("".join(counts))

    config = folder / "run.yaml"
    outdir = folder / "out"
    settings = {
        "project_id": "PROJECT",
        "inputs": {"madc": madc.name, "marker_lut": lookup.name, "allele_db": allele_db.name,
                   "match_counts": match_counts.name},
        "panel": {"design_length": trim_design or length, "sequencing_length": length},
        "filter": {"min_samples": 2, "min_reads": 10},
        "blast": {"min_coverage": 80, "min_identity": 85, "threads": 2},
        "output": {"dir": "out"},
    }
    config.write_text(yaml.safe_dump(settings))
    return Dataset(folder, config, madc, allele_db, match_counts, lookup, outdir, noise, new_alleles)


@pytest.fixture
def dataset(tmp_path):
    return make_dataset(tmp_path / "data")
