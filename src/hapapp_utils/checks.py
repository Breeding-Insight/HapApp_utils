"""Checks made before a run starts and after the allele database is updated."""

from __future__ import annotations

import math
import re
from collections import Counter
from dataclasses import dataclass

from .fasta import read_lengths
from .madc import read_layout, read_marker_ids
from .probe_design import read_probe_rows
from .tools import BUILD_DB_TOOLS, required_tools, tool_version

OK, WARNING, PROBLEM = "ok", "warning", "problem"


@dataclass
class Finding:
    level: str
    message: str


def read_fasta_ids(path):
    with open(path) as handle:
        return {line[1:].split()[0] for line in handle if line.startswith(">")}


def read_counts(path):
    """A match-count table: allele class at a marker -> number of IDs issued so far."""
    counts = {}
    with open(path) as handle:
        for line in handle:
            cells = line.strip().split("\t")
            if len(cells) >= 2:
                counts[cells[0]] = int(cells[1])
    return counts


def minted_ids_missing_from_fasta(old_counts_path, new_counts_path, fasta_path):
    """IDs issued by this run with no sequence in the new database.

    Only this run's IDs are checked, so a run on an older database with gaps can still finish
    and be judged on its own work. Older gaps are reported by counted_ids_missing_from_fasta.
    """
    old_counts = read_counts(old_counts_path)
    fasta_ids = read_fasta_ids(fasta_path)
    missing = []
    for key, count in read_counts(new_counts_path).items():
        for serial in range(old_counts.get(key, 0) + 1, count + 1):
            allele_id = f"{key}_{serial:04d}"
            if allele_id not in fasta_ids:
                missing.append(allele_id)
    return missing


def counted_ids_missing_from_fasta(counts_path, fasta_path):
    """All IDs a count table says were issued that are not in the database."""
    fasta_ids = read_fasta_ids(fasta_path)
    return [f"{key}_{serial:04d}"
            for key, count in read_counts(counts_path).items()
            for serial in range(1, count + 1)
            if f"{key}_{serial:04d}" not in fasta_ids]


def _check_tools(tools, findings):
    for tool in tools:
        version = tool_version(tool)
        if version is None:
            findings.append(Finding(PROBLEM, f"{tool} is not installed or not on PATH"))
        else:
            findings.append(Finding(OK, f"{tool}: {version}"))


def _check_inputs_exist(cfg, findings):
    all_found = True
    for name, path in vars(cfg.inputs).items():
        if path is None:
            continue
        if path.is_file():
            findings.append(Finding(OK, f"inputs.{name}: {path}"))
        else:
            findings.append(Finding(PROBLEM, f"inputs.{name} not found: {path}"))
            all_found = False
    return all_found


def _check_madc(cfg, findings):
    try:
        layout = read_layout(cfg.inputs.madc)
    except ValueError as err:
        findings.append(Finding(PROBLEM, str(err)))
        return None

    n_samples = len(layout.samples)
    first_column = cfg.panel.first_sample_column or layout.first_sample_column
    findings.append(Finding(OK, f"MADC: {n_samples} samples, first sample in column {first_column}"))
    if n_samples == 0:
        findings.append(Finding(PROBLEM, "MADC header has no sample columns"))
    if cfg.panel.first_sample_column and cfg.panel.first_sample_column != layout.first_sample_column:
        findings.append(Finding(WARNING, f"panel.first_sample_column is {cfg.panel.first_sample_column}, "
                                         f"but the header suggests {layout.first_sample_column}"))
    repeated = sorted(name for name, count in Counter(layout.samples).items() if count > 1)
    if repeated:
        findings.append(Finding(WARNING, f"Repeated sample names will get a numeric suffix: {repeated[:5]}"))
    effective = min(cfg.filter.min_samples, math.ceil(n_samples * 0.05))
    if effective < cfg.filter.min_samples:
        findings.append(Finding(OK, f"filter.min_samples is {cfg.filter.min_samples}, and the 5% rule "
                                    f"lowers it to {effective} for {n_samples} samples"))
    return layout


def _check_lookup_table(cfg, findings):
    with open(cfg.inputs.marker_lut, encoding="utf-8-sig") as handle:
        rows = [line.rstrip("\r\n").split(",") for line in handle]
    rows = [row for row in rows if len(row) >= 2]
    old_ids = {row[0] for row in rows}
    new_ids = {row[1] for row in rows}

    not_in_lookup = sorted(read_marker_ids(cfg.inputs.madc) - old_ids)
    if not_in_lookup:
        findings.append(Finding(WARNING, f"{len(not_in_lookup)} MADC marker(s) are not in the lookup table "
                                         f"and will be dropped: {not_in_lookup[:5]}"))
    column = cfg.panel.insert_size_column
    if column and (not rows or column not in rows[0]):
        findings.append(Finding(PROBLEM, f"panel.insert_size_column {column!r} is not in the lookup table header"))
    if cfg.inputs.dup_tags:
        with open(cfg.inputs.dup_tags, encoding="utf-8-sig") as handle:
            tags = {line.strip() for line in handle if line.strip()}
        if not tags:
            findings.append(Finding(WARNING, "inputs.dup_tags is empty, so no markers will be removed"))
        unknown = sorted(tags - new_ids)
        if unknown:
            findings.append(Finding(WARNING, f"{len(unknown)} marker(s) in inputs.dup_tags are not in the "
                                             f"lookup table's new IDs, so nothing removes them: {unknown[:5]}"))


def _check_database(cfg, findings):
    if not re.search(r"_v\d{1,3}(?=[_.]|$)", cfg.inputs.allele_db.name, flags=re.IGNORECASE):
        findings.append(Finding(WARNING, f"{cfg.inputs.allele_db.name} has no version tag such as _v003, "
                                         "so a new database would be named _v001"))
    try:
        gaps = counted_ids_missing_from_fasta(cfg.inputs.match_counts, cfg.inputs.allele_db)
    except ValueError as err:
        findings.append(Finding(PROBLEM, f"Could not read inputs.match_counts as a count table: {err}"))
        return
    if gaps:
        findings.append(Finding(WARNING, f"{len(gaps)} allele ID(s) in the count table have no sequence in the "
                                         f"database, from an earlier run: {gaps[:3]}"))
    else:
        findings.append(Finding(OK, "Count table and allele database agree"))


def _check_output_folder(cfg, findings):
    outdir = cfg.output.dir
    if outdir.exists() and any(outdir.iterdir()):
        findings.append(Finding(PROBLEM, f"Output folder is not empty: {outdir}"))
    else:
        findings.append(Finding(OK, f"Output folder: {outdir}"))


def preflight(cfg):
    """Everything that can be checked before starting. Returns (findings, MADC layout or None)."""
    findings = []
    _check_tools(required_tools(cfg), findings)
    layout = None
    if _check_inputs_exist(cfg, findings):
        layout = _check_madc(cfg, findings)
        _check_lookup_table(cfg, findings)
        _check_database(cfg, findings)
    _check_output_folder(cfg, findings)
    return findings, layout


def read_chromosome_map(path):
    """A two-column CSV of probe chromosome name, genome chromosome name."""
    if path is None:
        return {}
    names = {}
    with open(path, encoding="utf-8-sig") as handle:
        for line in handle:
            cells = [cell.strip() for cell in line.split(",")]
            if len(cells) >= 2 and cells[0] and cells[1]:
                names[cells[0]] = cells[1]
    return names


def _check_probe_and_genome(cfg, findings):
    try:
        _, rows = read_probe_rows(cfg.inputs.probe_design)
        layout = read_layout(cfg.inputs.madc)
    except ValueError as err:
        findings.append(Finding(PROBLEM, str(err)))
        return
    madc_markers = read_marker_ids(cfg.inputs.madc)
    used = [row for row in rows if row[0] in madc_markers]
    findings.append(Finding(OK, f"Probe design: {len(rows)} markers, {len(used)} of them in the MADC "
                                f"({len(madc_markers)} markers, {len(layout.samples)} samples)"))
    if not used:
        findings.append(Finding(PROBLEM, "No marker is in both the probe design file and the MADC"))
        return
    not_in_probe = sorted(madc_markers - {row[0] for row in rows})
    if not_in_probe:
        findings.append(Finding(WARNING, f"{len(not_in_probe)} MADC marker(s) are not in the probe design "
                                         f"file and will be left out: {not_in_probe[:5]}"))

    names = read_chromosome_map(cfg.inputs.chromosome_map)
    genome = read_lengths(cfg.inputs.genome)
    chromosomes = sorted({row[3] for row in used})
    unmatched = [name for name in chromosomes if names.get(name, name) not in genome]
    if unmatched:
        hint = "" if cfg.inputs.chromosome_map else " Set inputs.chromosome_map if the genome names them differently."
        findings.append(Finding(PROBLEM, f"Chromosome(s) not in the genome: {unmatched[:5]}, while the genome has "
                                         f"names like {list(genome)[:3]}.{hint}"))
    else:
        findings.append(Finding(OK, f"Genome: {len(genome)} sequences, including all {len(chromosomes)} "
                                    "chromosomes these markers use"))

    if cfg.database.flank_length < cfg.database.allele_length:
        findings.append(Finding(WARNING, f"database.flank_length {cfg.database.flank_length} is shorter than "
                                         f"database.allele_length {cfg.database.allele_length}, "
                                         "so some alleles may be cut short"))


def preflight_build_db(cfg):
    """Everything that can be checked before starting the core database workflow."""
    findings = []
    _check_tools(BUILD_DB_TOOLS, findings)
    if _check_inputs_exist(cfg, findings):
        _check_probe_and_genome(cfg, findings)
    _check_output_folder(cfg, findings)
    return findings
