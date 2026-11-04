"""Checks made before a run starts and after the allele database is updated."""

from __future__ import annotations

import os
import re
from collections import Counter
from dataclasses import dataclass

from .fasta import read_lengths
from .madc import longest_ref_read, read_layout, read_marker_ids
from .probe_design import probe_row_problem, read_probe_rows
from .steps.assign import step01_filter_missing_AND_ext_matchAlleles_from_madc_v1 as filter_alleles
from .tools import BUILD_DB_TOOLS, required_tools, tool_version

OK, WARNING, PROBLEM = "ok", "warning", "problem"
NOT_TEXT_ADVICE = "Save it as plain text or CSV, and decompress it if it is gzipped."


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


def _text_lines(path, setting, findings):
    """A text input's lines, or None with a problem recorded if it is not text."""
    try:
        with open(path, encoding="utf-8-sig") as handle:
            return handle.read().splitlines()
    except UnicodeDecodeError:
        findings.append(Finding(PROBLEM, f"{setting} is not a text file: {path}. {NOT_TEXT_ADVICE}"))
        return None


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
            reason = "is a folder, not a file" if path.exists() else "not found"
            findings.append(Finding(PROBLEM, f"inputs.{name} {reason}: {path}"))
            all_found = False
    return all_found


def _check_madc(cfg, findings):
    try:
        layout = read_layout(cfg.inputs.madc)
    except ValueError as err:
        findings.append(Finding(PROBLEM, str(err)))
        return None

    detected = layout.first_sample_column
    chosen = cfg.panel.first_sample_column
    first_column = chosen or detected
    n_samples = len(layout.samples) + detected - first_column
    findings.append(Finding(OK, f"MADC: {n_samples} samples, first sample in column {first_column}"))
    if n_samples <= 0:
        findings.append(Finding(PROBLEM, "MADC header has no sample columns"))
    if chosen and chosen < detected:
        findings.append(Finding(PROBLEM, f"panel.first_sample_column is {chosen}, but columns up to {detected - 1} of "
                                         "the MADC header are summary columns, not samples. Remove this setting to "
                                         "detect the column automatically."))
    elif chosen and chosen > detected:
        findings.append(Finding(WARNING, f"panel.first_sample_column is {chosen}, but the header suggests {detected}, "
                                         f"so the first {chosen - detected} sample column(s) will be left out"))
    read_length = longest_ref_read(cfg.inputs.madc)
    if read_length and read_length != cfg.panel.sequencing_length:
        findings.append(Finding(WARNING, f"panel.sequencing_length is {cfg.panel.sequencing_length}, "
                                         f"but the MADC's Ref alleles are {read_length} bp"))
    repeated = sorted(name for name, count in Counter(layout.samples).items() if count > 1)
    if repeated:
        findings.append(Finding(WARNING, f"Repeated sample names will get a numeric suffix: {repeated[:5]}"))
    if cfg.filter.min_samples is None:
        required = filter_alleles.required_samples(None, n_samples)
        findings.append(Finding(OK, f"filter.min_samples is not set, so the filter uses {required} samples: the smaller "
                                    f"of {filter_alleles.DEFAULT_MIN_SAMPLES} and 5% of {n_samples} samples"))
    elif cfg.filter.min_samples > n_samples:
        findings.append(Finding(PROBLEM, f"filter.min_samples is {cfg.filter.min_samples}, but the MADC has "
                                         f"{n_samples} samples, so no allele could pass"))
    return layout


def _check_lookup_table(cfg, findings):
    """Check the lookup table. Returns the new IDs of the markers this run will use, or None."""
    lines = _text_lines(cfg.inputs.marker_lut, "inputs.marker_lut", findings)
    if lines is None:
        return None
    too_short = [number for number, line in enumerate(lines, 1) if len(line.split(",")) < 2]
    if too_short:
        findings.append(Finding(PROBLEM, f"inputs.marker_lut needs at least two comma-separated columns, but line "
                                         f"{too_short[0]} does not ({len(too_short)} such line(s), blank lines count)"))
        return None
    rows = [line.split(",") for line in lines]
    new_id_of = {row[0]: row[1] for row in rows}

    madc_markers = read_marker_ids(cfg.inputs.madc)
    in_lookup = madc_markers & set(new_id_of)
    if not in_lookup:
        findings.append(Finding(PROBLEM, f"None of the MADC's markers are in the lookup table. Its first column should "
                                         f"hold the MADC's marker names, such as {sorted(madc_markers)[:1]}."))
        return None
    not_in_lookup = sorted(madc_markers - in_lookup)
    if not_in_lookup:
        findings.append(Finding(WARNING, f"{len(not_in_lookup)} MADC marker(s) are not in the lookup table "
                                         f"and will be dropped: {not_in_lookup[:5]}"))
    column = cfg.panel.insert_size_column
    if column and (not rows or column not in rows[0]):
        findings.append(Finding(PROBLEM, f"panel.insert_size_column {column!r} is not in the lookup table header"))

    tags = set()
    if cfg.inputs.dup_tags:
        tag_lines = _text_lines(cfg.inputs.dup_tags, "inputs.dup_tags", findings)
        tags = {line.strip() for line in tag_lines or [] if line.strip()}
        if tag_lines is not None and not tags:
            findings.append(Finding(WARNING, "inputs.dup_tags is empty, so no markers will be removed"))
        unknown = sorted(tags - set(new_id_of.values()))
        if unknown:
            findings.append(Finding(WARNING, f"{len(unknown)} marker(s) in inputs.dup_tags are not in the "
                                             f"lookup table's new IDs, so nothing removes them: {unknown[:5]}"))
    return {new_id_of[marker] for marker in in_lookup} - tags


def _check_database(cfg, findings, marker_ids):
    if not re.search(r"_v\d{1,3}(?=[_.]|$)", cfg.inputs.allele_db.name, flags=re.IGNORECASE):
        findings.append(Finding(WARNING, f"{cfg.inputs.allele_db.name} has no version tag such as _v003, "
                                         "so a new database would be named _v001"))
    try:
        fasta_ids = read_fasta_ids(cfg.inputs.allele_db)
    except UnicodeDecodeError:
        findings.append(Finding(PROBLEM, f"inputs.allele_db is not a text file: {cfg.inputs.allele_db}. {NOT_TEXT_ADVICE}"))
        return
    if not fasta_ids:
        findings.append(Finding(PROBLEM, f"inputs.allele_db has no FASTA records: {cfg.inputs.allele_db}"))
        return
    if marker_ids is not None:
        missing = sorted(marker for marker in marker_ids
                         if f"{marker}|Ref_0001" not in fasta_ids or f"{marker}|Alt_0002" not in fasta_ids)
        if missing:
            findings.append(Finding(PROBLEM, f"{len(missing)} of {len(marker_ids)} markers have no Ref_0001 and "
                                             f"Alt_0002 allele in inputs.allele_db, such as {missing[:3]}. Use the "
                                             "database built for this panel, or list these markers in inputs.dup_tags."))
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
    if outdir.exists():
        if not outdir.is_dir():
            findings.append(Finding(PROBLEM, f"Output folder is a file: {outdir}"))
            return
        if any(outdir.iterdir()):
            findings.append(Finding(PROBLEM, f"Output folder is not empty: {outdir}"))
            return
        nearest = outdir
    else:
        nearest = next(parent for parent in outdir.parents if parent.exists())
    if not os.access(nearest, os.W_OK):
        findings.append(Finding(PROBLEM, f"Cannot write to the output folder {outdir}: no permission in {nearest}"))
    else:
        findings.append(Finding(OK, f"Output folder: {outdir}"))


def preflight(cfg):
    """Everything that can be checked before starting. Returns (findings, MADC layout or None)."""
    findings = []
    _check_tools(required_tools(cfg), findings)
    layout = None
    if _check_inputs_exist(cfg, findings):
        layout = _check_madc(cfg, findings)
        if layout is not None:
            marker_ids = _check_lookup_table(cfg, findings)
            _check_database(cfg, findings, marker_ids)
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
    row_problems = [problem for problem in map(probe_row_problem, used) if problem]
    if row_problems:
        findings.append(Finding(PROBLEM, f"{len(row_problems)} marker(s) in the probe design file cannot be used: "
                                         f"{row_problems[0]}"))
    not_in_probe = sorted(madc_markers - {row[0] for row in rows})
    if not_in_probe:
        findings.append(Finding(WARNING, f"{len(not_in_probe)} MADC marker(s) are not in the probe design "
                                         f"file and will be left out: {not_in_probe[:5]}"))

    try:
        names = read_chromosome_map(cfg.inputs.chromosome_map)
    except UnicodeDecodeError:
        findings.append(Finding(PROBLEM, f"inputs.chromosome_map is not a text file: {cfg.inputs.chromosome_map}. "
                                         f"{NOT_TEXT_ADVICE}"))
        return
    try:
        genome = read_lengths(cfg.inputs.genome)
    except UnicodeDecodeError:
        findings.append(Finding(PROBLEM, f"inputs.genome is not a text file: {cfg.inputs.genome}. {NOT_TEXT_ADVICE}"))
        return
    if not genome:
        findings.append(Finding(PROBLEM, f"inputs.genome has no FASTA records: {cfg.inputs.genome}"))
        return
    chromosomes = sorted({row[3] for row in used})
    unmatched = [name for name in chromosomes if names.get(name, name) not in genome]
    if unmatched:
        hint = "" if cfg.inputs.chromosome_map else " Set inputs.chromosome_map if the genome names them differently."
        findings.append(Finding(PROBLEM, f"Chromosome(s) not in the genome: {unmatched[:5]}, while the genome has "
                                         f"names like {list(genome)[:3]}.{hint}"))
    else:
        findings.append(Finding(OK, f"Genome: {len(genome)} sequences, including all {len(chromosomes)} "
                                    "chromosomes these markers use"))

    length = cfg.database.allele_length
    if cfg.database.flank_length < length:
        findings.append(Finding(WARNING, f"database.flank_length {cfg.database.flank_length} is shorter than "
                                         f"database.allele_length {length}, so some alleles may be cut short"))
    read_length = longest_ref_read(cfg.inputs.madc)
    if read_length and read_length < length:
        findings.append(Finding(OK, f"Reads of {read_length} bp will be extended from the genome to {length} bp alleles"))
    elif read_length > length:
        findings.append(Finding(WARNING, f"database.allele_length {length} is shorter than the MADC's "
                                         f"{read_length} bp Ref alleles, so alleles will be shorter than the reads"))


def preflight_build_db(cfg):
    """Everything that can be checked before starting the core database workflow."""
    findings = []
    _check_tools(BUILD_DB_TOOLS, findings)
    if _check_inputs_exist(cfg, findings):
        _check_probe_and_genome(cfg, findings)
    _check_output_folder(cfg, findings)
    return findings
