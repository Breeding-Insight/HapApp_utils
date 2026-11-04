"""The MADC workflow: from a report to fixed allele IDs and a new database version.

This calls the step modules in order, as workflow/20261020-build02_madc_matchhaps_hapsearch.bash
does. Each run works on its own copies of the inputs inside a new folder, and the steps run from
inside that folder so the paths they rewrite never include the user's own directory names.
"""

from __future__ import annotations

import contextlib
import glob
import os
import re
import shutil
from pathlib import Path

from ..checks import PROBLEM, WARNING, minted_ids_missing_from_fasta, preflight, read_fasta_ids
from ..config import SCHEMA
from ..fasta import count_records, iter_fasta
from ..madc import count_match_rows, read_marker_ids
from ..steps.assign import (
    step00_check_allele_uniqueness as check_unique,
    step00_madc_update_snpID_v1_1 as update_ids,
    step00_rm_dupTag_closeTag_rawMADC_v1 as remove_tags,
    step01_filter_missing_AND_ext_matchAlleles_from_madc_v1 as filter_alleles,
    step03_check_cutadapt_allele_uniqueness_AND_update_tmp_rename_report_v1_1 as update_trimmed,
    step05_parse_madc_matchAlleles_blastn as assign_ids,
    step06_check_db_allele_uniqueness_v1 as check_db,
    step06_update_MADC_with_allele_uniqueness_v1 as update_madc_ids,
)
from ..steps.hap_pa import (
    util_hapDB_merge_presence_absence_PIhaps as merge_pa,
    util_haplotype_presence_absence_matrix as hap_pa,
)
from ..tools import required_tools, run_tool
from .common import BLAST_COLUMNS, CannotStart, PipelineError, call_step, run_stages

STAGES = {
    "update-ids": "Update marker IDs and check for duplicate alleles",
    "filter": "Filter alleles on read counts",
    "trim": "Remove adapters",
    "blast": "BLAST alleles against the database",
    "assign-ids": "Assign fixed allele IDs",
    "db-check": "Check the new database",
    "hap-pa": "Build presence/absence tables",
}


class DuplicateAllelesFound(PipelineError):
    """The MADC has the same sequence under different allele IDs."""


class _AssignRun:
    """One run. Paths are relative to the run folder, as the steps expect to rewrite them."""

    def __init__(self, cfg, first_sample_col, n_samples):
        self.cfg = cfg
        self.first_sample_col = str(first_sample_col)
        self.n_samples = n_samples
        self.report = f"report/{cfg.inputs.madc.stem}.csv"
        self.lut = f"lut/{cfg.inputs.marker_lut.name}"
        self.db = f"versions/{cfg.inputs.allele_db.stem}.fa"
        # step06 expects the count table to be named after the database
        self.counts = f"{self.db[:-3]}_matchCnt_lut.txt"
        self.readme = f"{assign_ids.bump_version_tag(self.db)[:-3]}_process.readme"
        self.no_new_allele_readme = f"{self.db[:-3]}_noNewAllele.readme"
        self.snpid = None
        self.base = None
        self.ids_report = None
        self.blast_query = None
        self.blast_out = None
        self.final_report = None
        self.madc_cleaned = None
        self.final_db = self.db
        self.final_counts = self.counts
        self.hap_pa_matrix = None
        self.hap_status = None
        self.warnings = []
        self.summary = {}

    def readme_header(self):
        cfg = self.cfg
        return [
            f"# Scripts directory: {Path(check_unique.__file__).parent}",
            f"# No-new-allele README: {self.no_new_allele_readme}",
            f"# Report: {cfg.inputs.madc}",
            f"# SNP ID LUT: {cfg.inputs.marker_lut}",
            f"# Allele DB input: {cfg.inputs.allele_db}",
            f"# Match-count LUT input: {cfg.inputs.match_counts}",
            f"# Duplicate tags file: {cfg.inputs.dup_tags or 'none'}",
            f"# Earlier presence/absence matrix: {cfg.inputs.previous_hap_pa or 'none'}",
            f"# First sample column: {self.first_sample_col}",
            f"# Panel design length: {cfg.panel.design_length}",
            f"# Sequencing length: {cfg.panel.sequencing_length}",
            f"# Insert size column: {cfg.panel.insert_size_column or 'none'}",
            f"# Adapter: {cfg.panel.adapter}",
            f"# Read filter: {self.filter_rule()}"
            + ("" if cfg.filter.min_samples is not None else " (default: 10 samples, or 5% of samples if smaller)"),
            f"# BLAST filtering thresholds: {cfg.blast.min_coverage}% coverage, {cfg.blast.min_identity}% identity",
            f"# Allele type ground truth: {cfg.blast.allele_type_truth}",
            f"# Reads for presence: {cfg.hap_pa.min_reads}",
            f"# Project descriptor: {cfg.hap_pa.descriptor or cfg.project_id}",
            f"# Code version: {cfg.output.code_version}",
        ]

    def copy_inputs(self):
        for folder in ("report", "versions", "lut"):
            Path(folder).mkdir(exist_ok=True)
        shutil.copyfile(self.cfg.inputs.madc, self.report)
        shutil.copyfile(self.cfg.inputs.marker_lut, self.lut)
        shutil.copyfile(self.cfg.inputs.allele_db, self.db)
        shutil.copyfile(self.cfg.inputs.match_counts, self.counts)

    def stage_update_ids(self):
        cfg = self.cfg
        call_step(update_ids, [self.lut, self.report])
        stem = self.report[:-4]
        self.snpid = f"{stem}_snpID.csv"
        if cfg.inputs.dup_tags and cfg.inputs.dup_tags.stat().st_size > 0:
            tags = f"lut/{cfg.inputs.dup_tags.name}"
            shutil.copyfile(cfg.inputs.dup_tags, tags)
            without_tags = f"{stem}_snpID_rmDup.csv"
            call_step(remove_tags, [tags, self.snpid, without_tags])
            self.snpid = without_tags
        self.base = self.snpid[:-4]
        call_step(check_unique, [self.snpid, self.first_sample_col])
        duplicates = f"{self.base}_dup.csv"
        if os.path.exists(duplicates):
            raise DuplicateAllelesFound(
                f"The MADC has the same sequence under different allele IDs. Review {duplicates} in the run folder.")
        before, after = len(read_marker_ids(self.report)), len(read_marker_ids(self.snpid))
        self.summary["markers"] = after
        kept = f"{after} of {before} markers kept" if after != before else f"{after} markers"
        return f"{kept}, no duplicate alleles"

    def filter_rule(self):
        """The rule the read filter applies in this run, in words."""
        samples = filter_alleles.required_samples(self.cfg.filter.min_samples, self.n_samples)
        reads = self.cfg.filter.min_reads
        return (f"at least {samples} sample{'' if samples == 1 else 's'} "
                f"with at least {reads} read{'' if reads == 1 else 's'}")

    def stage_filter(self):
        cfg = self.cfg
        # Without a chosen min_samples, the step applies its default rule
        chosen_samples = [] if cfg.filter.min_samples is None else ["--min_samples", str(cfg.filter.min_samples)]
        call_step(filter_alleles, [
            self.snpid, "--lut", self.lut,
            *chosen_samples,
            "--min_reads", str(cfg.filter.min_reads),
            "--first_sample_col", self.first_sample_col,
            "--split_by_colID", cfg.panel.insert_size_column or "",
        ])
        self.ids_report = f"{self.base}_tmp_rename.csv"
        self.blast_query = f"{self.base}_match.fa"
        self.final_report = f"{self.base}_rename.csv"
        pieces = glob.glob(f"{self.base}_match_*.fa") if cfg.panel.insert_size_column else [self.blast_query]
        kept = sum(count_records(piece) for piece in pieces)
        tested = count_match_rows(self.snpid)
        rule = self.filter_rule()
        if kept == 0:
            raise PipelineError(f"No allele passed the read filter of {rule}. "
                                "Lower filter.min_samples or filter.min_reads.")
        self.summary["alleles_tested"], self.summary["alleles_kept_by_filter"] = tested, kept
        self.summary["filter_samples_required"] = filter_alleles.required_samples(cfg.filter.min_samples, self.n_samples)
        return f"{kept} of {tested} alleles kept ({rule})"

    def _cutadapt(self, in_fasta, out_fasta):
        panel = self.cfg.panel
        run_tool(["cutadapt", "-a", panel.adapter, "-e", "0.2", "-n", "1", "--overlap", "5",
                  "-m", str(panel.design_length), "-o", out_fasta, in_fasta],
                 stdout_path=f"{out_fasta[:-3]}.log")

    def _trim_mixed_inserts(self, out_fasta):
        """Panels whose markers have different insert sizes: trim only those shorter than the read."""
        pieces = sorted(glob.glob(f"{self.base}_match_*.fa"))
        if not pieces:
            raise PipelineError("No FASTA files were written for panel.insert_size_column. "
                                "Check that column in the lookup table.")
        keep = []
        for piece in pieces:
            label = piece[len(f"{self.base}_match_"):-3]
            size = re.search(r"\d+", label)
            if size is None:
                raise PipelineError(f"Could not read an insert size from {label!r} in the lookup table")
            if int(size.group()) >= self.cfg.panel.sequencing_length:
                keep.append(piece)
            else:
                trimmed = f"{piece[:-3]}_cutadapt.fa"
                self._cutadapt(piece, trimmed)
                keep.append(trimmed)
        with open(out_fasta, "w") as combined:
            for piece in keep:
                with open(piece) as handle:
                    shutil.copyfileobj(handle, combined)

    def stage_trim(self):
        panel = self.cfg.panel
        if panel.design_length >= panel.sequencing_length:
            print("Design length is at least the sequencing length, so no adapters to remove.")
            return "skipped: reads are not longer than the design"
        trimmed = f"{self.base}_match_cutadapt.fa"
        if panel.insert_size_column:
            self._trim_mixed_inserts(trimmed)
        else:
            self._cutadapt(f"{self.base}_match.fa", trimmed)

        before = self.summary["alleles_kept_by_filter"]
        lengths = [len(sequence) for _, sequence in iter_fasta(trimmed)]
        after = len(lengths)
        if after == 0:
            raise PipelineError(f"Trimming left no alleles of at least {panel.design_length} bp. "
                                "Check panel.design_length and panel.adapter.")
        dropped = before - after
        if dropped:
            self.warnings.append(f"{dropped} of {before} alleles were dropped for being shorter than "
                                 f"panel.design_length ({panel.design_length}) after trimming. "
                                 "If that is a lot, check panel.design_length and panel.adapter.")
        if not any(length < panel.sequencing_length for length in lengths):
            self.warnings.append("No kept allele was shortened by trimming, so no adapter was found in them. "
                                 "Check panel.design_length and panel.adapter.")

        call_step(update_trimmed, [trimmed, self.db, self.ids_report])
        self.ids_report = f"{self.base}_tmp_rename_updatedSeq.csv"
        self.blast_query = f"{self.base}_match_cutadapt_unique.fa"
        self.final_report = f"{self.base}_rename_updatedSeq.csv"
        merged = after - count_records(self.blast_query)
        self.summary["alleles_dropped_by_trimming"], self.summary["alleles_merged_after_trimming"] = dropped, merged
        return f"{after} of {before} alleles kept, {merged} merged as duplicates"

    def stage_blast(self):
        run_tool(["makeblastdb", "-in", self.db, "-dbtype", "nucl"])
        self.blast_out = f"{self.blast_query}.alleledb.bn"
        run_tool(["blastn", "-task", "blastn-short", "-dust", "no", "-soft_masking", "false",
                  "-db", self.db, "-query", self.blast_query, "-out", self.blast_out,
                  "-evalue", "1e-5", "-num_threads", str(self.cfg.blast.threads),
                  "-max_target_seqs", "15", "-outfmt", BLAST_COLUMNS])
        with_hit = set()
        with open(self.blast_out) as handle:
            for line in handle:
                cells = line.split(None, 5)
                if cells[0].split("|")[0] == cells[4].split("|")[0]:
                    with_hit.add(cells[0])
        queried = count_records(self.blast_query)
        self.summary["alleles_blasted"] = queried
        if not with_hit:
            self.warnings.append("No allele had a BLAST hit at its own marker. Check that inputs.allele_db was "
                                 "built for this panel.")
        return f"{len(with_hit)} of {queried} alleles hit their own marker"

    def stage_assign_ids(self):
        cfg = self.cfg
        call_step(assign_ids, [
            self.counts, self.db, self.ids_report, self.blast_out, str(cfg.panel.sequencing_length),
            "--read_col", self.first_sample_col,
            "--cov_threshold", str(cfg.blast.min_coverage),
            "--iden_threshold", str(cfg.blast.min_identity),
            "--allele-type-ground-truth", cfg.blast.allele_type_truth,
        ])
        self.madc_cleaned = self.final_report
        with open(self.final_report) as handle:
            header = handle.readline().strip().split(",")
        samples = len(header) - header.index("AlleleSequence") - 1
        if samples != self.n_samples:
            raise PipelineError(f"The report has {samples} sample columns but the MADC has {self.n_samples}. "
                                f"Its columns start {header[3:6]}, so the first sample column may be wrong.")

        with open(self.final_report) as handle:
            next(handle)
            assigned = sum(1 for line in handle
                           if not line.split(",", 1)[0].endswith(("|Ref_0001", "|Alt_0002")))
        new_db = assign_ids.bump_version_tag(self.db)
        new = len(read_fasta_ids(new_db) - read_fasta_ids(self.db)) if os.path.exists(new_db) else 0
        existing = assigned - new
        removed = self.summary["alleles_blasted"] - assigned
        self.summary.update(alleles_existing=existing, alleles_new=new, alleles_removed_after_blast=removed)
        note = f"{existing} already in the database, {new} new, {removed} removed by the BLAST thresholds or merged"
        other_type = f"{self.db[:-3]}_mismatchedType_alleles.csv"
        if os.path.exists(other_type):
            with open(other_type) as handle:
                closer = sum(1 for _ in handle) - 1
            if closer > 0:
                self.summary["alleles_closer_to_other_type"] = closer
                note += f"; {closer} of the new ones are closest to the other allele type"
        return note

    def stage_db_check(self):
        new_db = assign_ids.bump_version_tag(self.db)
        new_counts = assign_ids.bump_version_tag(self.counts)
        if not os.path.exists(new_db):
            print("No new alleles found, so no new database was written.")
            self.readme = self.no_new_allele_readme
            self._add_code_version()
            self._tidy()
            return f"no new alleles, database stays at {len(read_fasta_ids(self.db))} sequences"

        missing = minted_ids_missing_from_fasta(self.counts, new_counts, new_db)
        if missing:
            raise PipelineError(f"{len(missing)} new allele ID(s) were issued with no sequence in {new_db}: "
                                f"{missing[:5]}. The database would be inconsistent, so the run was stopped.")

        run_tool(["makeblastdb", "-in", new_db, "-dbtype", "nucl"])
        call_step(check_db, [new_db])
        self.final_db, self.final_counts = new_db, new_counts
        duplicates = f"{new_db}.dup.csv"
        if os.path.exists(duplicates):
            self.final_db = assign_ids.bump_version_tag(new_db)
            self.final_counts = f"{self.final_db[:-3]}_matchCnt_lut.txt"
            run_tool(["makeblastdb", "-in", self.final_db, "-dbtype", "nucl"])
            call_step(update_madc_ids, [duplicates, self.madc_cleaned])
            without_duplicates = f"{self.madc_cleaned[:-4]}_rmDup.csv"
            if os.path.exists(without_duplicates):
                self.madc_cleaned = without_duplicates
        self._add_code_version()
        self._tidy()
        note = f"{len(read_fasta_ids(self.final_db))} sequences, was {len(read_fasta_ids(self.db))}"
        if self.final_db != new_db:
            collapsed = len(read_fasta_ids(new_db)) - len(read_fasta_ids(self.final_db))
            self.warnings.append(f"{collapsed} duplicate sequence(s) in the new database were collapsed, "
                                 f"which made another version: {self.final_db}")
        return note

    def _add_code_version(self):
        version = self.cfg.output.code_version
        self.final_report = f"{self.madc_cleaned[:-4]}_{version}.csv"
        with open(self.madc_cleaned) as source, open(self.final_report, "w") as target:
            target.write(f"Code_version,{source.readline()}")
            for line in source:
                target.write(f"{version},{line}")

    def _tidy(self):
        Path("report/matches").mkdir(exist_ok=True)
        for pattern in ("report/*match.fa*", "report/*_match_*"):
            for path in glob.glob(pattern):
                shutil.move(path, "report/matches/")
        with contextlib.suppress(FileNotFoundError):
            os.remove(f"{self.base}_tmp_rename.csv")
        with contextlib.suppress(FileNotFoundError):
            os.remove(f"{self.base}_tmp_rename_updatedSeq.csv")

    def stage_hap_pa(self):
        cfg = self.cfg
        call_step(hap_pa, [
            "--madc", self.madc_cleaned,
            "--descriptor", cfg.hap_pa.descriptor or cfg.project_id,
            "--read_threshold", str(cfg.hap_pa.min_reads),
        ])
        status = f"{self.madc_cleaned[:-4]}_hapStatus.csv"
        project_matrix = f"{self.madc_cleaned[:-4]}_hapPA.csv"
        Path("hapSearch/hapPA").mkdir(parents=True)
        Path("hapSearch/sample_hapStatus").mkdir(parents=True)
        self.hap_status = f"hapSearch/sample_hapStatus/{os.path.basename(status)}"
        shutil.copyfile(status, self.hap_status)
        db_name = os.path.basename(self.final_db)[:-3]
        self.hap_pa_matrix = f"hapSearch/hapPA/{db_name}_{cfg.project_id}_hapPA.csv"
        if cfg.inputs.previous_hap_pa:
            previous = f"hapSearch/hapPA/previous_{cfg.inputs.previous_hap_pa.name}"
            shutil.copyfile(cfg.inputs.previous_hap_pa, previous)
            call_step(merge_pa, ["--reports", f"{previous},{project_matrix}", "--output", self.hap_pa_matrix])
        else:
            shutil.copyfile(project_matrix, self.hap_pa_matrix)
        with open(project_matrix) as handle:
            present = sum(1 for _ in handle) - 1
        return f"{present} alleles present in at least one sample"

    def outputs(self):
        return {
            "report": self.final_report,
            "allele_db": self.final_db,
            "match_counts": self.final_counts,
            "hap_pa": self.hap_pa_matrix,
            "hap_status": self.hap_status,
        }

    def counts_summary(self):
        summary = {"database_before": len(read_fasta_ids(self.db)), **self.summary}
        if os.path.exists(self.final_db):
            summary["database_after"] = len(read_fasta_ids(self.final_db))
        if self.final_report and os.path.exists(self.final_report):
            with open(self.final_report) as handle:
                header = handle.readline().strip().split(",")
                summary["report_rows"] = sum(1 for _ in handle)
            summary["report_samples"] = len(header) - header.index("AlleleSequence") - 1
        return summary


def run_assign(cfg, until=None, verbose=False, command=None, console=None):
    """Run the workflow up to and including the stage named by until. Returns the saved manifest."""
    findings, layout = preflight(cfg)
    problems = [finding.message for finding in findings if finding.level == PROBLEM]
    if problems:
        raise CannotStart("Cannot start:\n  " + "\n  ".join(problems))
    first_sample_col = cfg.panel.first_sample_column or layout.first_sample_column
    n_samples = len(layout.samples) + layout.first_sample_column - first_sample_col
    details = {"first_sample_column": first_sample_col, "samples": n_samples}
    warnings = [finding.message for finding in findings if finding.level == WARNING]
    return run_stages(cfg, _AssignRun(cfg, first_sample_col, n_samples), STAGES, SCHEMA, required_tools(cfg), details,
                      warnings, until=until, verbose=verbose, command=command, console=console)
