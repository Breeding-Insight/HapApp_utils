"""The MADC workflow: from a report to fixed allele IDs and a new database version.

This calls the step modules in order, as workflow/20261020-build02_madc_matchhaps_hapsearch.bash
does. Each run works on its own copies of the inputs inside a new folder, and the steps run from
inside that folder so the paths they rewrite never include the user's own directory names.
"""

from __future__ import annotations

import contextlib
import glob
import hashlib
import json
import os
import platform
import re
import shutil
import sys
import time
from datetime import datetime
from pathlib import Path

import yaml

from .. import __version__
from ..checks import PROBLEM, minted_ids_missing_from_fasta, preflight, read_fasta_ids
from ..config import config_to_dict
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
from ..tools import required_tools, run_tool, tool_version

STAGES = {
    "update-ids": "Update marker IDs and check for duplicate alleles",
    "filter": "Filter alleles on read counts",
    "trim": "Remove adapters",
    "blast": "BLAST alleles against the database",
    "assign-ids": "Assign fixed allele IDs",
    "db-check": "Check the new database",
    "hap-pa": "Build presence/absence tables",
}
BLAST_COLUMNS = "6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue"


class PipelineError(Exception):
    """The run cannot continue, for a reason worded for the person running the tool."""


class DuplicateAllelesFound(PipelineError):
    """The MADC has the same sequence under different allele IDs."""


class _Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, text):
        for stream in self.streams:
            stream.write(text)
        return len(text)

    def flush(self):
        for stream in self.streams:
            stream.flush()


@contextlib.contextmanager
def _working_directory(path):
    previous = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _call_step(step, argv):
    """Run a step's main(), turning its argument errors into a normal failure."""
    try:
        step.main(argv)
    except SystemExit as err:
        if err.code not in (0, None):
            raise PipelineError(f"{step.__name__} rejected its arguments: {argv}") from err


class _AssignRun:
    """One run. Paths are relative to the run folder, as the steps expect to rewrite them."""

    def __init__(self, cfg, first_sample_col):
        self.cfg = cfg
        self.first_sample_col = str(first_sample_col)
        self.report = f"report/{cfg.inputs.madc.stem}.csv"
        self.lut = f"lut/{cfg.inputs.marker_lut.name}"
        self.db = f"versions/{cfg.inputs.allele_db.stem}.fa"
        # step06 expects the count table to be named after the database
        self.counts = f"{self.db[:-3]}_matchCnt_lut.txt"
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

    def copy_inputs(self):
        for folder in ("report", "versions", "lut"):
            Path(folder).mkdir()
        shutil.copyfile(self.cfg.inputs.madc, self.report)
        shutil.copyfile(self.cfg.inputs.marker_lut, self.lut)
        shutil.copyfile(self.cfg.inputs.allele_db, self.db)
        shutil.copyfile(self.cfg.inputs.match_counts, self.counts)

    def stage_update_ids(self):
        cfg = self.cfg
        _call_step(update_ids, [self.lut, self.report])
        stem = self.report[:-4]
        self.snpid = f"{stem}_snpID.csv"
        if cfg.inputs.dup_tags and cfg.inputs.dup_tags.stat().st_size > 0:
            tags = f"lut/{cfg.inputs.dup_tags.name}"
            shutil.copyfile(cfg.inputs.dup_tags, tags)
            without_tags = f"{stem}_snpID_rmDup.csv"
            _call_step(remove_tags, [tags, self.snpid, without_tags])
            self.snpid = without_tags
        self.base = self.snpid[:-4]
        _call_step(check_unique, [self.snpid, self.first_sample_col])
        duplicates = f"{self.base}_dup.csv"
        if os.path.exists(duplicates):
            raise DuplicateAllelesFound(
                f"The MADC has the same sequence under different allele IDs. Review {duplicates} in the run folder.")

    def stage_filter(self):
        cfg = self.cfg
        _call_step(filter_alleles, [
            self.snpid, "--lut", self.lut,
            "--min_samples", str(cfg.filter.min_samples),
            "--min_reads", str(cfg.filter.min_reads),
            "--first_sample_col", self.first_sample_col,
            "--split_by_colID", cfg.panel.insert_size_column or "",
        ])
        self.ids_report = f"{self.base}_tmp_rename.csv"
        self.blast_query = f"{self.base}_match.fa"
        self.final_report = f"{self.base}_rename.csv"

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
        _call_step(update_trimmed, [trimmed, self.db, self.ids_report])
        self.ids_report = f"{self.base}_tmp_rename_updatedSeq.csv"
        self.blast_query = f"{self.base}_match_cutadapt_unique.fa"
        self.final_report = f"{self.base}_rename_updatedSeq.csv"
        return None

    def stage_blast(self):
        run_tool(["makeblastdb", "-in", self.db, "-dbtype", "nucl"])
        self.blast_out = f"{self.blast_query}.alleledb.bn"
        run_tool(["blastn", "-task", "blastn-short", "-dust", "no", "-soft_masking", "false",
                  "-db", self.db, "-query", self.blast_query, "-out", self.blast_out,
                  "-evalue", "1e-5", "-num_threads", str(self.cfg.blast.threads),
                  "-max_target_seqs", "15", "-outfmt", BLAST_COLUMNS])

    def stage_assign_ids(self):
        cfg = self.cfg
        _call_step(assign_ids, [
            self.counts, self.db, self.ids_report, self.blast_out, str(cfg.panel.sequencing_length),
            "--read_col", self.first_sample_col,
            "--cov_threshold", str(cfg.blast.min_coverage),
            "--iden_threshold", str(cfg.blast.min_identity),
            "--allele-type-ground-truth", cfg.blast.allele_type_truth,
        ])
        self.madc_cleaned = self.final_report

    def stage_db_check(self):
        new_db = assign_ids.bump_version_tag(self.db)
        new_counts = assign_ids.bump_version_tag(self.counts)
        if not os.path.exists(new_db):
            print("No new alleles found, so no new database was written.")
            self._add_code_version()
            self._tidy()
            return "no new alleles"

        missing = minted_ids_missing_from_fasta(self.counts, new_counts, new_db)
        if missing:
            raise PipelineError(f"{len(missing)} new allele ID(s) were issued with no sequence in {new_db}: "
                                f"{missing[:5]}. The database would be inconsistent, so the run was stopped.")

        run_tool(["makeblastdb", "-in", new_db, "-dbtype", "nucl"])
        _call_step(check_db, [new_db])
        self.final_db, self.final_counts = new_db, new_counts
        duplicates = f"{new_db}.dup.csv"
        if os.path.exists(duplicates):
            self.final_db = assign_ids.bump_version_tag(new_db)
            self.final_counts = f"{self.final_db[:-3]}_matchCnt_lut.txt"
            run_tool(["makeblastdb", "-in", self.final_db, "-dbtype", "nucl"])
            _call_step(update_madc_ids, [duplicates, self.madc_cleaned])
            without_duplicates = f"{self.madc_cleaned[:-4]}_rmDup.csv"
            if os.path.exists(without_duplicates):
                self.madc_cleaned = without_duplicates
        self._add_code_version()
        self._tidy()
        return None

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
        _call_step(hap_pa, [
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
            _call_step(merge_pa, ["--reports", f"{previous},{project_matrix}", "--output", self.hap_pa_matrix])
        else:
            shutil.copyfile(project_matrix, self.hap_pa_matrix)

    def outputs(self):
        return {
            "report": self.final_report,
            "allele_db": self.final_db,
            "match_counts": self.final_counts,
            "hap_pa": self.hap_pa_matrix,
            "hap_status": self.hap_status,
        }

    def counts_summary(self):
        summary = {"database_before": len(read_fasta_ids(self.db))}
        if os.path.exists(self.final_db):
            summary["database_after"] = len(read_fasta_ids(self.final_db))
        if self.final_report and os.path.exists(self.final_report):
            with open(self.final_report) as handle:
                header = handle.readline().strip().split(",")
                summary["report_rows"] = sum(1 for _ in handle)
            summary["report_samples"] = len(header) - header.index("AlleleSequence") - 1
        return summary


def run_assign(cfg, until=None, verbose=False, command=None, console=None):
    """Run the workflow up to and including the stage named by until. Returns a summary."""
    console = console or sys.stdout
    findings, layout = preflight(cfg)
    problems = [finding.message for finding in findings if finding.level == PROBLEM]
    if problems:
        raise PipelineError("Cannot start:\n  " + "\n  ".join(problems))

    outdir = cfg.output.dir
    outdir.mkdir(parents=True, exist_ok=True)
    first_sample_col = cfg.panel.first_sample_column or layout.first_sample_column
    (outdir / "resolved_config.yaml").write_text(yaml.safe_dump(config_to_dict(cfg), sort_keys=False))
    manifest = {
        "hapapp_utils_version": __version__,
        "python": platform.python_version(),
        "command": command,
        "started": datetime.now().isoformat(timespec="seconds"),
        "status": "running",
        "first_sample_column": first_sample_col,
        "samples": len(layout.samples),
        "tools": {tool: tool_version(tool) for tool in required_tools(cfg)},
        "inputs": {name: {"path": str(path), "sha256": _sha256(path)}
                   for name, path in vars(cfg.inputs).items() if path is not None},
    }
    manifest_path = outdir / "run_manifest.json"

    def save_manifest():
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    save_manifest()
    run = _AssignRun(cfg, first_sample_col)
    with open(outdir / "pipeline.log", "w", buffering=1) as log:
        stream = _Tee(log, console) if verbose else log

        def say(text):
            print(text, file=console, flush=True)
            print(text, file=log)

        try:
            with _working_directory(outdir), contextlib.redirect_stdout(stream), contextlib.redirect_stderr(stream):
                run.copy_inputs()
                names = list(STAGES)
                for number, name in enumerate(names, 1):
                    say(f"[{number}/{len(names)}] {STAGES[name]}")
                    start = time.monotonic()
                    note = getattr(run, f"stage_{name.replace('-', '_')}")()
                    say(f"      {note or 'done'} ({time.monotonic() - start:.0f} s)")
                    if name == until:
                        break
                manifest["status"] = "finished" if until in (None, names[-1]) else f"stopped after {until}"
                manifest["outputs"] = {key: value for key, value in run.outputs().items()
                                       if value and os.path.exists(value)}
                manifest["counts"] = run.counts_summary()
        except Exception as err:
            manifest["status"] = "failed"
            manifest["error"] = str(err)
            print(f"[ERROR] {err}", file=log)
            raise
        finally:
            manifest["finished"] = datetime.now().isoformat(timespec="seconds")
            save_manifest()
    return manifest
