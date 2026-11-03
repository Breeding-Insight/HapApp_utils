"""The core database workflow: Ref and Alt alleles for a panel, built from the reference genome.

This calls the step modules in the same order as workflow/20260320-build01_species_refAlt_coreDB.bash.
The MADC cannot supply these alleles directly because its Ref and Alt sequences can contain IUPAC
codes, so it is only used to find where each amplicon sits on the genome.
"""

from __future__ import annotations

import os
import shutil
from collections import Counter
from pathlib import Path

from ..checks import PROBLEM, preflight_build_db, read_chromosome_map
from ..config import BUILD_DB_SCHEMA
from ..fasta import fetch_subsequences, iter_fasta, read_lengths, write_fasta
from ..madc import read_layout, read_marker_ids
from ..probe_design import read_probe_rows, write_probe_csv
from ..steps.core_db import (
    db00_prep_lut_from_probeDesign as prepare_lookup,
    db01_get_reference_sfetch_keys_from_snpID_lut_v1 as flank_keys,
    db01_prep_ref_alt_flankSeq_from_lut_v1 as ref_alt_flanks,
    db02_update_snpID_in_madc_v1 as update_ids,
    db03_ext_ref_alt_amp_from_madc_v1 as extract_amplicons,
    db05_determine_alleleOri_from_blast_AND_update_f180bp_v1 as orient_flanks,
    db07_generate_ref_alt_sfetch_keys_from_blast_v1_1 as allele_keys,
    db09_rm_dupTags_from_LUT_and_db_v001 as remove_markers,
)
from ..tools import BUILD_DB_TOOLS, run_tool
from .common import BLAST_COLUMNS, PipelineError, call_step, run_stages

STAGES = {
    "lookup": "Build the marker lookup table",
    "flanks": "Fetch flanking sequences from the genome",
    "orient": "Find the strand of each amplicon",
    "alleles": "Fetch Ref and Alt alleles",
    "duplicates": "Remove markers with near-identical Ref alleles",
    "finish": "Write the new database",
}


class _BuildDbRun:
    """One run. Paths are relative to the run folder, as the steps expect to rewrite them."""

    def __init__(self, cfg):
        self.cfg = cfg
        self.flank = cfg.database.flank_length
        self.length = cfg.database.allele_length
        self.report = f"report/{cfg.inputs.madc.stem}.csv"
        self.base = self.report[:-4]
        self.probe = f"probe/{cfg.inputs.probe_design.stem}.csv"
        self.chromosome_map = read_chromosome_map(cfg.inputs.chromosome_map)
        self.lookup = f"{self.probe[:-4]}_snpID_lut.csv"
        self.snpid = f"{self.base}_snpID.csv"
        self.counts = f"{self.base}_snpID_matchCnt_lut.txt"
        self.amplicons = f"{self.base}_snpID_ref_alt_amplicons.fa"
        self.ref_alt_flanks = f"{self.lookup[:-4]}_f{self.flank}bp_sfetchKeys_ref_alt.fa"
        self.oriented_flanks = f"{self.ref_alt_flanks[:-3]}_rev.fa"
        self.bottom_strand = f"{self.ref_alt_flanks[:-3]}.botloci"
        self.alleles = f"{self.amplicons}.f{self.flank}bp_rev_{self.length}bp_sfetch.fa"
        self.removed = f"{self.base}_markers_removed.tsv"
        self.final_lookup = None
        self.final_db = f"versions/{cfg.database.name}_v001.fa"
        self.final_counts = f"versions/{cfg.database.name}_v001_matchCnt_lut.txt"
        self.summary = {}

    def copy_inputs(self):
        for folder in ("probe", "report", "versions"):
            Path(folder).mkdir()
        header, rows = read_probe_rows(self.cfg.inputs.probe_design)
        write_probe_csv(header, rows, self.probe)
        shutil.copyfile(self.cfg.inputs.madc, self.report)

    def stage_lookup(self):
        call_step(prepare_lookup, [self.probe, "--madc", self.report])
        call_step(update_ids, [self.lookup, self.report])
        with open(self.lookup) as handle:
            markers = sum(1 for _ in handle) - 1
        if markers == 0:
            raise PipelineError("No marker is in both the probe design file and the MADC.")
        self.summary["markers_in_lookup"] = markers
        return f"{markers} markers"

    def stage_flanks(self):
        lengths = read_lengths(self.cfg.inputs.genome)
        names_by_genome = {genome: probe for probe, genome in self.chromosome_map.items()}
        lengths_file = "probe/genome_lengths.txt"
        with open(lengths_file, "w") as out:
            for name, length in lengths.items():
                out.write(f"{names_by_genome.get(name, name)}\t{length}\n")
        call_step(flank_keys, ["--lut", self.lookup, "--chr_len", lengths_file,
                               "--flankBP", str(self.flank), "--madc", self.snpid])
        keys = f"{self.lookup[:-4]}_f{self.flank}bp_sfetchKeys.txt"
        flanks = f"{keys[:-4]}.fa"
        try:
            fetch_subsequences(self.cfg.inputs.genome, keys, flanks, rename=self.chromosome_map)
        except ValueError as err:
            raise PipelineError(f"{err}. Set inputs.chromosome_map if the genome uses different names.") from err
        call_step(ref_alt_flanks, ["--snpID_lut", self.lookup, "--flankSeq", flanks, "--flank_len", str(self.flank)])
        built = sum(1 for _ in iter_fasta(self.ref_alt_flanks)) // 2
        if built == 0:
            raise PipelineError("No marker's Ref or Alt base matches the genome at its position. "
                                "Check that the genome is the assembly the probes were designed on, "
                                "and that its chromosome names are matched correctly.")
        return f"{built} of {self.summary['markers_in_lookup']} markers match the genome"

    def _blast(self, database, out):
        run_tool(["blastn", "-task", "blastn-short", "-dust", "no", "-soft_masking", "false",
                  "-db", database, "-query", self.amplicons, "-out", out,
                  "-evalue", "1e-5", "-num_threads", str(self.cfg.database.threads),
                  "-max_target_seqs", "10", "-outfmt", BLAST_COLUMNS])

    def stage_orient(self):
        call_step(extract_amplicons, [self.snpid])
        run_tool(["makeblastdb", "-in", self.ref_alt_flanks, "-dbtype", "nucl"])
        hits = f"{self.amplicons}.f{self.flank}bp.bn"
        self._blast(self.ref_alt_flanks, hits)
        call_step(orient_flanks, [hits, self.ref_alt_flanks])
        if not any(True for _ in iter_fasta(self.oriented_flanks)):
            raise PipelineError("No MADC amplicon aligns to its flanking sequence. Check that the genome is "
                                "the assembly the probes were designed on, and that its chromosome names "
                                "are matched correctly.")
        with open(self.bottom_strand) as handle:
            minus = sum(1 for _ in handle)
        self.summary["markers_on_minus_strand"] = minus
        return f"{minus} markers on the minus strand"

    def stage_alleles(self):
        run_tool(["makeblastdb", "-in", self.oriented_flanks, "-dbtype", "nucl"])
        hits = f"{self.amplicons}.f{self.flank}bp_rev.bn"
        self._blast(self.oriented_flanks, hits)
        call_step(allele_keys, [self.lookup, hits, str(self.length)])
        keys = f"{hits[:-3]}_{self.length}bp_sfetchKeys.txt"
        return f"{fetch_subsequences(self.oriented_flanks, keys, self.alleles)} alleles"

    def _duplicate_pairs(self):
        """Pairs of Ref alleles from different markers that are near-identical, from MMseqs2."""
        database = self.cfg.database
        ref_only = f"{self.alleles[:-3]}_ref_only.fa"
        write_fasta(((name, sequence) for name, sequence in iter_fasta(self.alleles)
                     if name.endswith("|Ref_0001")), ref_only)
        prefix = f"{ref_only[:-3]}_dupTag"
        temporary = "report/mmseqs_tmp"
        run_tool(["mmseqs", "easy-linclust", ref_only, prefix, temporary,
                  "--min-seq-id", str(database.dup_min_identity), "-c", str(database.dup_min_coverage),
                  "--cov-mode", "2", "-v", "0", "--threads", str(database.threads)])
        with open(f"{prefix}_cluster.tsv") as handle:
            pairs = [line.split() for line in handle if len(set(line.split())) == 2]
        shutil.rmtree(temporary)
        for leftover in (ref_only, f"{prefix}_all_seqs.fasta", f"{prefix}_rep_seq.fasta"):
            os.remove(leftover)
        return pairs

    def _markers_missing_alleles(self):
        """Markers in the lookup table without both alleles, which assign could not use."""
        names = {name for name, _ in iter_fasta(self.alleles)}
        with open(self.lookup) as handle:
            markers = [line.split(",")[1] for line in handle][1:]
        return [marker for marker in markers
                if f"{marker}|Ref_0001" not in names or f"{marker}|Alt_0002" not in names]

    def stage_duplicates(self):
        pairs = self._duplicate_pairs()
        missing = self._markers_missing_alleles()
        duplicates = sorted({pair[1].split("|")[0] for pair in pairs})
        self.summary["markers_removed_as_duplicates"] = duplicates
        self.summary["markers_removed_without_alleles"] = missing
        self.final_lookup = self.lookup
        self.current_counts, self.current_alleles = self.counts, self.alleles
        if not pairs and not missing:
            return "none removed"

        # db09 removes the marker named in the second column of each line
        with open(self.removed, "w") as out:
            for keep, remove in pairs:
                out.write(f"{keep}\t{remove}\n")
            for marker in missing:
                out.write(f"no_alleles\t{marker}\n")
        for path in (self.lookup, self.counts, self.alleles):
            call_step(remove_markers, [path, self.removed])
        self.final_lookup = self.lookup.replace("snpID_lut.csv", "rmDupTag_snpID_lut.csv")
        self.current_counts = self.counts.replace("matchCnt_lut.txt", "matchCnt_lut_rmDupTag.txt")
        self.current_alleles = self.alleles.replace("sfetch.fa", "sfetch_rmDupTag.fa")
        return f"{len(duplicates)} duplicate(s), {len(missing)} without both alleles"

    def stage_finish(self):
        shutil.copyfile(self.current_alleles, self.final_db)
        shutil.copyfile(self.current_counts, self.final_counts)
        run_tool(["makeblastdb", "-in", self.final_db, "-dbtype", "nucl"])
        alleles = dict(iter_fasta(self.final_db))
        if not alleles:
            raise PipelineError("No alleles were left to write to the new database.")
        self.summary["alleles"] = len(alleles)
        self.summary["markers_in_database"] = len({name.split("|")[0] for name in alleles})
        lengths = Counter(len(sequence) for sequence in alleles.values())
        self.summary["allele_lengths"] = {str(length): count for length, count in sorted(lengths.items())}

        # An allele can only fall short by running off the end of its flanking sequence
        with open(self.final_lookup) as handle:
            snps = {cells[1] for cells in (line.split(",") for line in handle) if cells[6].strip().upper() == "SNP"}
        short = sorted(name for name, sequence in alleles.items()
                       if name.split("|")[0] in snps and len(sequence) < self.length)
        self.summary["short_snp_alleles"] = short
        return f"{len(short)} SNP alleles shorter than {self.length} bp; raise flank_length" if short else None

    def outputs(self):
        return {
            "marker_lut": self.final_lookup,
            "allele_db": self.final_db,
            "match_counts": self.final_counts,
            "bottom_strand_markers": self.bottom_strand,
            "markers_removed": self.removed,
        }

    def counts_summary(self):
        return self.summary


def run_build_db(cfg, until=None, verbose=False, command=None, console=None):
    """Run the workflow up to and including the stage named by until. Returns the saved manifest."""
    findings = preflight_build_db(cfg)
    problems = [finding.message for finding in findings if finding.level == PROBLEM]
    if problems:
        raise PipelineError("Cannot start:\n  " + "\n  ".join(problems))
    details = {
        "madc_samples": len(read_layout(cfg.inputs.madc).samples),
        "madc_markers": len(read_marker_ids(cfg.inputs.madc)),
    }
    return run_stages(cfg, _BuildDbRun(cfg), STAGES, BUILD_DB_SCHEMA, BUILD_DB_TOOLS, details,
                      until=until, verbose=verbose, command=command, console=console)
