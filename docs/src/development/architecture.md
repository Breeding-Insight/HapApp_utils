# Architecture

## Layout

```text
src/hapapp_utils/
  cli.py                 the command line
  config.py              settings: schemas, templates, loading and validation
  checks.py              checks made before a run and after a database update
  madc.py                reads a MADC's layout without its counts
  probe_design.py        reads the probe file from Excel, CSV or text
  fasta.py               reads FASTA files and fetches subsequences
  tools.py               finds and runs external programs
  workflows/
    common.py            the run folder, process readme and manifest
    assign.py            the assign stages
    build_db.py          the build-db stages
  steps/
    assign/              steps from scripts/RefMatch_AltMatch_Other
    core_db/             steps from scripts/refAlt_coreDB
    hap_pa/              steps from scripts/hapPA
tests/                   the test suite
workflow/                the original bash workflows
docs/                    this documentation
```

The workflows decide the order of the steps and hold no analysis of their own. The analysis lives in the step
modules, and each module keeps the functions and arguments of its original script.

## Stages and steps: assign

| Stage | Steps |
|-------|-------|
| `update-ids` | `step00_madc_update_snpID_v1_1`, `step00_rm_dupTag_closeTag_rawMADC_v1`, `step00_check_allele_uniqueness` |
| `filter` | `step01_filter_missing_AND_ext_matchAlleles_from_madc_v1` |
| `trim` | cutadapt, `step03_check_cutadapt_allele_uniqueness_AND_update_tmp_rename_report_v1_1` |
| `blast` | `makeblastdb`, `blastn` |
| `assign-ids` | `step05_parse_madc_matchAlleles_blastn` |
| `db-check` | `step06_check_db_allele_uniqueness_v1`, `step06_update_MADC_with_allele_uniqueness_v1` |
| `hap-pa` | `util_haplotype_presence_absence_matrix`, `util_hapDB_merge_presence_absence_PIhaps` |

## Stages and steps: build-db

| Stage | Steps |
|-------|-------|
| `lookup` | `db00_prep_lut_from_probeDesign`, `db02_update_snpID_in_madc_v1` |
| `flanks` | `db01_get_reference_sfetch_keys_from_snpID_lut_v1`, `fasta.fetch_subsequences`, `db01_prep_ref_alt_flankSeq_from_lut_v1` |
| `orient` | `db03_ext_ref_alt_amp_from_madc_v1`, `makeblastdb`, `blastn`, `db05_determine_alleleOri_from_blast_AND_update_f180bp_v1` |
| `alleles` | `makeblastdb`, `blastn`, `db07_generate_ref_alt_sfetch_keys_from_blast_v1_1`, `fasta.fetch_subsequences` |
| `duplicates` | `mmseqs easy-linclust`, `db09_rm_dupTags_from_LUT_and_db_v001` |
| `finish` | `makeblastdb` |

## How a run works

`workflows/common.py` runs every workflow the same way:

1. It creates the output folder, and saves the settings and the input checksums.
2. It opens the process readme and writes the inputs and settings at the top.
3. It moves into the output folder, so the paths the steps rewrite never contain your own folder names.
4. It runs each stage and shows the stage's note and any new warnings.
5. It writes the outputs and counts to `run_manifest.json`, whether the run finished or not.

Each stage is a `stage_<name>` method on a run class. A step runs through its `main(argv)` function, the same way
the command line would call it. The workflow sends the step's printed output to the process readme.

## External programs

| Program | Used in | Why it stays external |
|---------|---------|----------------------|
| BLAST+ | both | Alignment |
| MMseqs2 | `build-db` | Clustering |
| cutadapt | `assign` | Adapter trimming. It installs with the package |

`fasta.py` replaces `esl-sfetch` and `seqkit`, since those two only cut ranges and pick records by name.
