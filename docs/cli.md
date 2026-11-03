# hapapp-utils

## Install

```bash
pip install -e ".[test]"
```

This installs the Python dependencies, including cutadapt. BLAST+ (`blastn`, `makeblastdb`) must also be on your PATH, and `build-db` also needs MMseqs2 (`mmseqs`). `hapapp-utils check` tells you if anything is missing.

## Commands

| Command | What it does |
|---|---|
| `init [file]` | Write a commented settings file, for `assign` or with `--workflow build-db` |
| `check -c file` | Check tools, inputs and settings without running anything |
| `build-db -c file` | Build a panel's first allele database from the reference genome |
| `assign -c file` | Turn a MADC into fixed allele IDs and a new database version |
| `hap-pa` | Build presence/absence tables from a report with fixed IDs |
| `merge-pa` | Merge presence/absence matrices from several projects |

`build-db`, `assign` and `check` also accept:

- `--set key=value` to override a setting, for example `--set filter.min_reads=30`. It can be repeated.
- `-o folder` to choose the output folder.

`build-db` and `assign` also accept:

- `--until stage` to stop after a stage.
  - For `build-db`: `lookup`, `flanks`, `orient`, `alleles`, `duplicates` or `finish`.
  - For `assign`: `update-ids`, `filter`, `trim`, `blast`, `assign-ids`, `db-check` or `hap-pa`.
- `-v` to show each step's output as well as saving it to the log.

To try different thresholds, run again with `--set` and a new `-o` folder. Earlier results and the input database are never changed.

Exit codes: 0 finished, 1 error, 2 bad command-line arguments, 3 duplicate alleles found in the MADC.

## Settings

Settings on the command line replace those in the file. Relative paths in the file are resolved from the file's own folder, and those on the command line from the current folder.

### build-db

| Setting | Variable in the bash workflow | Required | Default |
|---|---|---|---|
| `inputs.probe_design` | `PROBE_FILE` | yes | |
| `inputs.madc` | `REPORT` | yes | |
| `inputs.genome` | `REF_GENOME` | yes | |
| `inputs.chromosome_map` | none | no | none |
| `database.name` | `ALLELE_DB`, `MATCH_CNT` | yes | |
| `database.allele_length` | `REF_LEN` | yes | |
| `database.flank_length` | `FLANK_LEN` | no | 150 |
| `database.dup_min_identity` | written in the mmseqs command | no | 0.95 |
| `database.dup_min_coverage` | written in the mmseqs command | no | 0.99 |
| `database.threads` | written in the blastn and mmseqs commands | no | 6 |
| `output.dir` | none | yes | |

`inputs.probe_design` can be the DArT template as `.xlsx`, or `.csv` or tab-delimited `.txt`. Any notes above the `Marker Name` row are skipped. `CHR_LEN` is no longer needed because lengths are read from the genome.

The chromosome names in the probe file must match the genome. If they don't, either use a genome renamed to match, or give `inputs.chromosome_map`: a CSV with no header, one `probe name,genome name` per line. `check` reports any mismatch before a run starts.

### assign

| Setting | Variable in the bash workflow | Required | Default |
|---|---|---|---|
| `project_id` | `PROJECT_ID` | yes | |
| `inputs.madc` | `REPORT` | yes | |
| `inputs.marker_lut` | `MARKERID_LUT` | yes | |
| `inputs.allele_db` | `ALLELE_DB_DIR`, `ALLELE_DB` | yes | |
| `inputs.match_counts` | `MATCHCNT_LUT` | yes | |
| `inputs.dup_tags` | `DUP_TAGS` | no | none |
| `inputs.previous_hap_pa` | `LATEST_HAPPA_MATRIX` | no | none |
| `panel.design_length` | `DESIGN_LEN` | yes | |
| `panel.sequencing_length` | `SEQ_LEN` | yes | |
| `panel.first_sample_column` | `FIRST_SAMPLE_COL` | no | read from the MADC header |
| `panel.insert_size_column` | `LUT_INSERT_SIZE_COL` | no | none |
| `panel.adapter` | written in the cutadapt command | no | `ACCGATCTG` |
| `filter.min_samples` | `MIN_SAMPLES` | yes | |
| `filter.min_reads` | `MIN_READS` | yes | |
| `blast.min_coverage` | `COV` | yes | |
| `blast.min_identity` | `IDEN` | yes | |
| `blast.allele_type_truth` | `ALLELE_TYPE` | no | `report` |
| `blast.threads` | written in the blastn command | no | 6 |
| `hap_pa.min_reads` | `MIN_READS_FOR_PRESENCE` | no | 10 |
| `hap_pa.descriptor` | `PROJECT_DESCRIPTOR` | no | `project_id` |
| `output.dir` | none | yes | |
| `output.code_version` | `CODE_VER` | no | `v2`, the package's major version |

`filter.min_samples` is lowered to 5% of the samples, rounded up, if that is smaller. `check` warns when this happens.

## What a run writes

The output folder must be new or empty. Every run also writes:

- **A process readme** in `versions/`, named as in the bash workflows: `<database>_v004_process.readme` for the version being written, or `<database>_v003_noNewAllele.readme` if there were no new alleles. It opens with the input files, settings and tool versions, followed by everything each step printed.
- **`resolved_config.yaml`**, with the settings used.
- **`run_manifest.json`**, with the status, counts, tool versions, input checksums and output paths.

### build-db

| Path | Contents |
|---|---|
| `versions/*_v001.fa` | Ref and Alt alleles, all uppercase, with one line per sequence |
| `versions/*_v001_matchCnt_lut.txt` | The count table, starting at zero |
| `probe/*_snpID_lut.csv` | The lookup table, or `*_rmDupTag_snpID_lut.csv` if markers were removed |
| `probe/*.botloci` | Markers on the minus strand of the genome |
| `report/*_markers_removed.tsv` | Markers removed, if any: in each line, the second column is the marker removed, and the first is the marker it duplicates or `no_alleles` |

The allele database, count table and lookup table are the `allele_db`, `match_counts` and `marker_lut` inputs for `assign`. Markers can be removed for two reasons: their Ref allele is near-identical to another marker's, or they don't have both alleles because the probe bases don't match the genome or the amplicon doesn't align. They're removed from all three files, since `assign` could not use them.

### assign

It holds copies of the inputs, so every run stands alone.

| Path | Contents |
|---|---|
| `report/*_rename_v2.csv` | The report with fixed allele IDs and a `Code_version` column |
| `versions/*_v00N.fa` | The new allele database |
| `versions/*_v00N_matchCnt_lut.txt` | The new count table |
| `hapSearch/hapPA/` | The presence/absence matrix |
| `hapSearch/sample_hapStatus/` | Presence or absence of each allele in each sample |

## Working on the steps

| Folder | Holds | Used by |
|---|---|---|
| `src/hapapp_utils/steps/assign/` | `scripts/RefMatch_AltMatch_Other` | `assign`, build02 bash workflow |
| `src/hapapp_utils/steps/hap_pa/` | `scripts/hapPA` | `assign`, `hap-pa`, `merge-pa`, build02 bash workflow |
| `src/hapapp_utils/steps/core_db/` | `scripts/refAlt_coreDB` | `build-db`, build01 bash workflow |

Each step is one module with the same functions and arguments as its original script. Only `_v1.1` in filenames became `_v1_1`, because a dot cannot appear in a module name.

There are three ways to run a step, all with the same arguments:

```bash
python src/hapapp_utils/steps/assign/step05_parse_madc_matchAlleles_blastn.py --help
python -m hapapp_utils.steps.assign.step05_parse_madc_matchAlleles_blastn --help
```

or from Python, by calling that module's `main()` with a list of arguments.

With `pip install -e .`, edits take effect straight away. The bash workflows in `workflow/` still work: set `SCRIPTS_DIR` and `HAPPA_SCRIPT_DIR` to the folders above.

To add a step, write its module with a `main(argv=None)` function, then call it from the matching stage in `src/hapapp_utils/workflows/assign.py` or `build_db.py`. To add a whole stage, add a `stage_<name>` method there and list it in `STAGES`, and also in `ASSIGN_STAGES` or `BUILD_DB_STAGES` in `cli.py`.

Run the tests with `pytest`. They use a small synthetic dataset, and workflow tests are skipped if BLAST is not installed. To compare a real run with saved results:

```bash
HAPAPP_BASELINE_DIR=path/to/saved/run HAPAPP_BASELINE_CONFIG=path/to/run.yaml pytest -k real_data
```

## Differences from the bash workflows

### build-db

On the red clover panel, `build-db` and the build01 bash workflow produce the same alleles. It differs on purpose in these cases:

- `esl-sfetch` and `seqkit` aren't needed. Sequences are fetched in Python, and no index file is written next to the genome.
- A range that runs past the end of a sequence gives a shorter allele and a warning. `esl-sfetch` stops at that point, and the bash workflow carries on with what it has.
- The first database always has one line per sequence. Otherwise it would depend on whether markers were removed.
- The lookup table only holds markers in the MADC, because `db00` is given `--madc`. This follows HapApp's build01 workflow.
- A marker without both alleles is removed and reported.
- `flank_length` defaults to 150, as in HapApp's workflow. The bash workflow has 100.

### assign

On a panel with design length equal to sequencing length, `assign` output matches the bash workflow byte for byte, apart from the code version: `v2` in place of the workflow's `v1`. It differs on purpose in these cases:

- Duplicate alleles stop the run with exit code 3 instead of 0.
- Reads longer than the design are trimmed even when there is no insert size column. The bash workflow skipped trimming and continued.
- The insert size column can hold plain numbers such as `81`. The bash workflow only found files named `*_insert*bp.fa`.
- A missing or misspelled `dup_tags` path is an error instead of being skipped.
- A run stops if a new allele ID has no sequence in the new database.
- The report with the `Code_version` column is also written when there are no new alleles.

## What has been tested

`build-db` was checked on the red clover panel, which only has SNPs. Its 5,998 alleles match that panel's existing database, apart from letter case, and it removes the same marker. Indels have only been tested on synthetic data in `tests/`.
