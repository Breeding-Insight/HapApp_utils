# hapapp-utils

## Install

```bash
pip install -e ".[test]"
```

This installs the Python dependencies, including cutadapt. BLAST+ (`blastn`, `makeblastdb`) must also be on your PATH. `hapapp-utils check` tells you if anything is missing.

## Commands

| Command | What it does |
|---|---|
| `init [file]` | Write a commented settings file |
| `check -c file` | Check tools, inputs and settings without running anything |
| `assign -c file` | Turn a MADC into fixed allele IDs and a new database version |
| `hap-pa` | Build presence/absence tables from a report with fixed IDs |
| `merge-pa` | Merge presence/absence matrices from several projects |

`assign` and `check` also accept:

- `--set key=value` to override a setting, for example `--set filter.min_reads=30`. It can be repeated.
- `-o folder` to choose the output folder.

`assign` also accepts:

- `--until stage` to stop after a stage: `update-ids`, `filter`, `trim`, `blast`, `assign-ids`, `db-check` or `hap-pa`.
- `-v` to show each step's output as well as saving it to the log.

To try different thresholds, run again with `--set` and a new `-o` folder. Earlier results and the input database are never changed.

Exit codes: 0 finished, 1 error, 2 bad command-line arguments, 3 duplicate alleles found in the MADC.

## Settings

Settings on the command line replace those in the file. Relative paths in the file are resolved from the file's own folder, and those on the command line from the current folder.

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
| `output.code_version` | `CODE_VER` | no | `v1` |

`filter.min_samples` is lowered to 5% of the samples, rounded up, if that is smaller. `check` warns when this happens.

## What a run writes

The output folder must be new or empty. It holds copies of the inputs, so every run stands alone.

| Path | Contents |
|---|---|
| `report/*_rename_v1.csv` | The report with fixed allele IDs and a `Code_version` column |
| `versions/*_v00N.fa` | The new allele database |
| `versions/*_v00N_matchCnt_lut.txt` | The new count table |
| `hapSearch/hapPA/` | The presence/absence matrix |
| `hapSearch/sample_hapStatus/` | Presence or absence of each allele in each sample |
| `pipeline.log` | Everything each step printed, replacing `_process.readme` |
| `resolved_config.yaml` | The settings actually used |
| `run_manifest.json` | Status, tool versions, input checksums and output paths |

## Working on the steps

| Folder | Holds | Used by |
|---|---|---|
| `src/hapapp_utils/steps/assign/` | `scripts/RefMatch_AltMatch_Other` | `assign`, build02 bash workflow |
| `src/hapapp_utils/steps/hap_pa/` | `scripts/hapPA` | `assign`, `hap-pa`, `merge-pa`, build02 bash workflow |
| `src/hapapp_utils/steps/core_db/` | `scripts/refAlt_coreDB` | build01 bash workflow |

Each step is one module with the same functions and arguments as its original script. Only `_v1.1` in filenames became `_v1_1`, because a dot cannot appear in a module name.

There are three ways to run a step, all with the same arguments:

```bash
python src/hapapp_utils/steps/assign/step05_parse_madc_matchAlleles_blastn.py --help
python -m hapapp_utils.steps.assign.step05_parse_madc_matchAlleles_blastn --help
```

or from Python, by calling that module's `main()` with a list of arguments.

With `pip install -e .`, edits take effect straight away. The bash workflows in `workflow/` still work: set `SCRIPTS_DIR` and `HAPPA_SCRIPT_DIR` to the folders above.

To add a step to `assign`, write its module with a `main(argv=None)` function, then call it from the matching stage in `src/hapapp_utils/workflows/assign.py`. To add a whole stage, add a `stage_<name>` method there and list it in `STAGES`, and also in `ASSIGN_STAGES` in `cli.py`.

Run the tests with `pytest`. They use a small synthetic dataset, and workflow tests are skipped if BLAST is not installed. To compare a real run with saved results:

```bash
HAPAPP_BASELINE_DIR=path/to/saved/run HAPAPP_BASELINE_CONFIG=path/to/run.yaml pytest -k real_data
```

## Differences from the bash workflow

On a panel with design length equal to sequencing length, `assign` output matches the bash workflow byte for byte. It differs on purpose in these cases:

- Duplicate alleles stop the run with exit code 3 instead of 0.
- Reads longer than the design are trimmed even when there is no insert size column. The bash workflow skipped trimming and continued.
- The insert size column can hold plain numbers such as `81`. The bash workflow only found files named `*_insert*bp.fa`.
- A missing or misspelled `dup_tags` path is an error instead of being skipped.
- A run stops if a new allele ID has no sequence in the new database.
- The `_rename_v1.csv` report is also written when there are no new alleles.

## Not yet in the CLI

`build-db` (the build01 workflow) has not been ported. Its modules are in `steps/core_db/`, and `workflow/20260320-build01_species_refAlt_coreDB.bash` still runs them. This workflow has not been run since the move, because it needs a reference genome and probe design file to test.
