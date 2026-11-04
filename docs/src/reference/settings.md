# Settings

Settings live in a YAML file. [`init`](../commands/init.md) writes a commented one for you.

## Precedence

A value on the command line replaces the value in the file, and a value in the file replaces the default:

```bash
hapapp-utils assign -c run.yaml --set filter.min_reads=30 -o trials/min30
```

Relative paths in the file start from the file's own folder. Relative paths on the command line start from
your current folder.

## assign

The Bash variable column names the matching variable in the `build02` bash workflow.

### Project

| Setting | Required | Default | Bash variable | Description |
|---------|----------|---------|---------------|-------------|
| `project_id` | yes | | `PROJECT_ID` | Names the project in output files |

### Inputs for assign

| Setting | Required | Default | Bash variable | Description |
|---------|----------|---------|---------------|-------------|
| `inputs.madc` | yes | | `REPORT` | The MADC report |
| `inputs.marker_lut` | yes | | `MARKERID_LUT` | The marker lookup table |
| `inputs.allele_db` | yes | | `ALLELE_DB_DIR`, `ALLELE_DB` | The latest allele database |
| `inputs.match_counts` | yes | | `MATCHCNT_LUT` | The latest count table |
| `inputs.dup_tags` | no | none | `DUP_TAGS` | Markers to remove, one per line |
| `inputs.previous_hap_pa` | no | none | `LATEST_HAPPA_MATRIX` | Matrix to merge this project into |

### Panel

| Setting | Required | Default | Bash variable | Description |
|---------|----------|---------|---------------|-------------|
| `panel.design_length` | yes | | `DESIGN_LEN` | Length of the panel's insert. Reads longer than this go through adapter trimming |
| `panel.sequencing_length` | yes | | `SEQ_LEN` | Length of the reads |
| `panel.first_sample_column` | no | from the header | `FIRST_SAMPLE_COL` | Column number of the first sample. A value inside the summary columns stops the run |
| `panel.insert_size_column` | no | none | `LUT_INSERT_SIZE_COL` | Column in the lookup table holding each marker's insert size, for mixed panels |
| `panel.adapter` | no | `ACCGATCTG` | in the cutadapt command | Adapter sequence to trim |

### Filter

| Setting | Required | Default | Bash variable | Description |
|---------|----------|---------|---------------|-------------|
| `filter.min_samples` | no | 10 | `MIN_SAMPLES` | Samples that must reach `min_reads`. The tool uses 5% of the samples, rounded up, if that is smaller |
| `filter.min_reads` | no | 2 | `MIN_READS` | Reads a sample needs to count |

### BLAST

| Setting | Required | Default | Bash variable | Description |
|---------|----------|---------|---------------|-------------|
| `blast.min_coverage` | no | 90 | `COV` | Percent of the allele that must align to its closest match |
| `blast.min_identity` | no | 90 | `IDEN` | Percent identity of that alignment |
| `blast.allele_type_truth` | no | `report` | `ALLELE_TYPE` | Which side wins when the MADC and database disagree on allele type: `report` or `db` |
| `blast.threads` | no | 6 | in the blastn command | BLAST threads |

### Presence/absence

| Setting | Required | Default | Bash variable | Description |
|---------|----------|---------|---------------|-------------|
| `hap_pa.min_reads` | no | 10 | `MIN_READS_FOR_PRESENCE` | Reads a sample needs for an allele to count as present |
| `hap_pa.descriptor` | no | `project_id` | `PROJECT_DESCRIPTOR` | This project's column name in the matrix |

### Output for assign

| Setting | Required | Default | Bash variable | Description |
|---------|----------|---------|---------------|-------------|
| `output.dir` | yes | | none | Output folder, which must be new or empty |
| `output.code_version` | no | `v2` | `CODE_VER` | The value in the report's `Code_version` column and file name |

The filter and BLAST defaults come from the `step01` and `step05` scripts, and the filter's also match the
project's written method. The values at the top of the `build02` bash workflow differ (1 sample, 10 reads, 80%
coverage, 80% identity), because each project edits that block.

## build-db

The Bash variable column names the matching variable in the `build01` bash workflow.

### Inputs for build-db

| Setting | Required | Default | Bash variable | Description |
|---------|----------|---------|---------------|-------------|
| `inputs.probe_design` | yes | | `PROBE_FILE` | The DArT probe design template, as `.xlsx`, `.csv` or tab-delimited `.txt` |
| `inputs.madc` | yes | | `REPORT` | Any MADC report from the panel |
| `inputs.genome` | yes | | `REF_GENOME` | The reference genome FASTA, not gzipped |
| `inputs.chromosome_map` | no | none | none | Translates the probe file's chromosome names to the genome's |

### Database

| Setting | Required | Default | Bash variable | Description |
|---------|----------|---------|---------------|-------------|
| `database.name` | yes | | `ALLELE_DB`, `MATCH_CNT` | Start of the file names. The first version adds `_v001` |
| `database.allele_length` | yes | | `REF_LEN` | Length of the Ref and Alt alleles. This can exceed the read length |
| `database.flank_length` | no | 150 | `FLANK_LEN` | Bases to fetch on each side of a marker |
| `database.dup_min_identity` | no | 0.95 | in the mmseqs command | Markers whose Ref alleles reach this identity count as duplicates |
| `database.dup_min_coverage` | no | 0.99 | in the mmseqs command | Coverage those alleles must also reach |
| `database.threads` | no | 6 | in the blastn and mmseqs commands | Threads for BLAST and MMseqs2 |

### Output for build-db

| Setting | Required | Default | Bash variable | Description |
|---------|----------|---------|---------------|-------------|
| `output.dir` | yes | | `ALLELE_DB_DIR` | Output folder, which must be new or empty |

## Values you cannot set

These stay fixed in both the bash workflows and the CLI:

| Value | Setting |
|-------|---------|
| BLAST E-value | `1e-5` |
| BLAST maximum targets | 15 for `assign`, 10 for `build-db` |
| cutadapt error rate, repeat count and minimum overlap | `-e 0.2`, `-n 1`, `--overlap 5` |

## YAML tips

- Words such as `NO`, `yes` and `on` mean true or false in YAML, so put a value like that in quotes:
  `project_id: "NO"`. The tool stops with a message if it reads a text setting as true or false.
- A number with leading zeros loses them, so write `project_id: "007"` in quotes.
- Paths with spaces need no quotes.
