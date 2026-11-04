# hapapp-utils assign

Turn a MADC into fixed allele IDs and a new database version.

## Synopsis

```bash
hapapp-utils assign -c <FILE> [--set KEY=VALUE] [-o <DIR>] [--until <STAGE>] [-v]
```

## Description

`assign` runs in seven stages. [How Allele IDs Work](../reference/how-ids-work.md) covers the decisions
in detail.

| Stage | What it does |
|-------|--------------|
| `update-ids` | Replaces the MADC's marker names with your standard IDs, removes any duplicate tags, and stops if two allele IDs share a sequence |
| `filter` | Keeps alleles with enough reads in enough samples |
| `trim` | Removes adapters if reads run past the design length, and merges alleles that become identical |
| `blast` | Aligns the kept alleles to the database |
| `assign-ids` | Reuses IDs, issues new ones or removes each allele, then writes the report and the new database |
| `db-check` | Confirms each new ID has a sequence, and collapses any duplicate sequences |
| `hap-pa` | Builds the presence/absence tables |

The run copies its inputs into a new folder and works there, so your input database stays as it was.

## Arguments

| Argument | Description |
|----------|-------------|
| `-c`, `--config` | Settings file (required) |
| `--set` | Override a setting, as `KEY=VALUE`. You can repeat this option |
| `-o`, `--outdir` | Output folder, replacing `output.dir` |
| `--until` | Stop after this stage: `update-ids`, `filter`, `trim`, `blast`, `assign-ids`, `db-check` or `hap-pa` |
| `-v`, `--verbose` | Also show each step's output on screen, as well as in the process readme |

## Settings

[Settings](../reference/settings.md#assign) lists every setting. A minimal file:

```yaml
project_id: PROJECT01

inputs:
  madc: PROJECT01_MADC.csv
  marker_lut: species_snpID_lut.csv
  allele_db: species_allele_db_v003.fa
  match_counts: species_allele_db_v003_matchCnt_lut.txt

panel:
  design_length: 109
  sequencing_length: 109

output:
  dir: runs/PROJECT01
```

## Output

[Output Files](../reference/output-files.md#assign) lists every file.

## Examples

Run with the settings in the file:

```bash
hapapp-utils assign -c run.yaml
```

Try other thresholds in a new folder:

```bash
hapapp-utils assign -c run.yaml --set filter.min_samples=2 --set filter.min_reads=10 -o trials/2x10
```

Stop after filtering to see how many alleles pass:

```bash
hapapp-utils assign -c run.yaml --until filter -o trials/filter_only
```

The [Assign Allele IDs](../tutorials/assign-alleles.md) tutorial walks through a full run.
