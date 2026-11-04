# hapapp-utils build-db

Build a panel's first allele database from the reference genome.

## Synopsis

```bash
hapapp-utils build-db -c <FILE> [--set KEY=VALUE] [-o <DIR>] [--until <STAGE>] [-v]
```

## Description

`build-db` takes the Ref and Alt alleles from the genome, because the MADC's sequences can contain IUPAC
codes. It runs in six stages:

| Stage | What it does |
|-------|--------------|
| `lookup` | Builds the marker lookup table from the probe file, keeping markers in the MADC, and updates the MADC's marker IDs |
| `flanks` | Fetches the sequence around each marker from the genome, and writes a Ref and an Alt version |
| `orient` | Aligns the MADC's amplicons to those flanks to find each marker's strand, then turns the minus-strand flanks around |
| `alleles` | Aligns the amplicons again to find where each one starts, then cuts `allele_length` bases from that point |
| `duplicates` | Uses MMseqs2 to find markers with near-identical Ref alleles, and removes the second marker of each pair |
| `finish` | Writes the database and count table as `_v001`, and checks allele lengths |

The `alleles` stage is why an allele can be longer than the read. An 81 bp read and an `allele_length` of
109 give a 109 bp allele, extended from the genome.

The tool also removes any marker left without both alleles. This happens when the probe's bases do not
match the genome, or when the amplicon does not align to its flank. A warning lists these markers.

## Arguments

| Argument | Description |
|----------|-------------|
| `-c`, `--config` | Settings file (required) |
| `--set` | Override a setting, as `KEY=VALUE`. You can repeat this option |
| `-o`, `--outdir` | Output folder, replacing `output.dir` |
| `--until` | Stop after this stage: `lookup`, `flanks`, `orient`, `alleles`, `duplicates` or `finish` |
| `-v`, `--verbose` | Also show each step's output on screen, as well as in the process readme |
| `--debug` | Show the full Python error if something unexpected fails |

## Settings

[Settings](../reference/settings.md#build-db) lists every setting. A minimal file:

```yaml
inputs:
  probe_design: probe_design.xlsx
  madc: project_MADC.csv
  genome: reference_genome.fa

database:
  name: species_allele_db
  allele_length: 109

output:
  dir: runs/core_db
```

## Output

[Output Files](../reference/output-files.md#build-db) lists every file. The lookup table, database and count
table become the inputs for [`assign`](assign.md).

## Example

The [Build a Core Database](../tutorials/core-database.md) tutorial walks through a full run.
