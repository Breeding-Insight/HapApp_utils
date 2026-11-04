# Quick Start

You build a panel's database once, then process each new MADC against the latest version.

## 1. Build the first database

Skip this step if your panel already has an allele database.

Write a settings file:

```bash
hapapp-utils init build.yaml --workflow build-db
```

Open `build.yaml` and fill in your files:

```yaml
inputs:
  probe_design: probe_design.xlsx
  madc: DCl26-12123_MADC.csv
  genome: reference_genome.fa

database:
  name: redClover_allele_db
  allele_length: 109

output:
  dir: runs/core_db
```

Check it, then run it:

```bash
hapapp-utils check -c build.yaml
hapapp-utils build-db -c build.yaml
```

The command prints the paths of the three files you need next: the allele database, the count table
and the lookup table.

## 2. Assign allele IDs

Write a second settings file:

```bash
hapapp-utils init run.yaml
```

Point it at the MADC and the files from step 1. Copy those three files to a place you keep them, such
as a `databases` folder:

```yaml
project_id: DCl26-12123

inputs:
  madc: DCl26-12123_MADC.csv
  marker_lut: databases/redClover_snpID_lut.csv
  allele_db: databases/redClover_allele_db_v001.fa
  match_counts: databases/redClover_allele_db_v001_matchCnt_lut.txt

panel:
  design_length: 109
  sequencing_length: 109

output:
  dir: runs/DCl26-12123
```

Check it, then run it:

```bash
hapapp-utils check -c run.yaml
hapapp-utils assign -c run.yaml
```

## 3. Use the results

The report with fixed IDs is `report/*_rename_v2.csv` in the output folder. The new database is in
`versions/`. Point the next project's settings at that new version.

## Next steps

- [Choosing a Workflow](workflows.md) shows how the commands fit together.
- [Settings](../reference/settings.md) lists every setting and its default.
- [Warnings and Errors](../reference/messages.md) explains what the tool may tell you.
