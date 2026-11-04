# HapApp_utils

**hapapp-utils** gives each microhaplotype in a targeted genotyping panel a fixed ID, and keeps those
alleles in a versioned database across projects.

Panels such as DArTag report their results in a MADC table (Missing Allele Discovery Count). A MADC
gives many alleles at a marker the same name, such as `Chr01_000120735|RefMatch`, so you cannot follow
one haplotype from project to project. `hapapp-utils` compares each allele with your panel's database,
reuses its ID if the database already holds that sequence, and issues a new ID if it does not. Every
project then adds to the same database.

The output goes on to [BIGr](https://softwarehub.breedinginsight.org/software/bigr.html) and other tools for dosage calling and VCF
conversion. [Convert MADC data to a genotype-called VCF](https://softwarehub.breedinginsight.org/learn/madc-to-vcf.html) shows that workflow.

## Commands

| Command | Description |
|---------|-------------|
| [`init`](commands/init.md) | Write a commented settings file |
| [`check`](commands/check.md) | Check tools, inputs and settings without running anything |
| [`build-db`](commands/build-db.md) | Build a panel's first allele database from the reference genome |
| [`assign`](commands/assign.md) | Turn a MADC into fixed allele IDs and a new database version |
| [`hap-pa`](commands/hap-pa.md) | Build presence/absence tables from a report with fixed IDs |
| [`merge-pa`](commands/merge-pa.md) | Merge presence/absence matrices from several projects |

## What a run looks like

```text
$ hapapp-utils assign -c run.yaml
WARNING: 1 MADC marker(s) are not in the lookup table and will be dropped: ['Chr06_021528607']
[1/7] Update marker IDs and check for duplicate alleles
      2999 of 3000 markers kept, no duplicate alleles (0 s)
[2/7] Filter alleles on read counts
      9438 of 77491 alleles kept (at least 2 samples with at least 10 reads) (0 s)
[3/7] Remove adapters
      skipped: reads are not longer than the design (0 s)
[4/7] BLAST alleles against the database
      9438 of 9438 alleles hit their own marker (6 s)
[5/7] Assign fixed allele IDs
      708 already in the database, 8129 new, 601 removed by the BLAST thresholds or merged; 543 of the new ones are closest to the other allele type (0 s)
[6/7] Check the new database
      24529 sequences, was 16400 (1 s)
[7/7] Build presence/absence tables
      13176 alleles present in at least one sample (0 s)
```

## Where to start

- New to the tool: [Installation](getting-started/installation.md), then the
  [Quick Start](getting-started/quickstart.md).
- Starting a new panel: [Build a Core Database](tutorials/core-database.md).
- Processing a new order for a panel you already have:
  [Assign Allele IDs](tutorials/assign-alleles.md).
- Deciding on thresholds: [Compare Settings](tutorials/compare-settings.md).
