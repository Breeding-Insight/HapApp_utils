# Assign Allele IDs

This tutorial processes one MADC against a panel's database, using a red clover order with 94 samples.

## Write the settings

```bash
hapapp-utils init run.yaml
```

```yaml
project_id: DCl26-12123

inputs:
  madc: DCl26-12123_MADC.csv
  marker_lut: redClover_snpID_lut.csv
  allele_db: redClover_allele_db_v003.fa
  match_counts: redClover_allele_db_v003_matchCnt_lut.txt

panel:
  design_length: 109
  sequencing_length: 109

filter:
  min_samples: 2
  min_reads: 10

blast:
  min_coverage: 80
  min_identity: 85

output:
  dir: runs/DCl26-12123
```

These thresholds suit this multi-species trial. Leave out the `filter` and `blast` sections to use the
[defaults](../reference/settings.md#filter). [Compare Settings](compare-settings.md) shows how to choose.

## Check the settings

```bash
hapapp-utils check -c run.yaml
```

```text
  ok       MADC: 94 samples, first sample in column 17
  WARNING  1 MADC marker(s) are not in the lookup table and will be dropped: ['Chr06_021528607']
  ok       Count table and allele database agree
  ok       Output folder: runs/DCl26-12123

0 problem(s), 1 warning(s).
```

The warning is expected here, because the panel excluded this marker on purpose. A warning does not stop
a run. A problem does.

## Run it

```bash
hapapp-utils assign -c run.yaml
```

```text
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

## Read the counts

| Stage | Line | What it tells you |
|-------|------|-------------------|
| 1 | 2999 of 3000 markers kept | The lookup table left out one marker |
| 2 | 9438 of 77491 alleles kept | Most alleles in a MADC have few reads, and the filter removed them |
| 3 | skipped | The reads are no longer than the design, so they contain no adapter |
| 4 | 9438 of 9438 | Every allele aligned to its own marker in the database |
| 5 | 708 already in the database | These sequences kept their existing IDs |
| 5 | 8129 new | These sequences received new IDs |
| 5 | 601 removed | These failed the coverage or identity thresholds |
| 6 | 24529, was 16400 | The new database version holds the old alleles plus the new ones |

Stage 5 also counts new alleles whose closest match in the database has the other Ref or Alt type. Here that
is 543, and `versions/*_mismatchedType_alleles.csv` lists them.

## Look at the results

The report with fixed IDs has one column per sample:

```text
Code_version,AlleleID,CloneID,AlleleSequence,5236,5237,...
v2,Chr01_000120735|Ref_0001,Chr01_000120735,GGTGATGAGTGG...,0,12,...
v2,Chr01_000120735|RefMatch_0001,Chr01_000120735,GGTGATGAGTGA...,0,0,...
```

The process readme, `versions/redClover_allele_db_v004_process.readme`, opens with the input files and
settings for this run, followed by everything each step printed.

## Accept the new database

The run left your input database unchanged. If you accept these results, copy these files to where you
keep your databases, and use them for the next project:

- `versions/redClover_allele_db_v004.fa`
- `versions/redClover_allele_db_v004_matchCnt_lut.txt`

To add this project to a matrix of earlier projects, set
[`inputs.previous_hap_pa`](../reference/settings.md#inputs-for-assign) before the run, or use
[`merge-pa`](../commands/merge-pa.md) afterward.
