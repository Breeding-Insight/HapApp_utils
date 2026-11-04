# Build a Core Database

This tutorial builds the first allele database for a panel, using the red clover DArTag panel as the
example. It has 3,000 markers on the ARS_RC_1.1 assembly.

A MADC cannot supply these alleles directly, because its Ref and Alt sequences can contain IUPAC codes.
The red clover MADC has 334 markers with these codes. `build-db` uses the MADC only to find where each
amplicon sits, then takes the allele sequences from the genome.

## What you need

- The probe design file you submitted to DArT.
- The reference genome the probes came from.
- Any one MADC from the panel.

## Write the settings

```bash
hapapp-utils init build.yaml --workflow build-db
```

Fill in your files. This genome uses NCBI chromosome names and the probe file uses `Chr01`, so the
settings also point to a [chromosome map](../getting-started/input-files.md#chromosome-map):

```yaml
inputs:
  probe_design: Red_Clover_DArT_ProbeDesign.xlsx
  madc: DCl26-12123_MADC.csv
  genome: GCF_020283565.1_ARS_RC_1.1_genomic.fna
  chromosome_map: red_clover_chromosome_map.csv

database:
  name: redClover_allele_db
  allele_length: 109
  flank_length: 150

output:
  dir: runs/core_db
```

Set `allele_length` to at least your read length. A longer value makes the database ready for longer
reads later, because the tool extends each allele from the genome. It needs enough flanking sequence to
do this, so keep `flank_length` at or above `allele_length`.

## Check the settings

```bash
hapapp-utils check -c build.yaml
```

```text
  ok       makeblastdb: makeblastdb: 2.17.0+
  ok       blastn: blastn: 2.17.0+
  ok       mmseqs: 18-8cc5c
  ok       Probe design: 12377 markers, 3000 of them in the MADC (3000 markers, 94 samples)
  ok       Genome: 90 sequences, including all 7 chromosomes these markers use
  ok       Output folder: runs/core_db

0 problem(s), 0 warning(s).
```

## Run it

```bash
hapapp-utils build-db -c build.yaml
```

```text
[1/6] Build the marker lookup table
      3000 markers (1 s)
[2/6] Fetch flanking sequences from the genome
      3000 of 3000 markers match the genome (1 s)
[3/6] Find the strand of each amplicon
      1353 markers on the minus strand (5 s)
[4/6] Fetch Ref and Alt alleles
      6000 alleles (5 s)
[5/6] Remove markers with near-identical Ref alleles
      1 duplicate(s), 0 without both alleles (1 s)
[6/6] Write the new database
      5998 alleles from 2999 markers (1 s)
```

## Read the results

Stage 5 found one duplicate. `report/*_markers_removed.tsv` shows which marker the tool removed and why:

```text
Chr06_017793919|Ref_0001	Chr06_021528607|Ref_0001
```

The second marker's Ref allele matches the first's at 95% identity or more, so reads could not tell
them apart. The tool removed `Chr06_021528607` from the database, the count table and the lookup table.

You now have three files for `assign`:

- `versions/redClover_allele_db_v001.fa`
- `versions/redClover_allele_db_v001_matchCnt_lut.txt`
- `probe/*_rmDupTag_snpID_lut.csv`, or `probe/*_snpID_lut.csv` if the run removed no markers

Copy them to where you keep your databases, then continue with
[Assign Allele IDs](assign-alleles.md).
