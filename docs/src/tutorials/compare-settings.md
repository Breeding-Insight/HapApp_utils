# Compare Settings

The read filter and the BLAST thresholds decide which alleles get a permanent ID, so test them on your
own data. Each run gets its own folder and leaves the input database alone, so a trial run changes
nothing.

## Run the trials

Start from one settings file. Use `--set` to change a threshold for one run, and `-o` to give the run
its own folder:

```bash
hapapp-utils assign -c run.yaml --set filter.min_samples=2 --set filter.min_reads=10 -o trials/2x10
hapapp-utils assign -c run.yaml --set filter.min_samples=1 --set filter.min_reads=30 -o trials/1x30
hapapp-utils assign -c run.yaml --set blast.min_coverage=90 -o trials/cov90
```

Each run saves its counts in `run_manifest.json`, and its settings in `resolved_config.yaml`.

## An example

These runs used a red clover panel of 2,999 markers on 94 samples from four *Trifolium* species, starting
from a database of 16,400 sequences. Of the 77,491 alleles in the MADC, 41,119 had only one read across all
samples.

| Read filter | Coverage / identity | Kept by filter | New alleles | Removed at BLAST | Database size |
|-------------|---------------------|----------------|-------------|------------------|---------------|
| 5 samples, 2 reads (default) | 90 / 90 | 10,674 | 7,784 | 2,188 | 24,184 |
| 2 samples, 10 reads | 90 / 85 | 9,438 | 7,410 | 1,320 | 23,810 |
| 2 samples, 10 reads | 80 / 85 | 9,438 | 8,129 | 601 | 24,529 |
| 1 sample, 30 reads | 80 / 85 | 8,889 | 7,766 | 414 | 24,166 |

The default of 10 samples becomes 5 here, because the filter uses 5% of the samples when that is smaller.

## What each threshold trades off

The number of samples decides whether a single sample can add an allele. In this trial, one outlier sample
supplied 70% of the alleles that only the 1 sample rule kept, and that sample did not resemble any of the
four species. Requiring two samples keeps such a sample out of the database.

The number of reads sets how much evidence each sample must supply. A low value such as 2 keeps more alleles
that several low-depth samples share.

BLAST coverage decides how much of an allele must align. BLAST stops an alignment at an indel or at a
cluster of mismatches near the 3' end, so a real allele can fail at 90%. Lowering coverage to 80 here
recovered 719 alleles carrying 795,000 reads, and lost no allele that passed at 90.

BLAST identity decides how far an allele may differ from its closest match in the database. A panel you use
on related species may need a lower value.

## Choosing

Look at the counts on screen and in `run_manifest.json`, and ask these questions:

- Does most of the database growth come from one sample? Check the samples' read totals and shared alleles.
- Did the BLAST step remove alleles with many reads? Those are more likely real than alleles the read
  filter removed.
- How much do the usable calls per sample change? The `hapStatus` table shows this for each setting.

When you settle on values, write them in your settings file, so `resolved_config.yaml` and every later
run record them.
