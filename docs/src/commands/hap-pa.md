# hapapp-utils hap-pa

Build presence/absence tables from a report with fixed IDs.

## Synopsis

```bash
hapapp-utils hap-pa --madc <FILE> --descriptor <NAME> [--read-threshold <N>]
```

## Description

[`assign`](assign.md) runs this for you as its last stage. Use this command to rebuild the tables with
another threshold, without repeating the run.

It writes two files next to the input report:

| File | Contents |
|------|----------|
| `<report>_hapStatus.csv` | One row per allele and one column per sample, with 1 if the sample has at least the threshold of reads, and 0 if not |
| `<report>_hapPA.csv` | The alleles present in at least one sample, with one column named after `--descriptor` |

## Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--madc` | required | Report with fixed allele IDs, without the `Code_version` column |
| `--descriptor` | required | Column name for this project in the matrix |
| `--read-threshold` | 10 | Reads needed to count as present |

## Example

```bash
hapapp-utils hap-pa --madc report/DCl26-12123_MADC_snpID_rename.csv --descriptor DCl26-12123 --read-threshold 20
```
