# hapapp-utils hap-pa

Build presence/absence tables from a report with fixed IDs.

## Synopsis

```bash
hapapp-utils hap-pa --madc <FILE> --descriptor <NAME> [--min-reads <N>]
```

## Description

[`assign`](assign.md) runs this for you as its last stage. Use this command to rebuild the tables with
another read threshold, without repeating the run. `--min-reads` here does the same job as
[`hap_pa.min_reads`](../reference/settings.md#assign) in the settings file.

It writes two files next to the input report:

| File | Contents |
|------|----------|
| `<report>_hapStatus.csv` | One row per allele and one column per sample, with 1 if the sample has at least `--min-reads` reads, and 0 if not |
| `<report>_hapPA.csv` | The alleles present in at least one sample, with one column named after `--descriptor` |

## Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--madc` | required | Report with fixed allele IDs, without the `Code_version` column |
| `--descriptor` | required | Column name for this project in the matrix |
| `--min-reads` | 10 | Reads needed to count as present |
| `--debug` | off | Show the full Python error if something unexpected fails |

## Example

```bash
hapapp-utils hap-pa --madc report/PROJECT01_MADC_snpID_rename.csv --descriptor PROJECT01 --min-reads 20
```
