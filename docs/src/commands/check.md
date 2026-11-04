# hapapp-utils check

Check tools, inputs and settings without running anything.

## Synopsis

```bash
hapapp-utils check -c <FILE> [--set KEY=VALUE] [-o <DIR>]
```

## Description

`check` looks at the settings file's sections to tell whether the file is for `build-db` or `assign`, then
runs the checks for that command. Both commands run the same checks when they start, and refuse to run if
any finds a problem.

Each line reports `ok`, `WARNING` or `PROBLEM`. A warning lets the run go ahead. A problem stops it, and
`check` exits with code 1.

| Check | Applies to |
|-------|------------|
| The external programs exist on your PATH | both |
| Each input file exists and is readable text | both |
| The MADC has a header row and sample columns | both |
| The output folder is new or empty, and you can write to it | both |
| `panel.first_sample_column` fits the header | `assign` |
| `panel.sequencing_length` matches the MADC's Ref alleles | `assign` |
| The lookup table has the MADC's markers | `assign` |
| The database has `Ref_0001` and `Alt_0002` alleles for those markers | `assign` |
| The count table agrees with the database | `assign` |
| The probe file has the expected columns, and shares markers with the MADC | `build-db` |
| The genome has the chromosomes those markers use | `build-db` |
| `database.flank_length` and `database.allele_length` suit the reads | `build-db` |

[Warnings and Errors](../reference/messages.md) explains each message and what to do about it.

## Arguments

| Argument | Description |
|----------|-------------|
| `-c`, `--config` | Settings file (required) |
| `--set` | Override a setting, as `KEY=VALUE`. You can repeat this option |
| `-o`, `--outdir` | Output folder, replacing `output.dir` |
| `--debug` | Show the full Python error if something unexpected fails |

## Example

```bash
hapapp-utils check -c run.yaml
```
