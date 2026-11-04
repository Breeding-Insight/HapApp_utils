# Command Overview

```text
hapapp-utils [--version] [--debug] <command> [options]
```

| Command | Description |
|---------|-------------|
| [`init`](init.md) | Write a commented settings file |
| [`check`](check.md) | Check tools, inputs and settings without running anything |
| [`build-db`](build-db.md) | Build a panel's first allele database from the reference genome |
| [`assign`](assign.md) | Turn a MADC into fixed allele IDs and a new database version |
| [`hap-pa`](hap-pa.md) | Build presence/absence tables from a report with fixed IDs |
| [`merge-pa`](merge-pa.md) | Merge presence/absence matrices from several projects |

Run `hapapp-utils <command> --help` to see a command's options.

## Options for every command

| Option | Description |
|--------|-------------|
| `--version` | Print the version and stop |
| `--debug` | Show the full Python error if something unexpected fails |

These options go before the command, as in `hapapp-utils --debug assign -c run.yaml`.

## Exit codes

| Code | Meaning |
|------|---------|
| 0 | The command finished |
| 1 | An error stopped the command |
| 2 | The command line itself had a mistake |
| 3 | `assign` found the same sequence under different allele IDs in the MADC |
| 130 | You interrupted the run |
