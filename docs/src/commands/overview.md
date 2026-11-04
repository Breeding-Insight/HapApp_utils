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

Run `hapapp-utils --help` for this list and a typical workflow. Run `hapapp-utils <command> --help` for a command's
options, stages and examples.

## Options for every command

| Option | Description |
|--------|-------------|
| `--version` | Print the version and stop |
| `--debug` | Show the full Python error if something unexpected fails |

`--version` goes before the command. `--debug` can go before or after it, so `hapapp-utils --debug assign -c run.yaml` and
`hapapp-utils assign -c run.yaml --debug` do the same thing.

## Exit codes

| Code | Meaning |
|------|---------|
| 0 | The command finished |
| 1 | An error stopped the command |
| 2 | The command line itself had a mistake, or you gave no command |
| 3 | `assign` found the same sequence under different allele IDs in the MADC |
| 130 | You interrupted the run |
