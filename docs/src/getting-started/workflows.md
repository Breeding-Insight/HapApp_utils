# Choosing a Workflow

Your starting point decides which command you need.

| You have | You want | Use |
|----------|----------|-----|
| A probe design file, a reference genome and one MADC | The panel's first database | [`build-db`](../commands/build-db.md) |
| A MADC and the panel's database | Fixed allele IDs and an updated database | [`assign`](../commands/assign.md) |
| A MADC and a database, and you are unsure of the thresholds | Results to compare | [`assign`](../commands/assign.md) with `--set` and `-o` |
| A report with fixed IDs | Presence/absence tables only | [`hap-pa`](../commands/hap-pa.md) |
| Presence/absence matrices from several projects | One matrix | [`merge-pa`](../commands/merge-pa.md) |

## How the commands fit together

```text
probe design + genome + MADC
            |
        build-db
            |
   allele database v001 ---------+
   count table v001              |
   marker lookup table           |
                                 v
            first project's MADC -> assign -> report with fixed IDs
                                 |         -> presence/absence tables
                                 v
                      allele database v002 --+
                                             v
                 second project's MADC -> assign -> ...
```

You run `build-db` once per panel. You run `assign` once per project, and each run reads the database
version the last run wrote.

## Order matters

`assign` numbers new alleles in the order it finds them, and a sequence keeps its ID once it has one.
Run projects in the order you want their IDs issued, and always use the latest database. A new run on
an older version issues different IDs for the same sequences.

Each run leaves its input database unchanged and writes the new version to its own folder. To accept a
run, copy its new database and count table to where you keep them, and use them for the next project.
A trial run has no effect until you do this.

## Running a single step

Every step of both workflows is its own Python module. To rerun or debug one step, call it directly:

```bash
python -m hapapp_utils.steps.assign.step01_filter_missing_AND_ext_matchAlleles_from_madc_v1 --help
```

[Architecture](../development/architecture.md) lists which modules each stage runs.
