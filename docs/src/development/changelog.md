# Changelog

## 2.0.0

The first version as a command-line tool. Each script from the earlier pipeline is now a module with the same functions
and arguments.

### New

- The `hapapp-utils` command, with `init`, `check`, `build-db`, `assign`, `hap-pa` and `merge-pa`.
- One YAML settings file per run, with command-line overrides.
- A new folder for each run, which leaves input files unchanged and records settings, tool versions and input checksums.
- Checks before each run, plus warnings and counts at each stage.
- `build-db` reads the DArT Excel template and accepts a chromosome name map.
- `build-db` no longer needs `esl-sfetch` or `seqkit`.

### Fixed

These fixes also apply to the bash workflows.

- `step01` counted two MADC summary columns as samples, so alleles passed the read filter with too few real samples.
- `step03` dropped `Other` alleles when trimming found duplicates, which left the count table ahead of the database.
  (Scott Brainard, [PR #1](https://github.com/Breeding-Insight/HapApp_utils/pull/1))
- `step03` kept untrimmed sequences in the database when trimming found no duplicates.
- `step03` and `step06` crashed on a path with a space.
- `db00` misnamed the lookup table and wrote tabs for a `.txt` probe file.
- `db02` wrote the count table in a different order each run.
- `db05` crashed if an allele had no BLAST hit.
- `util_hapDB_merge_presence_absence_PIhaps` inflated `Presence_count` when merging a matrix that had already been
  through a merge.
- `build01` had two path errors, and `build02` looked for the wrong duplicate report when it used `DUP_TAGS`.

### Changed

- The `Code_version` column and report file name now say `v2`.
- `FIRST_SAMPLE_COL` now reaches `step01` and `step05`.

### Added to the repository

- `step00_rm_dupTag_closeTag_rawMADC_v1.py`, which `build02` called but the repository lacked.
