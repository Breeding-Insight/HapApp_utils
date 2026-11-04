# hapapp-utils merge-pa

Merge presence/absence matrices from several projects.

## Synopsis

```bash
hapapp-utils merge-pa <MATRIX> [<MATRIX> ...] -o <FILE>
```

## Description

`merge-pa` joins any number of `_hapPA.csv` files, or earlier merged matrices, into one table:

- It keeps every allele from every input, and fills gaps with 0.
- It adds a `CloneID` column and a `Presence_count` column, which counts the projects with each allele.
- It adds any missing `Ref_0001` and `Alt_0002` rows, so every marker has both.
- It sorts each marker's alleles as Ref, Alt, RefMatch, AltMatch, then any others.

You can merge an earlier merged matrix again. The tool recounts `Presence_count` each time.

[`assign`](assign.md) can do the merge for you if you set
[`inputs.previous_hap_pa`](../reference/settings.md#inputs-for-assign).

## Arguments

| Argument | Description |
|----------|-------------|
| `matrices` | The `_hapPA.csv` files or earlier merged matrices to merge |
| `-o`, `--output` | The merged matrix to write (required) |
| `--debug` | Show the full Python error if something unexpected fails |

## Example

```bash
hapapp-utils merge-pa earlier_hapPA.csv this_project_hapPA.csv -o all_projects_hapPA.csv
```
