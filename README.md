# HapApp_utils

[![License: Apache 2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://github.com/Breeding-Insight/HapApp_utils/blob/master/LICENSE)
[![Tests](https://github.com/Breeding-Insight/HapApp_utils/actions/workflows/tests.yml/badge.svg)](https://github.com/Breeding-Insight/HapApp_utils/actions/workflows/tests.yml)
[![Documentation](https://img.shields.io/badge/docs-readthedocs-blue.svg)](https://hapapp.readthedocs.io)

**Fixed microhaplotype IDs and versioned allele databases for targeted genotyping panels.**

`hapapp-utils` is a command-line tool for MADC reports from panels such as DArTag. A MADC gives many
alleles at a marker the same name, such as `Chr01_000120735|RefMatch`, so you cannot follow one
haplotype from project to project. `hapapp-utils` compares each allele with your panel's database,
reuses its ID if the database already holds that sequence, and issues a new ID if it does not. Every
project then adds to the same database.

## Highlights

- Builds a panel's first database from the reference genome and the DArT probe design template, and
  can make alleles longer than the reads.
- Gives each allele a permanent ID such as `Chr01_000120735|RefMatch_0002`, and writes each updated
  database as a new version.
- Takes its settings from one YAML file. Thresholds default to the values in the pipeline's step
  scripts, and you can change any of them for a single run from the command line.
- Checks tools, inputs and settings before it starts, then shows warnings and allele counts at each
  stage.
- Leaves your input files alone. Each run works on copies in its own folder, so you can compare
  settings side by side.
- Records the input checksums, tool versions and settings of every run.
- Keeps each pipeline step as its own Python module, which you can still run as a script.

## Installation

`hapapp-utils` needs Python 3.10 or newer and BLAST+. The `build-db` command also needs MMseqs2.
It is not on PyPI or Bioconda yet, so install it from source:

```bash
git clone https://github.com/Breeding-Insight/HapApp_utils.git
cd HapApp_utils
pip install .
```

Install the external tools with conda or Homebrew:

```bash
conda install -c conda-forge -c bioconda blast mmseqs2
```

```bash
brew install blast mmseqs2
```

`hapapp-utils check` tells you if it cannot find a tool.

## Quick Start

Build the panel's first database once:

```bash
hapapp-utils init build.yaml --workflow build-db   # write a settings file, then fill it in
hapapp-utils check -c build.yaml                   # check tools, inputs and settings
hapapp-utils build-db -c build.yaml                # build the database
```

Then process each MADC report against the latest database:

```bash
hapapp-utils init run.yaml                         # write a settings file, then fill it in
hapapp-utils check -c run.yaml
hapapp-utils assign -c run.yaml
```

A settings file for `assign` needs only your files and panel lengths:

```yaml
project_id: PROJECT01

inputs:
  madc: PROJECT01_MADC.csv
  marker_lut: species_snpID_lut.csv
  allele_db: species_allele_db_v003.fa
  match_counts: species_allele_db_v003_matchCnt_lut.txt

panel:
  design_length: 109
  sequencing_length: 109

output:
  dir: runs/PROJECT01
```

To try other thresholds, change them for one run and send the results to a new folder:

```bash
hapapp-utils assign -c run.yaml --set filter.min_samples=2 --set filter.min_reads=10 -o runs/trial_2x10
```

## Documentation

The full documentation is at [hapapp.readthedocs.io](https://hapapp.readthedocs.io):

- [Getting Started](https://hapapp.readthedocs.io/en/latest/getting-started/installation.html): installation, first runs and input files
- [Tutorials](https://hapapp.readthedocs.io/en/latest/tutorials/core-database.html): build a database, assign alleles and compare settings
- [Command Reference](https://hapapp.readthedocs.io/en/latest/commands/overview.html): every command and option
- [Settings](https://hapapp.readthedocs.io/en/latest/reference/settings.html): every setting and its default
- [How Allele IDs Work](https://hapapp.readthedocs.io/en/latest/reference/how-ids-work.html): how the filter and BLAST decide each allele's ID

## Commands

| Command | Description |
|---------|-------------|
| `hapapp-utils init` | Write a commented settings file |
| `hapapp-utils check` | Check tools, inputs and settings without running anything |
| `hapapp-utils build-db` | Build a panel's first allele database from the reference genome |
| `hapapp-utils assign` | Turn a MADC into fixed allele IDs and a new database version |
| `hapapp-utils hap-pa` | Build presence/absence tables from a report with fixed IDs |
| `hapapp-utils merge-pa` | Merge presence/absence matrices from several projects |

## Input Files

| File | Used by | Description |
|------|---------|-------------|
| MADC report | `build-db`, `assign` | Raw DArT report, or a report with only `AlleleID`, `CloneID`, `AlleleSequence` and samples |
| Probe design | `build-db` | The DArT template as `.xlsx`, `.csv` or tab-delimited `.txt` |
| Reference genome | `build-db` | FASTA. Its chromosome names must match the probe file, or you can supply a name map |
| Marker lookup table | `assign` | CSV of the MADC's marker names and your standard marker IDs. `build-db` writes one |
| Allele database and count table | `assign` | The latest versions from `build-db` or an earlier `assign` run |

## Output

`assign` writes these files to the run's folder:

| File | Contents |
|------|----------|
| `report/*_rename_v2.csv` | The report with fixed allele IDs, one column per sample |
| `versions/*_v00N.fa` | The new allele database, if the run found new alleles |
| `versions/*_v00N_matchCnt_lut.txt` | The new count table |
| `versions/*_v00N_process.readme` | The inputs, settings and everything each step printed |
| `hapSearch/hapPA/*_hapPA.csv` | Which alleles this project contains |
| `hapSearch/sample_hapStatus/*_hapStatus.csv` | Which samples contain each allele |
| `run_manifest.json` | Status, counts, warnings, tool versions and input checksums |

## Development

```bash
pip install -e ".[test]"
pytest
```

See [Contributing](https://hapapp.readthedocs.io/en/latest/development/contributing.html) for
the code layout and the documentation build.

## Citation

If you use HapApp_utils in your work, please cite:

> Zhao D, Lin M, Taniguti CH, Sandercock AM, et al. (2026). Standardized microhaplotype databases and frameworks for
> assessing and mining crop genetic diversity. *Theoretical and Applied Genetics* 139(9):244.
> https://doi.org/10.1007/s00122-026-05340-4

For the MADC format, see Breeding Insight's
[What is a MADC file?](https://softwarehub.breedinginsight.org/learn/what-is-a-madc-file.html).

**AI assistance:** AI assistance applies to version 2 of HapApp_utils and its documentation, not to version 1. Dongyan
Zhao wrote version 1, the step scripts and bash workflows. For version 2, Claude Opus 5.5 (Anthropic) helped convert those
workflows into the Python command-line tool, fix bugs in the step scripts, and write the tests and documentation. The
authors directed the work, reviewed and tested the results, and are ultimately responsible for the application.

## License

[Apache License 2.0](https://github.com/Breeding-Insight/HapApp_utils/blob/master/LICENSE)
