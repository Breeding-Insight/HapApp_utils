# Installation

`hapapp-utils` needs Python 3.10 or newer. It also uses these programs, which you install separately:

| Program | Needed for | conda package | Homebrew formula |
|---------|------------|---------------|------------------|
| BLAST+ (`blastn`, `makeblastdb`) | `build-db`, `assign` | `blast` | `blast` |
| MMseqs2 (`mmseqs`) | `build-db` | `mmseqs2` | `mmseqs2` |

Installing the package also installs cutadapt, which `assign` uses when reads run past the panel's
design length.

## From source

The package is not on PyPI or Bioconda yet.

```bash
git clone -b release https://github.com/Breeding-Insight/HapApp_utils.git
cd HapApp_utils
pip install .
```

Then install the external programs with one of these commands:

```bash
conda install -c conda-forge -c bioconda blast mmseqs2
```

```bash
brew install blast mmseqs2
```

## Check the installation

```bash
hapapp-utils --version
hapapp-utils --help
```

The [`check`](../commands/check.md) command reports any program it cannot find, so run it before your
first analysis.

## Working on the code

To edit the code and see your changes straight away, install it in editable mode with the test tools:

```bash
pip install -e ".[test]"
```

[Contributing](../development/contributing.md) covers the rest.
