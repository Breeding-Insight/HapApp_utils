# Contributing

Bug reports, documentation fixes and code are all welcome.

## Set up

```bash
git clone https://github.com/Breeding-Insight/HapApp_utils.git
cd HapApp_utils

python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[test]"
```

With an editable install, the `hapapp-utils` command runs your working copy. You also need BLAST+ and MMseqs2 for
the workflow tests. See [Installation](../getting-started/installation.md).

## Run the tests

```bash
pytest
```

The workflow tests skip themselves if BLAST+ or MMseqs2 is missing. The full suite takes about 8 minutes. To run part
of it, use `-k`:

```bash
pytest -k trimming
```

GitHub Actions runs the tests on every push and pull request, on Python 3.10 to 3.13 and on macOS
(`.github/workflows/tests.yml`). It also builds the package and this documentation.

`tests/test_docs.py` guards the documentation. It checks that each command page lists that command's options, that
[Settings](../reference/settings.md) lists every setting, and that internal links work. Update the pages when you
add or change an option or setting.

## Change a step

Edit the module in `src/hapapp_utils/steps/`. Keep its functions and arguments stable, because the workflows and
the bash scripts call them. Mark a behavior change with a `#### CHANGE:` comment that says why, as the existing code
does.

Changes that touch alleles or read counts need a check on real data. A refactor should reproduce saved results byte
for byte:

```bash
HAPAPP_BASELINE_DIR=path/to/saved/run HAPAPP_BASELINE_CONFIG=path/to/run.yaml pytest -k real_data
```

When you fix a bug, add a test that fails without the fix.

## Add a step

Write the module with a `main(argv=None)` function, then call it from the matching stage in
`src/hapapp_utils/workflows/assign.py` or `build_db.py` using `call_step`.

To add a whole stage, add a `stage_<name>` method to the run class and list the stage in `STAGES`. Then list it in
`ASSIGN_STAGES` or `BUILD_DB_STAGES` in `cli.py`, and add it to the command's page. A stage can return a short note
to show on screen, and can add to `self.warnings`.

## Add a setting

Add it to `SCHEMA` or `BUILD_DB_SCHEMA` in `config.py`, and to the matching template in that file. Then add it to
[Settings](../reference/settings.md).

## Build the documentation

This site uses [mdBook](https://rust-lang.github.io/mdBook/). Install it with `cargo install mdbook`, then build and
preview:

```bash
mdbook serve docs --open
```

Read the Docs builds the published site from `.readthedocs.yaml`. The colors in `docs/theme/breeding-insight.css`
come from the Breeding Insight palette.

## Writing style

Write plainly, in the active voice, and name who does what: the tool, the step or you. Give real numbers and real
output instead of general claims.
