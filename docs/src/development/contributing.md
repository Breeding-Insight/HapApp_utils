# Contributing

Bug reports, documentation fixes and code are all welcome.

## Set up

```bash
git clone -b release https://github.com/Breeding-Insight/HapApp_utils.git
cd HapApp_utils

python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[test]"
```

With an editable install, the `hapapp-utils` command runs your working copy. You also need BLAST+ and MMseqs2 for
the workflow tests. See [Installation](../getting-started/installation.md).

## Branches

The repository has two long-lived branches:

| Branch | What it holds | Where the code lives |
|--------|---------------|----------------------|
| `master` | The step scripts and bash workflows. Dongyan Zhao maintains this branch | `scripts/`, `workflow/` |
| `release` | The same steps and workflows, plus the `hapapp-utils` package, its tests and this documentation | `src/hapapp_utils/`, `workflow/`, `tests/`, `docs/` |

Base a change to the command, the tests or this documentation on `release`. Base a change to a step script or bash
workflow on `master`, and it reaches `release` at the next merge. [Publish a release](#publish-a-release) covers the
tags on `release` that publish the package.

`master` has no tests of its own, so GitHub Actions tests its changes once they merge into `release`.

### Merge `master` into `release`

Merges go one way. Merging `release` into `master` would bring the whole package layout with it.

`release` only accepts pull requests, so make the merge on a new branch and open a pull request from it into `release`:

```bash
git fetch upstream
git switch -c merge-master upstream/release
git merge upstream/master
```

Merge that pull request with **Create a merge commit**. A squash merge would lose the link to `master`, so the next
merge would report the same conflicts again.

Git follows the moved files, so most edits to existing scripts land in `src/hapapp_utils/steps/` with no conflict.
These changes need more from you:

| Change on `master` | What Git does | What to do |
|--------------------|---------------|------------|
| A new option in a bash workflow | Merges it with no conflict, but the command does not get the option | [Add a setting](#add-a-setting) and pass it to the step |
| A new argument at the end of a script's list | Reports a conflict at `parse_args` | Keep the new argument and `parse_args(argv)`, then add a setting if the command needs it |
| A new script | Reports `CONFLICT (file location)` and puts the file in the new folder | Run `git add` on the file, then [add a step](#add-a-step) |
| An edit to the README | Reports a conflict, because `release` has its own README | Keep the `release` text, and copy over anything it lacks |

A clean merge does not mean the command has the change, because the command passes its own settings to the steps. After
each merge, read what changed in the bash workflows, then run the tests:

```bash
git diff ORIG_HEAD HEAD -- workflow/
pytest
```

### Send a fix to `master`

A fix on `release` can go to `master` as a cherry-pick. Git follows the moved files, so most fixes apply cleanly:

```bash
git switch master
git cherry-pick -x <commit>
```

The next merge into `release` can then report a conflict on those lines. Keep the copy on `release`.

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

`tests/test_cli_help.py` guards the command-line help. It runs or loads every example that `--help` shows, so a
wrong example fails the tests. It also holds the help to 80 columns.

## Change a step

On `release`, edit the module in `src/hapapp_utils/steps/`. On `master`, the same script lives under `scripts/`. Keep
its functions and arguments stable, because the workflows and the bash scripts call them. Mark a behavior change with
a `#### CHANGE:` comment that says why, as the existing code does.

Changes that touch alleles or read counts need a check on real data. A refactor should reproduce saved results byte
for byte:

```bash
HAPAPP_BASELINE_DIR=path/to/saved/run HAPAPP_BASELINE_CONFIG=path/to/run.yaml pytest -k real_data
```

When you fix a bug, add a test that fails without the fix.

## Add a step

Write the module with a `main(argv=None)` function, then call it from the matching stage in
`src/hapapp_utils/workflows/assign.py` or `build_db.py` using `call_step`.

To add a whole stage, add a `stage_<name>` method to the run class and list the stage in `STAGES`. Then add it with
a one-line explanation to `ASSIGN_STAGE_HELP` or `BUILD_DB_STAGE_HELP` in `cli.py`, and add it to the command's page.
A stage can return a short note to show on screen, and can add to `self.warnings`.

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

## Publish a release

Pushing a version tag from `release` publishes the package to PyPI through `.github/workflows/release.yml`.

1. Set `__version__` in `src/hapapp_utils/__init__.py` and add the version to the [Changelog](changelog.md).
2. Merge that change into `release` with a pull request.
3. Tag the merged commit and push the tag, using the new version number:

```bash
git tag v2.0.1
git push upstream v2.0.1
```

The workflow stops unless the tag matches `__version__` and the commit is on `release`. It then runs the tests on Linux,
builds the package, publishes it to PyPI and creates a GitHub release.

If the workflow stops before the publish step, fix the problem, delete the tag and tag again. PyPI never accepts the same
version number twice, so a problem found after publishing needs a new version.

The first release needs a pending publisher on [PyPI's publishing page](https://pypi.org/manage/account/publishing/),
with the project `hapapp-utils`, the owner `Breeding-Insight`, the repository `HapApp_utils`, the workflow `release.yml`
and the environment `pypi`.

After PyPI has the release, update `recipes/hapapp-utils/meta.yaml` in
[bioconda-recipes](https://github.com/bioconda/bioconda-recipes), unless Bioconda's bot has already opened a pull request
for it.

## Writing style

Write plainly, in the active voice, and name who does what: the tool, the step or you. Give real numbers and real
output instead of general claims.
