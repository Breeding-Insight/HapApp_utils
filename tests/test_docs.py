"""Keep the documentation in step with the code. These tests skip if the docs are not next to the tests."""

import argparse
import re
from pathlib import Path

import pytest

from hapapp_utils.cli import ASSIGN_STAGES, BUILD_DB_STAGES, _build_parser
from hapapp_utils.config import BUILD_DB_SCHEMA, SCHEMA

DOCS = Path(__file__).parent.parent / "docs" / "src"
pytestmark = pytest.mark.skipif(not DOCS.exists(), reason="documentation is not available")


def subcommands():
    parser = _build_parser()
    action = next(item for item in parser._actions if isinstance(item, argparse._SubParsersAction))
    return action.choices


@pytest.mark.parametrize("command", ["init", "check", "build-db", "assign", "hap-pa", "merge-pa"])
def test_each_command_page_lists_its_options(command):
    page = (DOCS / "commands" / f"{command}.md").read_text()
    for action in subcommands()[command]._actions:
        if isinstance(action, argparse._HelpAction):
            continue
        names = action.option_strings or [action.dest]
        for name in names:
            assert f"`{name}`" in page, f"{command} page does not mention {name}"


def test_commands_all_have_a_page_and_appear_in_the_overview():
    overview = (DOCS / "commands" / "overview.md").read_text()
    for command in subcommands():
        assert (DOCS / "commands" / f"{command}.md").exists()
        assert f"`{command}`" in overview


def test_the_settings_page_lists_every_setting():
    page = (DOCS / "reference" / "settings.md").read_text()
    for key in list(SCHEMA) + list(BUILD_DB_SCHEMA):
        assert f"`{key}`" in page, f"settings page does not mention {key}"


def test_the_settings_page_shows_the_real_defaults():
    page = (DOCS / "reference" / "settings.md").read_text()
    for key, (_, default) in {**SCHEMA, **BUILD_DB_SCHEMA}.items():
        if isinstance(default, (int, float)) and not isinstance(default, bool):
            shown = int(default) if float(default).is_integer() else default
            row = next(line for line in page.splitlines() if line.startswith(f"| `{key}`"))
            assert f"| {shown} |" in row, f"{key} should show a default of {shown}: {row}"


@pytest.mark.parametrize("command,stages", [("assign", ASSIGN_STAGES), ("build-db", BUILD_DB_STAGES)])
def test_command_pages_describe_each_stage(command, stages):
    page = (DOCS / "commands" / f"{command}.md").read_text()
    for stage in stages:
        assert f"| `{stage}` |" in page, f"{command} page does not describe {stage}"


def pages():
    return sorted(DOCS.rglob("*.md"))


def anchors(path):
    found = set()
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            title = line.lstrip("#").strip().lower()
            found.add(re.sub(r"[^\w\- ]", "", title).replace(" ", "-"))
    return found


@pytest.mark.parametrize("page", pages(), ids=lambda path: str(path.relative_to(DOCS)))
def test_internal_links_work(page):
    for target in re.findall(r"\]\((?!https?://|<)([^)\s]+)\)", page.read_text()):
        path, _, anchor = target.partition("#")
        linked = (page.parent / path).resolve() if path else page
        assert linked.exists(), f"{page.name} links to a missing file: {target}"
        if anchor:
            assert anchor in anchors(linked), f"{page.name} links to a missing heading: {target}"


def test_every_page_is_in_the_summary():
    summary = (DOCS / "SUMMARY.md").read_text()
    for page in pages():
        relative = str(page.relative_to(DOCS))
        if relative != "SUMMARY.md":
            assert f"({relative})" in summary, f"{relative} is not in SUMMARY.md"
