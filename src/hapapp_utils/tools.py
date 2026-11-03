"""Finding and running the external programs the pipeline depends on."""

from __future__ import annotations

import shutil
import subprocess

VERSION_FLAGS = {"blastn": "-version", "makeblastdb": "-version", "cutadapt": "--version"}


class ToolError(Exception):
    """An external program failed."""


def required_tools(cfg):
    tools = ["makeblastdb", "blastn"]
    if cfg.panel.design_length < cfg.panel.sequencing_length:
        tools.append("cutadapt")
    return tools


def tool_version(name):
    """The program's version text, or None if it is not installed."""
    path = shutil.which(name)
    if path is None:
        return None
    result = subprocess.run([path, VERSION_FLAGS[name]], capture_output=True, text=True)
    return (result.stdout or result.stderr).strip().splitlines()[0]


def run_tool(command, stdout_path=None):
    """Run a program, writing its output to the current log. Raises ToolError if it fails."""
    print("$ " + " ".join(command))
    result = subprocess.run(command, capture_output=True, text=True)
    if stdout_path is not None:
        with open(stdout_path, "w") as handle:
            handle.write(result.stdout)
    else:
        print(result.stdout, end="")
    print(result.stderr, end="")
    if result.returncode != 0:
        tail = result.stderr.strip()[-500:]
        raise ToolError(f"{command[0]} exited with code {result.returncode}: {tail}")
