"""What every workflow shares: its run folder, log, saved settings and manifest."""

from __future__ import annotations

import contextlib
import hashlib
import json
import os
import platform
import sys
import time
import traceback
from datetime import datetime

import yaml

from .. import __version__
from ..config import config_to_dict
from ..tools import ToolError, tool_version

BLAST_COLUMNS = "6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue"


class PipelineError(Exception):
    """The run cannot continue, for a reason worded for the person running the tool."""


class CannotStart(PipelineError):
    """A problem found before the run began, so there is no process readme to read."""


class _Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, text):
        for stream in self.streams:
            stream.write(text)
        return len(text)

    def flush(self):
        for stream in self.streams:
            stream.flush()


@contextlib.contextmanager
def _working_directory(path):
    previous = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def call_step(step, argv):
    """Run a step's main(), turning any failure into a message that names the step."""
    name = step.__name__.rsplit(".", 1)[-1]
    try:
        step.main(argv)
    except SystemExit as err:
        if err.code not in (0, None):
            raise PipelineError(f"{name} stopped with code {err.code}: {argv}") from err
    except Exception as err:
        traceback.print_exc()
        raise PipelineError(f"{name} stopped with {type(err).__name__}: {err}") from err


def _write_readme_header(log, run, manifest, outdir):
    """The start of the process readme: which inputs and settings produced what follows."""
    print(time.strftime("%c"), file=log)
    print("\n --- Scripts, files, and parameters used ---", file=log)
    print(f"# hapapp-utils version: {manifest['hapapp_utils_version']}", file=log)
    print(f"# Command: {manifest['command']}", file=log)
    for tool, version in manifest["tools"].items():
        print(f"# Tool version: {tool} {str(version).removeprefix(f'{tool}: ')}", file=log)
    print(f"# Output folder: {outdir}", file=log)
    print(f"# Process README: {run.readme}", file=log)
    for line in run.readme_header():
        print(line, file=log)
    print(file=log)


def run_stages(cfg, run, stages, schema, tools, details=None, until=None, verbose=False, command=None, console=None):
    """Run a workflow's stages inside its output folder. Returns the saved manifest.

    run provides copy_inputs(), outputs(), counts_summary(), readme_header(), a readme path that
    a stage may change, and a stage_<name>() method for each stage. A stage may return a short note
    to show in place of "done".
    """
    console = console or sys.stdout
    outdir = cfg.output.dir
    try:
        outdir.mkdir(parents=True, exist_ok=True)
    except OSError as err:
        raise CannotStart(f"Cannot create the output folder {outdir}: {err.strerror}") from err
    (outdir / "resolved_config.yaml").write_text(yaml.safe_dump(config_to_dict(cfg, schema), sort_keys=False))
    manifest = {
        "hapapp_utils_version": __version__,
        "python": platform.python_version(),
        "command": command,
        "started": datetime.now().isoformat(timespec="seconds"),
        "status": "running",
        **(details or {}),
        "tools": {tool: tool_version(tool) for tool in tools},
        "inputs": {name: {"path": str(path), "sha256": _sha256(path)}
                   for name, path in vars(cfg.inputs).items() if path is not None},
    }
    manifest_path = outdir / "run_manifest.json"

    def save_manifest():
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    save_manifest()
    opened_as = run.readme
    (outdir / opened_as).parent.mkdir(parents=True, exist_ok=True)
    try:
        with open(outdir / opened_as, "w", buffering=1) as log:
            stream = _Tee(log, console) if verbose else log

            def say(text):
                print(text, file=console, flush=True)
                print(text, file=log)

            _write_readme_header(log, run, manifest, outdir)
            try:
                with _working_directory(outdir), contextlib.redirect_stdout(stream), contextlib.redirect_stderr(stream):
                    run.copy_inputs()
                    names = list(stages)
                    for number, name in enumerate(names, 1):
                        say(f"[{number}/{len(names)}] {stages[name]}")
                        start = time.monotonic()
                        note = getattr(run, f"stage_{name.replace('-', '_')}")()
                        say(f"      {note or 'done'} ({time.monotonic() - start:.0f} s)")
                        if name == until:
                            break
                    manifest["status"] = "finished" if until in (None, names[-1]) else f"stopped after {until}"
                    manifest["outputs"] = {key: value for key, value in run.outputs().items()
                                           if value and os.path.exists(value)}
                    manifest["counts"] = run.counts_summary()
            except KeyboardInterrupt:
                manifest["status"] = "interrupted"
                print("[STOPPED] Interrupted by the user", file=log)
                raise
            except Exception as err:
                manifest["status"] = "failed"
                manifest["error"] = str(err)
                if not isinstance(err, (PipelineError, ToolError)):
                    traceback.print_exc(file=log)
                print(f"[ERROR] {err}", file=log)
                raise
    finally:
        if run.readme != opened_as:
            os.replace(outdir / opened_as, outdir / run.readme)
        manifest["process_readme"] = run.readme
        manifest["finished"] = datetime.now().isoformat(timespec="seconds")
        save_manifest()
    return manifest
