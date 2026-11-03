"""What every workflow shares: its run folder, log, saved settings and manifest."""

from __future__ import annotations

import contextlib
import hashlib
import json
import os
import platform
import sys
import time
from datetime import datetime

import yaml

from .. import __version__
from ..config import config_to_dict
from ..tools import tool_version

BLAST_COLUMNS = "6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue"


class PipelineError(Exception):
    """The run cannot continue, for a reason worded for the person running the tool."""


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
    """Run a step's main(), turning its argument errors into a normal failure."""
    try:
        step.main(argv)
    except SystemExit as err:
        if err.code not in (0, None):
            raise PipelineError(f"{step.__name__} stopped with code {err.code}: {argv}") from err


def run_stages(cfg, run, stages, schema, tools, details=None, until=None, verbose=False, command=None, console=None):
    """Run a workflow's stages inside its output folder. Returns the saved manifest.

    run provides copy_inputs(), outputs(), counts_summary() and a stage_<name>() method for each
    stage. A stage may return a short note to show in place of "done".
    """
    console = console or sys.stdout
    outdir = cfg.output.dir
    outdir.mkdir(parents=True, exist_ok=True)
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
    with open(outdir / "pipeline.log", "w", buffering=1) as log:
        stream = _Tee(log, console) if verbose else log

        def say(text):
            print(text, file=console, flush=True)
            print(text, file=log)

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
        except Exception as err:
            manifest["status"] = "failed"
            manifest["error"] = str(err)
            print(f"[ERROR] {err}", file=log)
            raise
        finally:
            manifest["finished"] = datetime.now().isoformat(timespec="seconds")
            save_manifest()
    return manifest
