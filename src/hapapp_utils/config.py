"""Run settings: YAML loading, command-line overrides and validation."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import yaml

from . import __version__

REQUIRED = object()

# Every setting: dotted key -> (kind, default). Thresholds that decide which alleles get a
# permanent ID have no default, so a run never depends on a value nobody chose.
SCHEMA = {
    "project_id": ("str", REQUIRED),
    "inputs.madc": ("path", REQUIRED),
    "inputs.marker_lut": ("path", REQUIRED),
    "inputs.allele_db": ("path", REQUIRED),
    "inputs.match_counts": ("path", REQUIRED),
    "inputs.dup_tags": ("path", None),
    "inputs.previous_hap_pa": ("path", None),
    "panel.design_length": ("int", REQUIRED),
    "panel.sequencing_length": ("int", REQUIRED),
    "panel.first_sample_column": ("int", None),
    "panel.insert_size_column": ("str", None),
    "panel.adapter": ("str", "ACCGATCTG"),
    "filter.min_samples": ("int", REQUIRED),
    "filter.min_reads": ("int", REQUIRED),
    "blast.min_coverage": ("float", REQUIRED),
    "blast.min_identity": ("float", REQUIRED),
    "blast.allele_type_truth": ("str", "report"),
    "blast.threads": ("int", 6),
    "hap_pa.min_reads": ("int", 10),
    "hap_pa.descriptor": ("str", None),
    "output.dir": ("path", REQUIRED),
    "output.code_version": ("str", f"v{__version__.split('.')[0]}"),
}

BUILD_DB_SCHEMA = {
    "inputs.probe_design": ("path", REQUIRED),
    "inputs.madc": ("path", REQUIRED),
    "inputs.genome": ("path", REQUIRED),
    "inputs.chromosome_map": ("path", None),
    "database.name": ("str", REQUIRED),
    "database.allele_length": ("int", REQUIRED),
    "database.flank_length": ("int", 150),
    "database.dup_min_identity": ("float", 0.95),
    "database.dup_min_coverage": ("float", 0.99),
    "database.threads": ("int", 6),
    "output.dir": ("path", REQUIRED),
}

TEMPLATE = """\
# Settings for: hapapp-utils assign -c this_file.yaml
# Relative paths are resolved from the folder holding this file.

project_id: MY_PROJECT

inputs:
  madc: MY_PROJECT_MADC.csv
  marker_lut: marker_id_lookup.csv         # old marker ID, new marker ID, ...
  allele_db: species_allele_db_v001.fa
  match_counts: species_allele_db_v001_matchCnt_lut.txt
  # dup_tags: markers_to_remove.txt        # one marker ID per line, in the new ID format
  # previous_hap_pa: earlier_hapPA.csv     # presence/absence matrix to merge this project into

panel:
  design_length: 109
  sequencing_length: 109                   # adapters are trimmed if this exceeds design_length
  # first_sample_column: 17                # detected from the MADC header unless set here
  # insert_size_column: insert_size        # column in marker_lut, for mixed insert size panels
  # adapter: ACCGATCTG

filter:
  # Keep an allele if this many samples each have at least min_reads reads.
  # min_samples is lowered to 5% of the samples (rounded up) if that is smaller.
  min_samples: 2
  min_reads: 10

blast:
  min_coverage: 80
  min_identity: 85
  # allele_type_truth: report              # report or db, when the two disagree on allele class
  # threads: 6

hap_pa:
  # min_reads: 10                          # reads needed for an allele to count as present
  # descriptor: MY_PROJECT                 # column name in the matrix, defaults to project_id

output:
  dir: runs/MY_PROJECT                     # must be new or empty
  # code_version: v2                       # recorded in the report's Code_version column
"""

BUILD_DB_TEMPLATE = """\
# Settings for: hapapp-utils build-db -c this_file.yaml
# Relative paths are resolved from the folder holding this file.

inputs:
  probe_design: probe_design.xlsx          # the DArT template, as .xlsx, .csv or tab-delimited .txt
  madc: MY_PROJECT_MADC.csv                # any report from this panel
  genome: reference_genome.fa
  # chromosome_map: chromosome_names.csv   # probe name,genome name; if the genome names differ

database:
  name: species_allele_db                  # the first version is written as species_allele_db_v001.fa
  allele_length: 109                       # can be longer than the reads, to allow for longer reads later
  # flank_length: 150                      # bases fetched either side of each marker
  # dup_min_identity: 0.95                 # markers with Ref alleles this similar count as duplicates
  # dup_min_coverage: 0.99
  # threads: 6

output:
  dir: runs/core_db                        # must be new or empty
"""


class ConfigError(Exception):
    """A problem with the run settings, worded for the person running the tool."""


def _flatten(mapping, prefix=""):
    flat = {}
    for key, value in mapping.items():
        dotted = f"{prefix}{key}"
        if isinstance(value, dict):
            flat.update(_flatten(value, f"{dotted}."))
        else:
            flat[dotted] = value
    return flat


def _convert(key, kind, value, base):
    if kind == "path":
        if not isinstance(value, str):
            raise ConfigError(f"{key} must be a file path, got {value!r}")
        return (base / Path(value).expanduser()).resolve()
    if kind == "int":
        if isinstance(value, bool) or not isinstance(value, int):
            raise ConfigError(f"{key} must be a whole number, got {value!r}")
        return value
    if kind == "float":
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ConfigError(f"{key} must be a number, got {value!r}")
        return float(value)
    return str(value)


def _validate(values):
    for key in ("panel.design_length", "panel.sequencing_length", "panel.first_sample_column",
                "filter.min_samples", "filter.min_reads", "blast.threads", "hap_pa.min_reads",
                "database.allele_length", "database.flank_length", "database.threads"):
        if values.get(key) is not None and values[key] < 1:
            raise ConfigError(f"{key} must be at least 1, got {values[key]}")
    for key in ("blast.min_coverage", "blast.min_identity"):
        if key in values and not 0 < values[key] <= 100:
            raise ConfigError(f"{key} must be above 0 and at most 100, got {values[key]}")
    for key in ("database.dup_min_identity", "database.dup_min_coverage"):
        if key in values and not 0 < values[key] <= 1:
            raise ConfigError(f"{key} must be above 0 and at most 1, got {values[key]}")
    if values.get("blast.allele_type_truth", "report") not in ("report", "db"):
        raise ConfigError("blast.allele_type_truth must be report or db, "
                          f"got {values['blast.allele_type_truth']!r}")
    for key in ("project_id", "database.name"):
        if key in values and ("/" in values[key] or "\\" in values[key]):
            raise ConfigError(f"{key} cannot contain a slash, got {values[key]!r}")


def workflow_of(path):
    """Which workflow a settings file is for, going by its sections."""
    with open(path) as handle:
        try:
            raw = yaml.safe_load(handle)
        except yaml.YAMLError as err:
            raise ConfigError(f"Could not read {path}: {err}") from err
    return "build-db" if isinstance(raw, dict) and "database" in raw else "assign"


def load_config(path, overrides=(), outdir=None, schema=None):
    """Read a YAML file and apply overrides: command line first, then the file, then defaults."""
    schema = SCHEMA if schema is None else schema
    path = Path(path)
    if not path.is_file():
        raise ConfigError(f"Settings file not found: {path}")
    with open(path) as handle:
        try:
            raw = yaml.safe_load(handle)
        except yaml.YAMLError as err:
            raise ConfigError(f"Could not read {path}: {err}") from err
    if not isinstance(raw, dict):
        raise ConfigError(f"{path} must contain settings as key: value pairs")

    file_base = path.resolve().parent
    # A section whose settings are all commented out reads as an empty value
    section_names = {key.partition(".")[0] for key in schema if "." in key}
    supplied = {key: (value, file_base) for key, value in _flatten(raw).items()
                if not (key in section_names and value is None)}
    for item in overrides:
        key, separator, text = item.partition("=")
        if not separator:
            raise ConfigError(f"--set needs key=value, got {item!r}")
        supplied[key.strip()] = (yaml.safe_load(text), Path.cwd())
    if outdir is not None:
        supplied["output.dir"] = (str(outdir), Path.cwd())

    unknown = sorted(set(supplied) - set(schema))
    if unknown:
        raise ConfigError(f"Unknown setting(s): {', '.join(unknown)}. Valid settings: {', '.join(schema)}")

    values = {}
    missing = []
    for key, (kind, default) in schema.items():
        value, base = supplied.get(key, (None, file_base))
        if value is None:
            if default is REQUIRED:
                missing.append(key)
            values[key] = None if default is REQUIRED else default
        else:
            values[key] = _convert(key, kind, value, base)
    if missing:
        raise ConfigError(f"Missing required setting(s): {', '.join(missing)}")
    _validate(values)

    sections = {}
    for key, value in values.items():
        section, dot, name = key.partition(".")
        if dot:
            sections.setdefault(section, {})[name] = value
        else:
            sections[key] = value
    return SimpleNamespace(**{key: SimpleNamespace(**value) if isinstance(value, dict) else value
                              for key, value in sections.items()})


def config_to_dict(cfg, schema=None):
    """The settings as plain data, for saving alongside a run."""
    result = {}
    for key in SCHEMA if schema is None else schema:
        section, dot, name = key.partition(".")
        if dot:
            value = getattr(getattr(cfg, section), name)
            result.setdefault(section, {})[name] = str(value) if isinstance(value, Path) else value
        else:
            result[key] = getattr(cfg, key)
    return result
