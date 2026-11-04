"""Reading the layout of a MADC report without loading its read counts."""

from __future__ import annotations

from dataclasses import dataclass

# Columns that can come before the samples in a raw DArT MADC
METADATA_COLUMNS = {
    "AlleleID", "CloneID", "AlleleSequence", "ClusterConsensusSequence", "CallRate", "OneRatioRef",
    "OneRatioSnp", "FreqHomRef", "FreqHomSnp", "FreqHets", "PICRef", "PICSnp", "AvgPIC",
    "AvgCountRef", "AvgCountSnp", "RatioAvgCountRefAvgCountSnp",
}

# In a full layout the samples follow this column, whatever other tools call the statistics columns before it
LAST_METADATA_COLUMN = "RatioAvgCountRefAvgCountSnp"


@dataclass
class MadcLayout:
    first_sample_column: int  # 1-based
    samples: list
    filler_columns: int  # the number of leading * cells, or 0 if the file has no * rows


def read_layout(path):
    """Find the header row and work out which column the samples start in."""
    filler_columns = 0
    try:
        with open(path, encoding="utf-8-sig") as handle:
            for line in handle:
                cells = line.rstrip("\r\n").split(",")
                if cells[0] == "*" and not filler_columns:
                    filler_columns = cells.count("*")
                if cells[0] == "AlleleID":
                    for required in ("CloneID", "AlleleSequence"):
                        if required not in cells:
                            raise ValueError(f"{path}: header row has no {required} column")
                    if LAST_METADATA_COLUMN in cells:
                        first = cells.index(LAST_METADATA_COLUMN) + 1
                    else:
                        first = next((i for i, name in enumerate(cells) if name not in METADATA_COLUMNS), len(cells))
                    return MadcLayout(first + 1, cells[first:], filler_columns)
    except UnicodeDecodeError as err:
        raise ValueError(f"{path} is not a text file. Export the MADC as a .csv file, not Excel, "
                         "and decompress it if it is gzipped.") from err
    raise ValueError(f"{path}: no header row starting with AlleleID. Is this a MADC report?")


def longest_ref_read(path):
    """The length of the longest Ref allele in the report, or 0 if it has none."""
    longest = 0
    with open(path, encoding="utf-8-sig") as handle:
        for line in handle:
            cells = line.split(",", 3)
            if len(cells) > 2 and cells[0].endswith("|Ref"):
                longest = max(longest, len(cells[2]))
    return longest


def count_match_rows(path):
    """The number of alleles in the report other than Ref and Alt."""
    count = 0
    header_seen = False
    with open(path, encoding="utf-8-sig") as handle:
        for line in handle:
            if not header_seen:
                header_seen = line.startswith("AlleleID,")
            elif not line.split(",", 1)[0].endswith(("|Ref", "|Alt")):
                count += 1
    return count


def read_marker_ids(path):
    """The markers named in the report, taken from the part of each AlleleID before the |."""
    markers = set()
    header_seen = False
    with open(path, encoding="utf-8-sig") as handle:
        for line in handle:
            if not header_seen:
                header_seen = line.startswith("AlleleID,")
                continue
            markers.add(line.split(",", 1)[0].rsplit("|", 1)[0])
    return markers
