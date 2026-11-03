"""Reading a DArT probe design file from Excel, comma-delimited or tab-delimited form."""

from __future__ import annotations

import csv
from pathlib import Path

HEADER_NAMES = ("marker name", "markername")
MIN_COLUMNS = 7  # Marker Name through Marker_Type, the columns db00 reads


def _cell_text(value):
    if value is None:
        return ""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value).strip()


def _spreadsheet_tables(path):
    import openpyxl

    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    try:
        for sheet in workbook.worksheets:
            yield [[_cell_text(value) for value in row] for row in sheet.iter_rows(values_only=True)]
    finally:
        workbook.close()


def _text_table(path, delimiter):
    with open(path, encoding="utf-8-sig", newline="") as handle:
        return [[cell.strip() for cell in row] for row in csv.reader(handle, delimiter=delimiter)]


def read_probe_rows(path):
    """Return (header, marker rows), skipping any notes above the Marker Name header row."""
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in (".xlsx", ".xlsm"):
        tables = _spreadsheet_tables(path)
    elif suffix in (".txt", ".tsv"):
        tables = [_text_table(path, "\t")]
    else:
        tables = [_text_table(path, ",")]

    for table in tables:
        for index, row in enumerate(table):
            if row and row[0].lower() in HEADER_NAMES:
                rows = [line for line in table[index + 1:] if line and line[0]]
                too_short = next((line[0] for line in rows if len(line) < MIN_COLUMNS), None)
                if too_short is not None:
                    raise ValueError(f"{path.name}: marker {too_short} has fewer than {MIN_COLUMNS} columns")
                return row, rows
    raise ValueError(f"{path.name}: no header row starting with Marker Name")


def write_probe_csv(header, rows, path):
    """Write the plain CSV that db00 reads. Its line splitting cannot handle a comma inside a cell."""
    with open(path, "w") as out:
        for row in [["Marker Name"] + header[1:]] + rows:
            out.write(",".join(cell.replace(",", ";") for cell in row) + "\n")
