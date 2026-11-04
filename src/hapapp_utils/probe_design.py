"""Reading a DArT probe design file from Excel, comma-delimited or tab-delimited form."""

from __future__ import annotations

import csv
import re
import zipfile
from pathlib import Path

HEADER_NAMES = ("marker name", "markername")
MIN_COLUMNS = 7  # Marker Name through Marker_Type, the columns db00 reads
EXPECTED_COLUMNS = {3: "chrom", 4: "pos", 5: "variant"}
VARIANT = re.compile(r"\[[^\[\]/]+/[^\[\]/]+\]")


def _cell_text(value):
    if value is None:
        return ""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value).strip()


def _spreadsheet_tables(path):
    import openpyxl
    from openpyxl.utils.exceptions import InvalidFileException

    try:
        workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    except (InvalidFileException, zipfile.BadZipFile, KeyError, OSError) as err:
        raise ValueError(f"{path.name} could not be read as an Excel file. Save it again as .xlsx or .csv.") from err
    try:
        for sheet in workbook.worksheets:
            yield [[_cell_text(value) for value in row] for row in sheet.iter_rows(values_only=True)]
    finally:
        workbook.close()


def _text_table(path, delimiter):
    try:
        with open(path, encoding="utf-8-sig", newline="") as handle:
            return [[cell.strip() for cell in row] for row in csv.reader(handle, delimiter=delimiter)]
    except UnicodeDecodeError as err:
        raise ValueError(f"{path.name} is not a text file. Save it as .xlsx, .csv or tab-delimited .txt.") from err


def _check_header(name, header):
    cells = [cell.lower() for cell in header]
    if len(cells) < MIN_COLUMNS or not all(word in cells[index] for index, word in EXPECTED_COLUMNS.items()):
        raise ValueError(f"{name} does not look like a DArT probe design file. It should have the columns Marker Name, "
                         "Target Sequence, Reference Genome, Chrom, Chrom Pos, VariantAllelesDef and Marker_Type, "
                         f"but its header starts {header[:6]}.")


def read_probe_rows(path):
    """Return (header, marker rows), skipping any notes above the Marker Name header row."""
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".xls":
        raise ValueError(f"{path.name} is in the older Excel format. Save it as .xlsx or .csv.")
    if suffix in (".xlsx", ".xlsm"):
        tables = _spreadsheet_tables(path)
    elif suffix in (".txt", ".tsv"):
        tables = [_text_table(path, "\t")]
    else:
        tables = [_text_table(path, ",")]

    for table in tables:
        for index, row in enumerate(table):
            if row and row[0].lower() in HEADER_NAMES:
                _check_header(path.name, row)
                rows = [line for line in table[index + 1:] if line and line[0]]
                too_short = next((line[0] for line in rows if len(line) < MIN_COLUMNS), None)
                if too_short is not None:
                    raise ValueError(f"{path.name}: marker {too_short} has fewer than {MIN_COLUMNS} columns")
                return row, rows
    raise ValueError(f"{path.name}: no header row starting with Marker Name")


def probe_row_problem(row):
    """Why db00 could not use this marker's row, or None if it can."""
    if not row[4].isdigit():
        return f"{row[0]} has the position {row[4]!r}, which should be a whole number"
    if not VARIANT.fullmatch(row[5]):
        return f"{row[0]} has the variant {row[5]!r}, which should look like [A/T]"
    return None


def write_probe_csv(header, rows, path):
    """Write the plain CSV that db00 reads. Its line splitting cannot handle a comma inside a cell."""
    with open(path, "w") as out:
        for row in [["Marker Name"] + header[1:]] + rows:
            out.write(",".join(cell.replace(",", ";") for cell in row) + "\n")
