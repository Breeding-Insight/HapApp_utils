#!/usr/bin/env python3

###############################################################################
###	 @author: Dongyan Zhao
###  @email: dongyan.zhao@ufl.edu
###  @date: 03/30/2026
###  Usage: db00_prep_lut_from_probeDesign.py [-h] [--madc MADC] probe
###############################################################################
"""
Prepare a SNP‑lookup table (LUT) from a probe design file.

Features:
• Optional MADC file – if omitted, LUT is built from all markers in the probe.
• APP’s output file automatically named: <probe>_snpID_lut.csv
• Simplified argparse at the bottom of this file.
• Small, well‑documented functions for easy debugging.
• Progress messages prefix with [INFO].
"""

from typing import Optional


# ----------------------------------------------------------------------
def get_panel_marker_ids(madc_path: str):
    """
    Read a MADC report and return a list of panel marker IDs.
    Lines beginning with '#', '*', or the header 'AlleleID' are ignored.
    """
    panel_markers = []
    with open(madc_path, 'r', encoding='utf-8') as fp:
        for line in fp:
            if line.startswith('#') or line.startswith('*') or \
               line.startswith('AlleleID'):
                continue
            panel_markers.append(line.strip().split(',')[1])
    print('[INFO] Number of markers in the panel:', len(panel_markers))
    print('[INFO] First 5 markers:', panel_markers[:5])
    return panel_markers


# ----------------------------------------------------------------------
def write_lut_header(out_fp: str, delimiter: str):
    """Write the LUT header, specialized for the expected CSV format."""
    header = [ "Panel_markerID", "Marker_ID", "Chr",
               "Pos", "Ref", "Alt", "Type",
               "Indel_pos", "Priority", "Note" ]
    with open(out_fp, 'w', encoding='utf-8') as fp:
        fp.write(delimiter.join(header) + '\n')
    return out_fp


# ----------------------------------------------------------------------
def prepare_lut(
        probe_path: str,
        panel_markers: Optional[list]
    ) -> int:
    """
    Build the LUT from <probe_path>.

    Returns the number of entries written.
    """
    out_path = probe_path.replace('.txt', '_snpID_lut.csv') \
        .replace('.csv', '_snpID_lut.csv')
    delimiter = '\t' if probe_path.endswith('.txt') else ','

    write_lut_header(out_path, delimiter)

    entry_count = 0
    with open(probe_path, 'r', encoding='utf-8-sig') as probe_fp,\
         open(out_path, 'a', encoding='utf-8') as out_fp:
        for line in probe_fp:
            if line.startswith('Marker Name') or line.startswith('#'):
                continue

            parts = line.strip().split(delimiter)
            marker_id = parts[0].strip()
            if panel_markers is not None and marker_id not in panel_markers:
                continue

            # Build Bi‑marker and positional fields
            bi_marker_id = f"{parts[3].strip()}_{parts[4].strip().zfill(9)}"
            chrom = parts[3].strip()
            pos   = parts[4].strip()

            # Variant definition: REF/ALT
            var_def   = parts[5].strip().strip('[]').split('/')
            ref, alt  = var_def[0].strip(), var_def[1].strip()

            var_type = parts[6].strip()
            priority = parts[7].strip()
            note     = parts[8].strip() if len(parts) > 8 else ''

            indel_pos = 'check_indel_pos' if 'indel' in var_type.lower() else ''

            # Write the row
            out_fp.write(
                delimiter.join([marker_id, bi_marker_id, chrom,
                                pos, ref, alt, var_type,
                                indel_pos, priority, note]) + '\n'
            )
            entry_count += 1

    print(f'[INFO] Prepared LUT "{out_path}" with {entry_count} entries')
    return entry_count


# ----------------------------------------------------------------------
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Build a SNP LUT from probe design file")
    parser.add_argument(
        "--madc", "-m",
        help="Raw MADC report (optional). If omitted, all markers in probe are used."
    )
    parser.add_argument("probe", help="Probe design file (.txt or .csv)")
    args = parser.parse_args()

    panel_markers = get_panel_marker_ids(args.madc) if args.madc else None
    prepare_lut(args.probe, panel_markers)