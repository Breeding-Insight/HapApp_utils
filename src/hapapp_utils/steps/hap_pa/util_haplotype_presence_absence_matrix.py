#!/usr/bin/env python3

import os
import pandas as pd
import argparse


def get_project_haplotypes(madc, project_description, read_threshold):
    """
    Reads a MADC CSV file, converts allele read counts into presence/absence
    binary values based on a threshold, and outputs two summary CSV files.
    """

    print(f"[INFO] Loading data from {madc}")

    # 1. Load the data and clean up unnecessary columns
    df = pd.read_csv(madc, index_col='AlleleID')

    # Use errors='ignore' so it doesn't crash if these columns are already missing
    df = df.drop(['CloneID', 'AlleleSequence'], axis=1, errors='ignore')

    # 2. Coerce all remaining data to numeric values
    print("[INFO] Converting data to numeric formats")
    df = df.apply(pd.to_numeric, errors='coerce')

    # 3. Vectorized thresholding to avoid fragmented DataFrames
    # This evaluates the whole dataframe at once and converts True/False to 1/0
    print(f"[INFO] Applying read threshold of {read_threshold}")
    df = (df >= read_threshold).astype(int)

    # 4. Generate the first output file name safely using os.path
    base_name, ext = os.path.splitext(madc)
    outf_status = f"{base_name}_hapStatus.csv"

    print(f"[INFO] Writing haplotype status to {outf_status}")
    df.to_csv(outf_status, index_label='AlleleID')

    # 5. Filter for rows (AlleleIDs) that are present in at least one sample
    df_present = df[df.sum(axis=1) >= 1]

    # 6. Create the Presence/Absence (PA) summary DataFrame natively in Pandas
    # This avoids manual file writing and looping
    pa_df = pd.DataFrame({project_description: 1}, index=df_present.index)

    outf_pa = f"{base_name}_hapPA.csv"
    print(f"[INFO] Writing haplotype presence/absence summary to {outf_pa}")
    pa_df.to_csv(outf_pa, index_label='AlleleID')

    print("[INFO] Processing complete.")


if __name__ == '__main__':
    # Simplified argparse block at the bottom of the script
    parser = argparse.ArgumentParser(description="Assign haplotypes based on read count thresholds.")

    parser.add_argument('--madc', help='Path to the input MADC CSV file')
    parser.add_argument('--descriptor', help='Name of the PI/Project for the output summary in the format of P01_PI_project')
    parser.add_argument('--read_threshold', type=int, default=10,
                        help='Read count threshold for allele presence (default: 10)')

    args = parser.parse_args()

    get_project_haplotypes(args.madc, args.descriptor, args.read_threshold)