def concat_madc(report_paths, outfile):
    """
    report_paths: list of MADC report file paths to merge (order preserved).
    #### CHANGE: refactored from taking a path to a file-of-paths into taking
    the already-resolved list of report paths directly, so the CLI layer can
    build that list from either input mode (a file-of-paths via
    --report_list_file, or a comma-separated --reports argument) without
    duplicating the merge logic.
    """
    import pandas as pd
    import os

    # --- OPTIMIZATION 1: Load all files into a list first ---
    dfs_to_join = []

    for curr_file in report_paths:
        curr_file = curr_file.strip()
        if not curr_file:
            continue

        df_temp = pd.read_csv(curr_file, index_col='AlleleID')

        # Extract order_id for the current file
        order_id = os.path.basename(curr_file).split("_")[0]
        print(f"Loading: {order_id}")

        # Append the order_id to all sample names (columns) in the current dataframe to ensure uniqueness
        #df_temp.columns = [f"{col}_{order_id}" for col in df_temp.columns]

        # Add to our list instead of joining immediately
        dfs_to_join.append(df_temp)

    if not dfs_to_join:
        print("Error: Input file list is empty.")
        return

    print("Merging dataframes...")
    # Perform a single outer join on all dataframes at once using concat axis=1
    df = pd.concat(dfs_to_join, axis=1, join='outer')

    # --- Sort and Format ---
    print("Processing and sorting groups...")
    df.fillna(0, inplace=True)
    df['CloneID'] = df.index.str.split('|', n=1).str[0]

    # Move CloneID to the first column
    first_col = df.pop('CloneID')
    df.insert(0, 'CloneID', first_col)

    # --- OPTIMIZATION 2 (rewritten): add missing Alt_0002 baseline rows and
    # order rows within each marker, both as single vectorized operations
    # instead of a Python-level loop over every CloneID group.
    #
    # #### CHANGE: the previous per-group loop built `reindex` from only
    # [Ref_0001, Alt_0002, sorted RefMatch..., sorted AltMatch...] and then
    # called clone_df.reindex(reindex) -- which silently DROPPED any row
    # whose label wasn't in that list (e.g. 'Other'-class alleles: 257 of
    # them in the real hemp MADC file). This rewrite keeps every row and
    # places unmatched/'Other' classes after AltMatch within their marker.

    # 1) bulk-add any missing Ref_0001 / Alt_0002 baseline rows per marker
    # (no per-marker loop).
    #
    # #### CHANGE: originally (and in my first rewrite) only Alt_0002 was
    # backfilled. But in a single order's own report Ref_0001 and Alt_0002
    # are always paired 1:1 (verified against the real hemp MADC file: 3040
    # markers, 3040 Ref_0001, 3040 Alt_0002) -- so whatever can cause a
    # marker's Alt_0002 to be absent from every order being merged (e.g. a
    # marker added to the panel/db after some orders ran) can equally cause
    # Ref_0001 to be absent. Both baseline anchors now get the same
    # protection, so every marker always has an explicit (possibly
    # all-zero) Ref_0001 and Alt_0002 row in the merged output.
    clone_ids = pd.Index(df['CloneID'].unique())
    for suffix in ('|Ref_0001', '|Alt_0002'):
        expected_ids = clone_ids + suffix
        missing_ids = expected_ids.difference(df.index)
        if len(missing_ids) > 0:
            print(f"  Adding {len(missing_ids)} missing {suffix.lstrip('|')} baseline rows.")
            missing_rows = pd.DataFrame(0, index=missing_ids, columns=df.columns)
            missing_rows['CloneID'] = missing_rows.index.str.split('|', n=1).str[0]
            df = pd.concat([df, missing_rows])

    # 2) bulk classify each row (Ref < Alt < RefMatch < AltMatch < everything
    # else, e.g. Other) and sort once across the whole table -- RefMatch/
    # AltMatch are listed before Ref/Alt in the regex alternation so they
    # match first (Ref/Alt would otherwise also match the RefMatch/AltMatch
    # prefix). Unmatched classes (Other, etc.) fall through to rank 4.
    allele_class = df.index.str.split('|', n=1).str[1].str.extract(
        r'^(RefMatch|AltMatch|Ref|Alt)', expand=False)
    class_rank = allele_class.map(
        {'Ref': 0, 'Alt': 1, 'RefMatch': 2, 'AltMatch': 3}).fillna(4)
    df['_class_rank'] = class_rank.values
    df['_allele_id'] = df.index  # secondary sort key = full AlleleID string,
                                  # same tie-break the old sorted() used

    print("Finalizing file...")
    df_ordered = df.sort_values(['CloneID', '_class_rank', '_allele_id'], kind='mergesort')
    df_ordered = df_ordered.drop(columns=['_class_rank', '_allele_id'])

    df_ordered.fillna(0, inplace=True)

    # Add a `Presence_count` column
    df_ordered['Presence_count'] = df_ordered.drop(columns=['CloneID']).sum(axis=1)

    # Reset the index and give it a column name
    df_ordered = df_ordered.reset_index().rename(columns={'index': 'AlleleID'})

    # --- NEW CODE: Convert numeric columns to integers ---
    # Drop the string columns to isolate just the numeric sample columns and Presence_count
    cols_to_convert = df_ordered.columns.drop(['AlleleID', 'CloneID'])
    # Convert those specific columns to integers
    df_ordered[cols_to_convert] = df_ordered[cols_to_convert].astype(int)
    # -----------------------------------------------------

    df_ordered.to_csv(outfile, index=False)
    print("Complete!")


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Generate stats from missing allele count")

    #### CHANGE: comma-separated report paths directly on the command line.
    parser.add_argument('--reports', default=None,
                        help='Comma-separated list of MADC report paths (e.g. "a.csv,b.csv,c.csv"). '
                             'Mutually exclusive with --report_list_file.')

    #### CHANGE: a file containing one report path per line, as the
    #### alternative input mode to --reports.
    parser.add_argument('--report_list_file', default=None,
                        help='A file containing MADC reports with allele name reformatted and unique sample names, '
                             'one file per line. Mutually exclusive with --reports.')

    parser.add_argument('--output',
                        help='Output file name')

    args = parser.parse_args()

    #### CHANGE: validate that exactly one input mode was given, then build
    #### the resolved list of report paths for concat_madc().
    if args.reports and args.report_list_file:
        parser.error("Provide either --reports or --report_list_file, not both.")
    if not args.reports and not args.report_list_file:
        parser.error("Provide either --reports or --report_list_file.")

    if args.reports:
        report_paths = [p.strip() for p in args.reports.split(',') if p.strip()]
    else:
        with open(args.report_list_file, 'r') as inp:
            report_paths = [line.strip() for line in inp if line.strip()]

    concat_madc(report_paths, args.output)
