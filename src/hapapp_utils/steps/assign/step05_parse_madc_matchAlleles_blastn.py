#!/usr/bin/python3

###############################################################################
###  @author: Dongyan Zhao
###  @email: dongyan.zhao@ufl.edu
###  @date: 03/30/2026
###  Usage: step05_parse_madc_allele81bp_blastn_v1.py [-h] [--read_col READ_COL] [--cov_threshold COV_THRESHOLD] [--iden_threshold IDEN_THRESHOLD] [--allele-type-ground-truth {db,report}] db_allele_cnt_inp db_allele_fasta report blast seq_len
###############################################################################

import re
import os
import csv


#### CHANGE (item 1): shared version-bump helper used by both the LUT and FASTA
#### updaters. Finds a "_v<digits>" tag anywhere in the filename (not just right
#### before the extension) and increments it, always zero-padding to 3 digits.
#### Previously generate_new_db_lut located the tag positionally
#### (split('_')[-3]) but substituted with a regex that only matched 3-digit
#### tags (r'v\d{3}') -- a filename using "v1" instead of "v001" would fail to
#### match, so the "updated" path was silently identical to the input path and
#### the LUT got overwritten in place.
def bump_version_tag(path):
    dirname, filename = os.path.split(path)
    m = re.search(r'_v(\d{1,3})(?=[_.]|$)', filename, flags=re.IGNORECASE)
    if m:
        new_tag = f"_v{int(m.group(1)) + 1:03d}"
        new_filename = filename[:m.start()] + new_tag + filename[m.end():]
    else:
        base, ext = os.path.splitext(filename)
        new_filename = f"{base}_v001{ext}"
    return os.path.join(dirname, new_filename)


def get_allele_class(allele_id):
    """
    Extract the LITERAL class token (Ref, Alt, RefMatch, AltMatch, Other) from
    an allele ID of the form 'marker|Class_suffix' or 'marker|Class_tmp_suffix'.
    Returns None if the ID doesn't contain '|'. Note this is NOT the same as
    the Ref/Alt "lineage" base used by get_locus_lineage_base() below --
    this returns the class exactly as written (e.g. 'RefMatch' stays
    'RefMatch', it is not normalized down to 'Ref').
    """
    if '|' not in allele_id:
        return None
    suffix = allele_id.split('|', 1)[1]
    return suffix.split('_', 1)[0]


def get_locus_lineage_bases(qid, sseqid):
    """
    Reproduces the ORIGINAL script's marker/lineage-base comparison (from its
    Case-2 branch) used to decide whether a temp allele is "at the same
    locus" as its best hit. Bug fix: the branch is chosen ONCE, based only on
    whether the TEMP allele's ID contains 'Match', and applied to derive BOTH
    bases -- not decided independently per ID. This matters for:
      - a RefMatch/AltMatch temp allele hitting a baseline Ref_0001/Alt_0002
        record: both sides normalize to the same 'marker|Ref' or 'marker|Alt'
        lineage base, so this is correctly treated as a MATCH (not a
        cross-lineage mismatch) -- RefMatch ~ Ref, AltMatch ~ Alt.
      - an 'Other'-class temp allele: both sides reduce to the bare marker
        regardless of what class the hit is, so 'Other' temp alleles are
        never flagged as lineage-mismatched, matching original behavior.
    (An earlier version of this helper computed each side's Match-branch
    independently, which broke both of the cases above.)
    """
    if 'Match' in qid:
        tmp_base = '_'.join(qid.rsplit('_')[:-2]).replace('Match', '')
        db_base = '_'.join(sseqid.rsplit('_')[:-1]).replace('Match', '')
    else:
        tmp_base = qid.rsplit('|', 1)[0]
        db_base = sseqid.rsplit('|', 1)[0]
    return tmp_base, db_base


def get_db_allele_fasta(db_allele_fasta, seq_len):
    db_fasta = {}
    seq_id = ''
    seq_chunks = []
    with open(db_allele_fasta) as inp:
        for raw in inp:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                # flush previous
                if seq_id.endswith('Ref_0001') or seq_id.endswith('Alt_0002'):
                    db_fasta[seq_id] = ''.join(seq_chunks)[:seq_len]
                seq_id = line[1:]  # drop '>'
                seq_chunks = []
            else:
                seq_chunks.append(line)
        # flush last
        if seq_id.endswith('Ref_0001') or seq_id.endswith('Alt_0002'):
            db_fasta[seq_id] = ''.join(seq_chunks)[:seq_len]
    return db_fasta


def get_db_allele_counts(db_allele_cnt_inp):
    inp = open(db_allele_cnt_inp)
    line = inp.readline()
    db_allele_cnt = {}
    while line:
        line_array = line.strip().split('\t')
        db_allele_cnt[line_array[0]] = int(line_array[1])
        line = inp.readline()
    inp.close()
    return db_allele_cnt


def get_tmp_rename_report(report):
    inp = open(report)
    tmp_rename_report = {}
    line = inp.readline()  # first data line
    while line:
        line_array = line.strip().split(',')
        tmp_rename_report[line_array[0]] = line_array[1:]
        line = inp.readline()
    inp.close()
    return tmp_rename_report


def get_unique_blast_hits(blast):
    """
    # BLAST results of only RefMatch and AltMatch alleles
    chr3.1_078721100|RefMatch_tmp_0001 81  1   59  chr3.1_078721100|AltMatch_0006 81   1   59      59    73  94.915  3.30e-21
    [qseqid                          qlen qstart qend             sseqid        slen sstart send length qcovs pident  evalue]
    [0                                 1     2     3               4               5   6      7   8        9      10  11]

    Keep the best hit per query where the LOCUS (marker) matches between
    query and subject -- the allele CLASS (Ref/Alt/RefMatch/AltMatch/Other)
    is intentionally NOT part of this filter (#### CHANGE item 2): a temp
    allele of one class is allowed to best-hit a db allele of a different
    class at the same marker. determine_allele_status() is responsible for
    deciding what to do when the classes disagree.

    #### CHANGE: best-hit ranking is now evalue-first (smallest evalue wins),
    #### not query-coverage-first. The old (qcov_len, pident, ...) ranking
    #### let a longer-but-lower-identity alignment beat a shorter alignment
    #### that was a PERFECT match, purely because it covered slightly more
    #### of the query -- e.g. a 100%-identity, 78bp alignment to the true
    #### Ref_0001 losing out to a 98.8%-identity, 81bp alignment to an
    #### unrelated AltMatch allele, just because 81 > 78. BLAST's own evalue
    #### already jointly weighs alignment length, identity, and statistical
    #### significance, so it's the more principled primary criterion.
    #### pident, qcov_len, and closeness of slen to qlen remain as
    #### tie-breakers (evalue can tie or bottom out at very strong hits).
    Best = smallest evalue, then max percent identity, then max query covered
    length, then smallest |slen-qlen|.
    Expects columns (0..11): qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue
    """
    blast_unique = {}
    with open(blast) as inp:
        for raw in inp:
            line = raw.strip()
            if not line:
                continue
            cols = line.split()
            qseqid = cols[0]
            qlen = int(cols[1])
            qstart, qend = int(cols[2]), int(cols[3])
            sseqid = cols[4]
            slen = int(cols[5])
            pident = float(cols[10])
            evalue = float(cols[11])

            #### CHANGE (item 2): base = marker only (everything before '|'),
            #### for both query and subject. No more class-aware branching.
            qbase = qseqid.split('|', 1)[0]
            sbase = sseqid.split('|', 1)[0]

            if qbase != sbase:
                continue

            qcov_len = abs(qend - qstart) + 1
            #### CHANGE: evalue-first score tuple (negated so the existing
            #### "higher tuple wins" comparison below still works -- lower
            #### evalue -> higher -evalue -> wins).
            score = (-evalue, pident, qcov_len, -abs(slen - qlen))

            prev = blast_unique.get(qseqid)
            if prev is None:
                blast_unique[qseqid] = (score, cols)
            else:
                if score > prev[0]:
                    # strict ">" -> first-seen hit wins on an exact tie
                    blast_unique[qseqid] = (score, cols)

    result = {k: v[1] for k, v in blast_unique.items()}
    print('\n[INFO] Number of unique BLAST queries written out:', len(result))
    return result


def determine_allele_status(db_allele_cnt, blast_unique, cov_threshold, iden_threshold,
                             tmp_rename_report, allele_type_ground_truth='report'):
    """
    Classify each temp allele's best BLAST hit as:
      - exact match to an existing allele of the SAME class (cov==100% & id==100% & slen==qlen)
          -> reuse that existing allele's ID,
      - exact match to an existing allele of a DIFFERENT (non-baseline) class
          -> handling controlled by `allele_type_ground_truth` (#### CHANGE item 3):
             'report' relabels the db allele under the temp's own class (one copy of the
             sequence is kept; the old ID is retired); 'db' trusts the db's existing
             label and just merges the temp allele into it, no relabeling,
      - exact match to a baseline Ref_0001/Alt_0002 -> drop with a warning (unchanged),
      - new allele at same locus (cov >= cov_threshold & id >= iden_threshold,
          including "exact" hits that are shorter than the matched db allele --
          #### CHANGE item 4: those no longer count as "exact"). This now also
          covers hits whose Ref/Alt LINEAGE disagrees with the temp allele's own
          (#### CHANGE, this round): previously discarded with a
          "Not aligned to the correct marker locus" warning, these are now kept
          and minted as new alleles; `allele_type_ground_truth` controls whether
          the new allele's class comes from the temp ('report') or from the
          best-hit db allele ('db'),
      - or discard (coverage/identity below threshold).

    Returns:
      updated_db_allele_cnt, dbVStmp_alleles_lut, tmpVSdb_alleles_lut, new_alleles,
      retired_db_ids, corrections, cross_lineage_rows
    """

    print('\n[INFO] Analyze BLAST results of the non-Ref/Alt alleles in this report:')
    print('  # allele_type_ground_truth =', allele_type_ground_truth)

    tmpVSdb_alleles_lut = {}   # temp -> fixed (existing or new)
    dbVStmp_alleles_lut = {}   # fixed -> [temp] (exact matches and optionally new)
    new_alleles = {}           # temp -> new fixed
    updated_db_allele_cnt = dict(db_allele_cnt)
    discarded = set()

    retired_db_ids = set()       # #### CHANGE item 3: old (wrong-class) db IDs to drop from the FASTA
    corrections = []             # #### CHANGE item 3: audit rows for the exact-match type-mismatch CSV
    relabel_cache = {}           # #### CHANGE item 3: sseqid -> corrected ID (exact-match case)
    cross_lineage_rows = []      # #### CHANGE (this round): audit rows for the Case-2 lineage-mismatch CSV

    cov_thr = float(cov_threshold)
    id_thr = float(iden_threshold)

    def assign_new(tmp_id, locus_class_base):
        # locus_class_base: 'marker|Class' -- the LUT key whose counter supplies
        # the next serial number. Callers choose this base according to
        # allele_type_ground_truth (temp's own class, or the best hit's class).
        updated_db_allele_cnt[locus_class_base] = updated_db_allele_cnt.get(locus_class_base, 0) + 1
        fixed = f"{locus_class_base}_{updated_db_allele_cnt[locus_class_base]:04d}"
        new_alleles[tmp_id] = fixed
        tmpVSdb_alleles_lut[tmp_id] = fixed
        return fixed

    for qid, cols in blast_unique.items():
        qlen = int(cols[1])
        qstart, qend = int(cols[2]), int(cols[3])
        sseqid = cols[4]
        pident = float(cols[10])
        slen = int(cols[5])

        qcov_len = abs(qend - qstart) + 1
        cov_len_thr = qlen * cov_thr / 100.0
        #### CHANGE (item 4): require slen == qlen too, so a temp allele that is a
        #### perfect but SHORTER substring of an existing (longer) db allele no
        #### longer counts as "exact" -- it falls through to the new-allele path
        #### below and is kept as its own distinct allele.
        is_exact = (qcov_len == qlen) and (pident == 100.0) and (slen == qlen)
        passes_thresholds = (qcov_len >= cov_len_thr) and (pident >= id_thr)

        # Case 1: exact 100%/100%/same-length match to an existing allele
        if is_exact:
            if sseqid.endswith('Ref_0001') or sseqid.endswith('Alt_0002'):
                print('  [WARNING] RefMatch or AltMatch matches Ref_0001 or Alt_0002, ignore:', cols)
                discarded.add(qid)
                continue

            q_class = get_allele_class(qid)
            s_class = get_allele_class(sseqid)

            if q_class is not None and s_class is not None and q_class != s_class:
                if allele_type_ground_truth == 'db':
                    #### CHANGE (this round): trust the db's existing class label;
                    #### merge the temp allele into it as-is, no relabeling.
                    dbVStmp_alleles_lut.setdefault(sseqid, []).append(qid)
                    tmpVSdb_alleles_lut[qid] = sseqid
                else:
                    #### CHANGE (item 3): the matched db allele is exactly this
                    #### sequence but is labeled under the WRONG class. Relabel it
                    #### under the temp allele's own class (LUT count + 1) instead
                    #### of keeping the sequence twice under two IDs.
                    if sseqid in relabel_cache:
                        corrected = relabel_cache[sseqid]
                    else:
                        corrected = assign_new(qid, qid.rsplit('_', 2)[0])
                        relabel_cache[sseqid] = corrected
                        retired_db_ids.add(sseqid)
                        seq = tmp_rename_report.get(qid, [None, None])[1] if qid in tmp_rename_report else None
                        corrections.append({
                            'db_allele_ID': sseqid,
                            'db_sequence': seq,   # identical to temp's sequence by construction (100%/100% hit)
                            'report_allele_ID': qid,
                            'report_sequence': seq,
                            'corrected_db_allele_ID': corrected,
                        })
                    dbVStmp_alleles_lut.setdefault(corrected, []).append(qid)
                    tmpVSdb_alleles_lut[qid] = corrected
            else:
                dbVStmp_alleles_lut.setdefault(sseqid, []).append(qid)
                tmpVSdb_alleles_lut[qid] = sseqid
            continue

        # Case 2: passes coverage and identity thresholds -> create a NEW allele at same locus
        elif passes_thresholds:
            tmp_lineage_base, db_lineage_base = get_locus_lineage_bases(qid, sseqid)

            if tmp_lineage_base == db_lineage_base:
                fixed = assign_new(qid, qid.rsplit('_', 2)[0])
                dbVStmp_alleles_lut.setdefault(fixed, []).append(qid)
            else:
                #### CHANGE (this round): previously this branch printed
                #### "[WARNING] Not aligned to the correct marker locus" and
                #### discarded qid entirely. These are near-exact (threshold-
                #### passing but not 100%/100%) hits at the correct MARKER whose
                #### Ref/Alt lineage differs from the temp allele's own -- they
                #### are now kept and minted as new alleles instead.
                is_baseline_hit = sseqid.endswith('Ref_0001') or sseqid.endswith('Alt_0002')
                if allele_type_ground_truth == 'db' and not is_baseline_hit:
                    locus_class_base = sseqid.rsplit('_', 1)[0]   # marker|HitClass
                else:
                    # 'report' ground truth, OR the best hit is a baseline record
                    # (never mint a pseudo Ref_000N/Alt_000N off of a baseline hit)
                    locus_class_base = qid.rsplit('_', 2)[0]      # marker|TempClass
                fixed = assign_new(qid, locus_class_base)
                dbVStmp_alleles_lut.setdefault(fixed, []).append(qid)

                print(f'  [INFO] Cross-lineage near-match at correct marker (kept as new allele, '
                      f'ground_truth={allele_type_ground_truth}):', cols)
                cross_lineage_rows.append({
                    'report_allele_ID': qid,
                    'report_class': get_allele_class(qid),
                    'report_sequence': tmp_rename_report.get(qid, [None, None])[1] if qid in tmp_rename_report else None,
                    'best_hit_db_ID': sseqid,
                    'best_hit_class': get_allele_class(sseqid),
                    'pident': pident,
                    'qcov_pct': round(100.0 * qcov_len / qlen, 3),
                    'new_allele_ID': fixed,
                })
            continue

        # Case 3: discard (either coverage or identity failed)
        else:
            if qcov_len >= cov_len_thr and pident < id_thr:
                print(f'[INFO] >={cov_threshold}% coverage and <{iden_threshold}% identity', cols)
            discarded.add(qid)
            continue

    print('  # Number of non-Ref/Alt alleles with fixed ID assigned:', len(tmpVSdb_alleles_lut))
    print('    - Number already existing in the database:', len(tmpVSdb_alleles_lut) - len(new_alleles))
    print('    - Number of NEW alleles found in this report:', len(new_alleles))
    print('  # Number of alleles discarded because of low BLAST sequence coverage/identity:', len(discarded))
    if corrections:
        print('  # Number of db alleles relabeled to the temp allele\'s class (exact-match type correction):',
              len(relabel_cache))
        print('    - Number of temp alleles involved:', len(corrections))
    if cross_lineage_rows:
        print('  # Number of cross-lineage near-matches kept as new alleles (previously discarded):',
              len(cross_lineage_rows))

    return (updated_db_allele_cnt, dbVStmp_alleles_lut, tmpVSdb_alleles_lut, new_alleles,
            retired_db_ids, corrections, cross_lineage_rows)


def generate_report_with_fixed_alleleID(report, tmp_rename_report, dbVStmp_alleles_lut, tmpVSdb_alleles_lut, new_alleles, db_fasta, read_col):
    print('\n[INFO] Assign fixed allele IDs to non-Ref/Alt alleles in this report:')
    import pandas as pd
    dup = []
    dup_count = 0
    uniq_count = 0
    cols = tmp_rename_report['AlleleID']
    # tmp_rename_report: {alfalfaRep2vsXJDY1_shared_1029546|AltMatch_tmp_0001: [alfalfaRep2vsXJDY1_shared_1029546,CTTTCAGGATTGTCGATTTCCAAGCTGTTAGATTCACCACAGTGCATAATTAAAGTACTTCAAAACCACCAAATTTTAAAA,3230,0,25...]
    # dbVStmp_alleles_lut: {'chr3.1_078721100|RefMatch_0002': ['chr3.1_078721100|RefMatch_tmp_0003'] ...}
    for i in dbVStmp_alleles_lut:
        if len(dbVStmp_alleles_lut[i]) > 1:
            same_allele_seq = {}
            for j in dbVStmp_alleles_lut[i]:
                if j in tmp_rename_report:
                    same_allele_seq[j] = tmp_rename_report[j][:2] + list(map(float, tmp_rename_report[j][2:]))
                    dup.append([i, j] + tmp_rename_report[j])
                else:
                    print(' [INFO] This allele does not exist in temporary report: ', j)

            if not same_allele_seq:
                continue

            dup_count += len(same_allele_seq)
            uniq_count += 1

            #### CHANGE (item 5, bug A): pick the metadata (CloneID/AlleleSequence)
            #### deliberately from the FIRST present temp allele, not whichever
            #### happened to be last in the loop. Numeric (read-count) columns
            #### are summed across the whole group regardless of order.
            canonical_id = next(iter(same_allele_seq))
            meta_info = tmp_rename_report[canonical_id][:2]
            df = pd.DataFrame.from_dict(same_allele_seq, orient='index', columns=cols)
            combined = meta_info + list(map(str, df.sum()[2:]))

            for j in same_allele_seq:
                del tmp_rename_report[j]
                del tmpVSdb_alleles_lut[j]

            tmp_rename_report[canonical_id] = combined
            #### CHANGE (item 5, bug B): restore the fixed-ID mapping for the
            #### surviving row. Previously this was deleted in the loop above
            #### and never re-added, so the combined row silently fell into
            #### "allele_discarded" in the write loop below and never made it
            #### into the output report.
            tmpVSdb_alleles_lut[canonical_id] = i
            dup.append([i, 'combined'] + combined)
        else:
            pass

    if dup_count > 0:
        print('[WARNING] Number of non-Ref/Alt alleles having the same sequences:')
        print('   # Number of non-Ref/Alt alleles aligned to the same alleles in the database: ', dup_count)
        print('   # Number of non-Ref/Alt alleles after COMBINING those with the same sequences: ', uniq_count)

    #### CHANGE (item 6): consolidate the two near-identical output-path branches
    # DAl22-7011_MADC_Report_Part1_tmp_rename_cutadapt.csv
    if 'updatedSeq' in report:
        report_outf = report.replace('tmp_rename_updatedSeq', 'rename_updatedSeq')
        dup_outf = report.replace('tmp_rename_updatedSeq', 'rename_updatedSeq_dup')
    else:
        report_outf = report.replace('tmp_rename', 'rename')
        dup_outf = report.replace('tmp_rename', 'rename_dup')

    #### CHANGE (item 6): dup file handle is now always closed (context manager)
    if len(dup) > 0:
        with open(dup_outf, 'w') as outp_dup:
            outp_dup.write('AlleleID' + ',' + ','.join(tmp_rename_report['AlleleID']) + '\n')
            for i in dup:
                outp_dup.write(','.join(i) + '\n')

    new_alleles_fasta = {}
    allele_cnt = 0
    allele_discarded = []
    with open(report_outf, 'w') as outp_report:
        for i in tmp_rename_report:
            if i == 'AlleleID':
                outp_report.write(i + ',' + ','.join(tmp_rename_report[i][:2] + tmp_rename_report[i][read_col - 2:]) + '\n')
            elif i.endswith('Ref_0001') or i.endswith('Alt_0002'):
                outp_report.write(i + ',' + tmp_rename_report[i][0] + ',' + db_fasta[i].upper() + ',' + ','.join(tmp_rename_report[i][read_col - 2:]) + '\n')
            else:
                if i in tmpVSdb_alleles_lut:
                    outp_report.write(tmpVSdb_alleles_lut[i] + ',' + tmp_rename_report[i][0] + ',' + tmp_rename_report[i][1].upper() + ',' + ','.join(tmp_rename_report[i][read_col - 2:]) + '\n')
                    allele_cnt += 1
                else:
                    allele_discarded.append(i)

                # new_alleles = {'chr3.1_078721100|AltMatch_tmp_001': 'chr3.1_078721100|AltMatch_012', ...}
                if i in new_alleles:
                    if '>' + new_alleles[i] not in new_alleles_fasta:
                        new_alleles_fasta['>' + new_alleles[i]] = tmp_rename_report[i][1].upper()
                    else:
                        print('[WARNING] Allele already exists: ', i, new_alleles[i])
                else:
                    pass

    print('\n[INFO] Report with fixed allele IDs assigned to non-Ref/Alt alleles:')
    print('  # Number of non-Ref/Alt alleles assigned fixed IDs: ', allele_cnt)
    print('  # Number of non-Ref/Alt alleles discarded because of low coverage/identity or same sequences: ', len(allele_discarded), '\n\n')
    return new_alleles_fasta


def write_type_correction_csv(db_allele_fasta, corrections):
    """#### CHANGE (item 3): audit CSV for exact-match db alleles relabeled to the temp allele's class.
    #### CHANGE (this round): named off the db FASTA path (with a distinguishing
    #### suffix), matching write_cross_lineage_csv's convention, rather than off
    #### the blast path.
    """
    if not corrections:
        return None
    outf = db_allele_fasta.replace('.fa', "") + '_exact_match_type_corrections.csv'
    fieldnames = ['db_allele_ID', 'db_sequence', 'report_allele_ID', 'report_sequence', 'corrected_db_allele_ID']
    with open(outf, 'w', newline='') as cf:
        writer = csv.DictWriter(cf, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(corrections)
    print('\n[INFO] Wrote exact-match db-allele type-mismatch corrections to:', outf)
    print('  # Number of correction rows:', len(corrections))
    return outf


def write_cross_lineage_csv(db_allele_fasta, cross_lineage_rows):
    """#### CHANGE (this round): audit CSV for Case-2 lineage-mismatched hits kept as new alleles."""
    if not cross_lineage_rows:
        return None
    outf = db_allele_fasta.replace('.fa', "") + '_mismatchedType_alleles.csv'
    fieldnames = ['report_allele_ID', 'report_class', 'report_sequence', 'best_hit_db_ID',
                  'best_hit_class', 'pident', 'qcov_pct', 'new_allele_ID']
    with open(outf, 'w', newline='') as cf:
        writer = csv.DictWriter(cf, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(cross_lineage_rows)
    print('\n[INFO] Wrote db vs. report type-mismatched alleles minting log to:', outf)
    print('  # Number of rows:', len(cross_lineage_rows))
    return outf


def generate_new_db_lut(db_alleleCnt_lut_file, updated_db_allele_cnt, db_allele_cnt):
    #### CHANGE (item 1): use the shared, robust version-bump helper
    outf = bump_version_tag(db_alleleCnt_lut_file)
    with open(outf, 'w') as outp_lut:
        for i in updated_db_allele_cnt.keys():
            if i in db_allele_cnt:
                outp_lut.write('\t'.join([i, str(updated_db_allele_cnt[i]), str(db_allele_cnt[i])]) + '\n')
            else:
                outp_lut.write('\t'.join([i, str(updated_db_allele_cnt[i]), '0']) + '\n')

    print('[INFO] Update match allele cnt LUT with new allele counts for RefMatch and AltMatch:')
    print('  # Existing allele COUNT database: ', db_alleleCnt_lut_file)
    print('  # Updated allele COUNT database: ', outf)


def update_db_allele_fasta(db_allele_fasta, new_alleles_fasta, retired_db_ids=None):
    """
    Merge new allele sequences into the existing FASTA, drop any db IDs that
    were relabeled (#### CHANGE item 3: retired_db_ids -- keeps exactly one
    copy of a relabeled sequence, under its corrected ID), write an updated
    FASTA with a bumped version, and sort sequences per marker in this order:
      Ref_* (Ref_0001 first),
      Alt_* (Alt_0002 first),
      RefMatch_*,
      AltMatch_*,
      Other_*,
    with numeric suffix ascending within each class.
    # Report how many sequences are in the old and new FASTA.
    """
    if retired_db_ids is None:
        retired_db_ids = set()

    def parse_fasta(path):
        records = {}
        seq_id = None
        seq_chunks = []
        with open(path) as inp:
            for raw in inp:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if seq_id is not None:
                        records[seq_id] = "".join(seq_chunks)
                    seq_id = line[1:]  # store without '>'
                    seq_chunks = []
                else:
                    seq_chunks.append(line)
            if seq_id is not None:
                records[seq_id] = "".join(seq_chunks)
        return records

    def normalize_new(new_dict):
        out = {}
        for k, v in new_dict.items():
            key = k[1:] if isinstance(k, str) and k.startswith(">") else k
            out[key] = v
        return out

    class_rank = {"Ref": 0, "Alt": 1, "RefMatch": 2, "AltMatch": 3, "Other": 4}
    unknown_rank = len(class_rank)

    def header_sort_key(header):
        try:
            marker, cls = header.split("|", 1)
        except ValueError:
            return (header, unknown_rank, float("inf"))
        parts = cls.split("_", 1)
        cls_name = parts[0]
        try:
            num = int(parts[1]) if len(parts) > 1 else float("inf")
        except ValueError:
            num = float("inf")
        rank = class_rank.get(cls_name, unknown_rank)
        return (marker, rank, num)

    # Load, drop retired (relabeled) IDs, normalize, merge
    existing = parse_fasta(db_allele_fasta)
    n_retired = 0
    for rid in retired_db_ids:
        if existing.pop(rid, None) is not None:
            n_retired += 1

    normalized_new = normalize_new(new_alleles_fasta)
    merged = dict(existing)
    merged.update(normalized_new)

    #### CHANGE (item 1): shared version-bump helper
    outf = bump_version_tag(db_allele_fasta)

    with open(outf, "w") as outp:
        for header in sorted(merged.keys(), key=header_sort_key):
            outp.write(f">{header}\n{merged[header]}\n")

    print("\n[INFO] Update allele SEQUENCE database (FASTA) with new non-Ref/Alt alleles:")
    print("  # Existing allele SEQUENCE database:", db_allele_fasta)
    print("    - Number of allele sequences in existing database:", len(existing) + n_retired)
    if n_retired:
        print("    - Number of relabeled (type-corrected) entries removed:", n_retired)
    print("  # Updated allele SEQUENCE database:", outf)
    print("    - Number of allele sequences in the updated database:", len(merged))


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('db_allele_cnt_inp',
                         help='Counts of RefMatch and AltMatch in allele db')

    parser.add_argument('db_allele_fasta',
                         help='Current fasta db containing existing alleles, where new alleles will be added into')

    parser.add_argument('report',
                         help='DArTag report with alleles assigned temporary names and allele sequences with adapter trimmed')

    parser.add_argument('blast',
                         help='BLASTN results of the allele sequences to the allele DB')

    parser.add_argument('seq_len', type=int,
                         help='The maximum length of amplicons in MADC')

    parser.add_argument('--read_col', type=int, default=17,
                         help='Column number of read counts in the DArTag report')

    parser.add_argument('--cov_threshold', type=float, default=90.0,
                         help='Minimum coverage threshold for BLASTN hits to be considered as RefMatch/AltMatch')

    parser.add_argument('--iden_threshold', type=float, default=90.0,
                         help='Minimum identity threshold for BLASTN hits to be considered as RefMatch/AltMatch')

    #### CHANGE (this round): new argument controlling which side's class label
    #### is authoritative when a temp allele's class disagrees with the class
    #### of the db allele it best matches.
    parser.add_argument('--allele-type-ground-truth', dest='allele_type_ground_truth',
                         choices=['db', 'report'], default='report',
                         help="Which side's class label wins on a class disagreement: "
                              "'report' trusts the class encoded in the DArTag report/temp "
                              "allele ID (default); 'db' trusts the class already recorded "
                              "in the allele database.")

    args = parser.parse_args()

    import datetime
    now = datetime.datetime.now()
    nowf = now.strftime("%Y-%m-%d, %H:%M:%S")
    print('\n---------\n[INFO] Processed on ' + nowf)
    print('  # DArTag report with alleles assigned temporary names: ' + args.report)
    print('  # Counts of non-Ref/Alt alleles: ' + args.db_allele_cnt_inp)
    print('  # Current fasta db containing existing alleles: ' + args.db_allele_fasta)

    #### CHANGE (item 6): removed dead "if args.cov_threshold is None" branch
    #### (the argument always has a default, so that branch could never fire)
    print('  # Coverage and identity thresholds for BLASTN hits to be retained:')
    print('    - Minimum coverage threshold: ', args.cov_threshold)
    print('    - Minimum identity threshold: ', args.iden_threshold)
    print('  # Allele-type ground truth on class disagreement: ', args.allele_type_ground_truth)

    db_allele_cnt = get_db_allele_counts(args.db_allele_cnt_inp)
    db_fasta = get_db_allele_fasta(args.db_allele_fasta, args.seq_len)
    tmp_rename_report = get_tmp_rename_report(args.report)
    blast_unique = get_unique_blast_hits(args.blast)

    #### CHANGE (item 6): unique-hits file handle now closed via context manager
    unique_outf = args.blast + '.unique.csv'
    with open(unique_outf, 'w') as outp:
        for key, value in blast_unique.items():
            outp.write(','.join(value) + '\n')

    #### CHANGE (item 3 / this round): determine_allele_status now takes
    #### tmp_rename_report and allele_type_ground_truth, and returns two extra
    #### values: retired_db_ids/corrections (exact-match) and cross_lineage_rows
    #### (Case-2 lineage-mismatch, kept as new alleles instead of discarded)
    (updated_db_allele_cnt, dbVStmp_alleles_lut, tmpVSdb_alleles_lut,
     new_alleles, retired_db_ids, corrections, cross_lineage_rows) = determine_allele_status(
        db_allele_cnt, blast_unique, args.cov_threshold, args.iden_threshold,
        tmp_rename_report, args.allele_type_ground_truth)

    write_type_correction_csv(args.db_allele_fasta, corrections)
    write_cross_lineage_csv(args.db_allele_fasta, cross_lineage_rows)

    new_alleles_fasta = generate_report_with_fixed_alleleID(
        args.report, tmp_rename_report, dbVStmp_alleles_lut, tmpVSdb_alleles_lut,
        new_alleles, db_fasta, args.read_col)

    if len(new_alleles.keys()) > 0:
        generate_new_db_lut(args.db_allele_cnt_inp, updated_db_allele_cnt, db_allele_cnt)
        #### CHANGE (item 3): pass retired_db_ids through so relabeled entries
        #### are dropped from the FASTA rather than duplicated
        update_db_allele_fasta(args.db_allele_fasta, new_alleles_fasta, retired_db_ids)
    else:
        print('\n[INFO] # No new alleles found in this project\n')
