#!/bin/bash

###############################################################################
###	 @author: Dongyan Zhao
###  @email: dongyan.zhao@ufl.edu
###  @date: 03/30/2026
###  Usage: bash 20260320-build01_species_refAlt_coreDB.bash
##############################################################################

# Note: There are often IUPAC codes in Ref and Alt alleles from MADC report, therefore, cannot be used directly to build core DB.
#       Some DArTag panel creators also use IUPAC codes in the probe design file, therefore, cannot be used directly either.
#       The safest way to build core Ref/Alt DB is to use the reference genome sequences

# Required file
   # 1. Reference genome FASTA
   # 2. Chromosome lengths of the reference genome (format: chr 2300000)
   # 3. Panel probe design file
   # 4. One MADC report

# Workflow
   #    1). Prepare 180-300 bp flanking sequences of Ref and Alt alleles from probe design file submitted to DArT
   #    2). Extract amplicon sequences of Ref and Alt alleles from DArTag report
   #    3). BLAST Ref and Alt amplicon sequences to the 180-300 bp flanking sequences
   #    4). Determine alignment orientation between amplicon vs. flanking sequences, based on which update the flanking sequences to the same orientation of DArTag amplicons
   #    5). Rerun BLAST against the updated flanking sequences
   #    6). Get sfetch keys for the amplicons (because the 3' end of some amplicons are inaccurate)


# ===============
# Make changes based on custom file locations
# ===============
# Keep a readme file of the running process ######
exec &> ./HapApp_utils/data/species_allele_db_v001_process.readme

# scripts and files ######
SCRIPTS_DIR='./HapApp_utils/scripts/refAlt_coreDB'
PROBE_FILE='./HapApp_utils/data/flankseq/species_DArTag-probe-design.txt'
CHR_LEN='./HapApp_utils/data/ref/final_genome.fasta.len'
REF_GENOME='./HapApp_utils/data/ref/final_genome.fasta'
REPORT='./HapApp_utils/data/genotyping_report/species_MADC.csv'

ALLELE_DB_DIR='./HapApp_utils/data'
ALLELE_DB='species_allele_db_v001.fa'
MATCH_CNT='species_allele_db_v001_matchCnt_lut.txt'


# Parameters
REF_LEN=81 # Often to prepare Ref and Alt alleles longer than the amplicon length to accommodate for future longer amplicons
FLANK_LEN=100

# Start processing
now=$(date)
printf "%s\n" "$now"
printf "\n --- Scripts, files, and parameters used ---\n"
printf "# Report: $REPORT\n"


# =======
# Prepare flanking sequences from the reference genome
# =======
printf "\n# 1). Preparing marker ID LUT\n"
python3 $SCRIPTS_DIR/db00_prep_lut_from_probeDesign.py $PROBE_FILE
SNPID_LUT=${PROBE_FILE%????}'_snpID_lut.csv'


printf '\n# 2). Update snpIDs in DArTag report\n'
python3 $SCRIPTS_DIR/db02_update_snpID_in_madc_v1.py $SNPID_LUT $REPORT
REPORT_ID=${REPORT%????}'_snpID.csv'
REPORT_ID_MATCHCNT=${REPORT%????}'_snpID_matchCnt_lut.txt'


printf '\n# 3). Prepare sfetch coordinates of flanking sequences based on panel probe design file'
python3 $SCRIPTS_DIR/db01_get_reference_sfetch_keys_from_snpID_lut_v1.py \
--lut $SNPID_LUT \
--chr_len $CHR_LEN \
--flankBP $FLANK_LEN \
--madc $REPORT_ID
FLANK_SFETCH=${SNPID_LUT%????}'_f'$FLANK_LEN'bp_sfetchKeys.txt'


printf '\n# 4). Fetch flanking sequences'
FLANK_SEQ=${SNPID_LUT%????}'_f'$FLANK_LEN'bp_sfetchKeys.fa'
esl-sfetch -Cf $REF_GENOME $FLANK_SFETCH > $FLANK_SEQ



printf '\n# 5). Prepare Ref and Alt allele sequences - longer than the amplicon size\n'
python3 $SCRIPTS_DIR/db01_prep_ref_alt_flankSeq_from_lut_v1.py \
--snpID_lut $SNPID_LUT \
--flankSeq $FLANK_SEQ \
--flank_len $FLANK_LEN
REF_ALT_FLANK=${SNPID_LUT%????}'_f'$FLANK_LEN'bp_sfetchKeys_ref_alt.fa'
makeblastdb -in $REF_ALT_FLANK -dbtype nucl


printf '\n# 6). Extract amplicon sequences of Ref and Alt alleles from MADC report\n'
python3 $SCRIPTS_DIR/db03_ext_ref_alt_amp_from_madc_v1.py $REPORT
REF_ALT=${REPORT%????}'_ref_alt_amplicons.fa'



printf '\n# 7). BLAST Ref and Alt amplicon sequences from MADC report to flanking sequences\n'
REF_ALT_BLAST=${REPORT%????}'_snpID_ref_alt_amplicons.fa.f'$FLANK_LEN'bp.bn'
blastn -task blastn-short -dust no -soft_masking false \
-db $REF_ALT_FLANK \
-query $REF_ALT \
-out $REF_ALT_BLAST \
-evalue 1e-5 -num_threads 6 -max_target_seqs 10 \
-outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'
### 2025.3.21 I had to change the max_target_seqs to 10 because some marker are close to each other and the amplicons are similar 


printf '\n# 8). Determine alignment orientation between amplicon vs. flanking sequences,'
printf '# based on which update flanking sequences to the amplicon orientation\n'
python3 $SCRIPTS_DIR/db05_determine_alleleOri_from_blast_AND_update_f180bp_v1.py $REF_ALT_BLAST $REF_ALT_FLANK


printf '\n# 9). Rerun BLAST against the updated flanking sequences\n'
REF_ALT_FLANK_REV=${REF_ALT_FLANK%???}'_rev.fa'
makeblastdb -in $REF_ALT_FLANK_REV -dbtype nucl
REF_ALT_BLAST_REV=${REPORT%????}'_snpID_ref_alt_amplicons.fa.f'$FLANK_LEN'bp_rev.bn'
blastn -task blastn-short -dust no -soft_masking false \
-db $REF_ALT_FLANK_REV \
-query $REF_ALT \
-out $REF_ALT_BLAST_REV \
-evalue 1e-5 -num_threads 6 -max_target_seqs 10 \
-outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'


printf '\n# 10). Get sfetch keys for the amplicons (because the 3 end of some amplicons are inaccurate)\n'
python3 $SCRIPTS_DIR/db07_generate_ref_alt_sfetch_keys_from_blast_v1.1.py $SNPID_LUT $REF_ALT_BLAST_REV $REF_LEN


printf '\n# 11). Get the amplicon sequences from the flanking sequences\n'
esl-sfetch --index $REF_ALT_FLANK_REV
REF_ALT_BLAST_REV_SFETCH=${REPORT%????}'_snpID_ref_alt_amplicons.fa.f'$FLANK_LEN'bp_rev_'$REF_LEN'bp_sfetchKeys.txt'
SFETCH_FA=${REPORT%????}'_snpID_ref_alt_amplicons.fa.f'$FLANK_LEN'bp_rev_'$REF_LEN'bp_sfetch.fa'
esl-sfetch -Cf $REF_ALT_FLANK_REV $REF_ALT_BLAST_REV_SFETCH > $SFETCH_FA


printf '\n# 12). Check if there are identical Ref sequences from different marker loci'
printf '\n  # Get Ref allele sequences only'
SFETCH_FA_REF=${REPORT%????}'_snpID_ref_alt_amplicons.fa.f'$FLANK_LEN'bp_rev_'$REF_LEN'bp_sfetch_ref_only.fa'
seqkit grep -nrp '\|Ref_0001\b' \
$SFETCH_FA > $SFETCH_FA_REF


printf '\n  # Check sequence identity and output those with >95%% identity and >99%% coverage using MMseqs'
# --cov-mode 2: Both sides. Require both query and target coverage ≥ -c (i.e., min(qcov, tcov) ≥ -c).
# mmseqs easy-linclust [input.fa] [output_prefix] [tmp_dir]
SFETCH_FA_REF_MMSEQS=${REPORT%????}'_snpID_ref_alt_amplicons.fa.f'$FLANK_LEN'bp_rev_'$REF_LEN'bp_sfetch_ref_only.fa'
MMSEQS_OUT=${REPORT%????}'_snpID_ref_alt_amplicons.fa.f'$FLANK_LEN'bp_rev_'$REF_LEN'bp_sfetch_ref_only_dupTag'
TMP="$(dirname "$REPORT")/tmp"
mmseqs easy-linclust \
$SFETCH_FA_REF \
$MMSEQS_OUT \
$TMP \
--min-seq-id 0.95 -c 0.99 --cov-mode 2 -v 0 --threads 8

MMSEQS_OUT_ALLSEQS=${REPORT%????}'_snpID_ref_alt_amplicons.fa.f'$FLANK_LEN'bp_rev_'$REF_LEN'bp_sfetch_ref_only_dupTag_all_seqs.fasta'
MMSEQS_OUT_REPSEQS=${REPORT%????}'_snpID_ref_alt_amplicons.fa.f'$FLANK_LEN'bp_rev_'$REF_LEN'bp_sfetch_ref_only_dupTag_rep_seq.fasta'


printf '\n  # Extract highly similar non-self-pairing pairs'
MMSEQS_OUT_CLUSTER=${REPORT%????}'_snpID_ref_alt_amplicons.fa.f'$FLANK_LEN'bp_rev_'$REF_LEN'bp_sfetch_ref_only_dupTag_cluster.tsv'
MMSEQS_OUT_CLUSTER_NONSELF=${REPORT%????}'_snpID_ref_alt_amplicons.fa.f'$FLANK_LEN'bp_rev_'$REF_LEN'bp_sfetch_ref_only_dupTag_cluster_nonSelf.tsv'
awk '$1!=$2' $MMSEQS_OUT_CLUSTER > $MMSEQS_OUT_CLUSTER_NONSELF

printf '\n  # Check is the non-self-pairing file is empty'
if [ ! -s "$MMSEQS_OUT_CLUSTER_NONSELF" ]; then
    printf "# File is empty."
    rm $MMSEQS_OUT_CLUSTER_NONSELF

    printf '\n# 13) Copy allele fasta and match count lut to db directory\n'
    cp $SFETCH_FA $ALLELE_DB_DIR/$ALLELE_DB
    cp $MATCHCNT_LUT $ALLELE_DB_DIR/$MATCH_CNT


    printf '\n# 14) Make blastdb\n'
    makeblastdb -in $ALLELE_DB_DIR/$ALLELE_DB -dbtype nucl
else
    printf "\n  # File contains data."
    printf "\n  # Remove duplicate markers from the snpID_lut.csv\n"
    python3 $SCRIPTS_DIR/db09_rm_dupTags_from_LUT_and_db_v001.py \
    $SNPID_LUT \
    $MMSEQS_OUT_CLUSTER_NONSELF
    SNPID_LUT_RMDUPTAG=${PROBE_FILE%????}'_rmDupTag_snpID_lut.csv'

    printf "\n  # Remove duplicate markers from the matchCnt lut\n"
    python3 $SCRIPTS_DIR/db09_rm_dupTags_from_LUT_and_db_v001.py \
    $REPORT_ID_MATCHCNT \
    $MMSEQS_OUT_CLUSTER_NONSELF
    REPORT_ID_MATCHCNT_RMDUPTAG=${REPORT%????}'_snpID_matchCnt_lut_rmDupTag.txt'

    printf "\n  # Remove duplicate markers from the Ref and Alt allele sequence file\n"
    python3 $SCRIPTS_DIR/db09_rm_dupTags_from_LUT_and_db_v001.py \
    $SFETCH_FA \
    $MMSEQS_OUT_CLUSTER_NONSELF
    SFETCH_FA_RMDUPTAG=${REPORT%????}'_snpID_ref_alt_amplicons.fa.f'$FLANK_LEN'bp_rev_'$REF_LEN'bp_sfetch_rmDupTag.fa'

    printf '\n# 13) Copy allele fasta and match count lut to db directory\n'
    cp $SFETCH_FA_RMDUPTAG $ALLELE_DB_DIR/$ALLELE_DB
    cp $REPORT_ID_MATCHCNT_RMDUPTAG $ALLELE_DB_DIR/$MATCH_CNT

    printf '\n# 14) Make blastdb\n'
    makeblastdb -in $ALLELE_DB_DIR/$ALLELE_DB -dbtype nucl

    rm -rf $TMP
    rm -f $MMSEQS_OUT_ALLSEQS
    rm -f $MMSEQS_OUT_REPSEQS
fi




