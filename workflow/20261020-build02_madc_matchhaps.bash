#!/bin/bash

###############################################################################
###	 @author: Dongyan Zhao
###  @email: dongyan.zhao@ufl.edu
###  @date: 10/20/2026
###  Usage: bash 20260320-build02_species_madc_haps.bash
###  Base db was already generated for REF and ALT alleles
###############################################################################


######################################################################
# ###### Note: USER INPUT information ######
# The readme file should be one version above the current db version
# If no new alleles are found, give the readme file another name

#### 1. Parameters and files for fixed allele ID assignment
PROCESS_README='/Users/dongyan.zhao/PycharmProjects/BI/chestnut_dartag_v2_00_microhaplotype_db/data/DCnutv2_allele_db_v003_process.readme'
NO_NEW_ALLELE_README='/Users/dongyan.zhao/PycharmProjects/BI/chestnut_dartag_v2_00_microhaplotype_db/data/DCnutv2_allele_db_v002_DCnut26-11629_noNewAllele.readme'
exec > >(tee "$PROCESS_README") 2>&1

# Scripts and files
CUTADAPT_PATH="/Library/Frameworks/Python.framework/Versions/3.9/bin/cutadapt"
SCRIPTS_DIR='/Users/dongyan.zhao/PycharmProjects/BI/HapApp_utils/scripts/RefMatch_AltMatch_Other'
CODE_VER='v1'

MARKERID_LUT='/Users/dongyan.zhao/Downloads/hap/chestnut_DArTag-probe-design-v2_type_snpID_lut.csv'
ALLELE_DB_DIR='/Users/dongyan.zhao/PycharmProjects/BI/chestnut_dartag_v2_00_microhaplotype_db/data'
ALLELE_DB='DCnutv2_allele_db_v002.fa'
MATCHCNT_LUT='DCnutv2_allele_db_v002_matchCnt_lut.txt'
REPORT='/Users/dongyan.zhao/Downloads/hap/DCnut26-11629_MADC.csv'

# First sample column in the MADC file, which is the first column of sample read counts
FIRST_SAMPLE_COL=17

# If the panel is a mixed insert size panel,
# the column name in the markerID_lut that contains the insert size information
# If left empty, the script will not split the MADC file by insert size and will process all microhaplotypes together
LUT_INSERT_SIZE_COL='insert_size'

# Initial filter of microhaplotypes with missing data
MIN_SAMPLES=1
MIN_READS=10

# Marker panel information
DESIGN_LEN=109
SEQ_LEN=109
# For microhaplotypes that bear different type assignments between the db and the report,
# which type should be used as the ground truth for the BLAST filtering step? Options: "db" or "report"
ALLELE_TYPE="report"

# BLAST thresholds
COV=80
IDEN=80

#### 2. Parameters and files for generating haplotype presence-absence matrix
HAPPA_SCRIPT_DIR="/Users/dongyan.zhao/PycharmProjects/BI/HapApp_utils/scripts/hapPA"
PROJECT_DESCRIPTOR='P04_DCnut26-11629'
MIN_READS_FOR_PRESENCE=10

LATEST_HAPPA_MATRIX="${ALLELE_DB_DIR}/hapSearch/hapPA/${PROJECT_DESCRIPTOR}_hapPA.csv"


# ###### End of user input information ######
######################################################################

#### CHANGE: A failed step now halts immediately with a clear [FATAL] message naming the line and exit code.
set -Eeuo pipefail
shopt -s nullglob
trap 'echo "[FATAL] Pipeline aborted at line $LINENO (exit $?). See $PROCESS_README for the full log." >&2' ERR

#### CHANGE: shared helper for bumping a "...vNNN...fa" allele-db filename.
bump_db_version() {
    local fasta="$1" increment="$2" ver new_ver
    ver=$(echo "$fasta" | grep -o 'v[0-9]\{3\}' | cut -c2- | sed 's/^0*//')
    if [ -z "$ver" ]; then
        echo "[ERROR] Could not parse a 'vNNN' version from filename: $fasta" >&2
        return 1
    fi
    new_ver=$(printf '%03d' $((ver + increment)))
    echo "$fasta" | sed "s/v[0-9]\{3\}/v${new_ver}/"
}


DUP_TAGS=""
CUTADAPT=""
CUTADAPT_UNI=""
TMP_RENAME_UPDATED=""

# Processing starts here
now=$(date)
printf "%s\n" "$now"
printf "\n --- Scripts, files, and parameters used ---\n"
printf "# Report: $REPORT\n"
printf "# First sample column: $FIRST_SAMPLE_COL\n"
printf "# Panel design length: $DESIGN_LEN\n"
printf "# Sequencing length: $SEQ_LEN\n"
printf '# BLAST filtering thresholds: %s%% coverage, %s%% identity\n' "$COV" "$IDEN"



printf "\n# 0). Update snpIDs to Chr_00xxxxxxx format in MADC"
printf "\n  # This is necessary because marker IDs do not always follow the Chr_00xxxxxxx format\n"
REPORT_SNPID="${REPORT%????}_snpID.csv"
python3 "$SCRIPTS_DIR/step00_madc_update_snpID_v1.1.py" "$MARKERID_LUT" "$REPORT"



# Remove duplicate tags and keep only unique tags
printf "\n###### Start processing MADC file ######\n"
printf "\n# 1). Remove duplicate tags."
# ! -s: File does not exist or is empty
if [[ ! -s "${DUP_TAGS}" ]]; then
    printf "\n  # No duplicate tags file provided. Skipping duplicate tag removal step.\n"
else
    printf "\n  # Duplicate tags file: $DUP_TAGS\n"
    REPORT_SNPID_RMDUP="${REPORT%????}_snpID_rmDup.csv"
    python3 "$SCRIPTS_DIR/step00_rm_dupTag_closeTag_rawMADC_v1.py" "$DUP_TAGS" "$REPORT_SNPID" "$REPORT_SNPID_RMDUP"
    printf "\n  # MADC file with duplicate tags removed: $REPORT_SNPID_RMDUP\n"
    REPORT_SNPID="$REPORT_SNPID_RMDUP"
fi


printf "\n# 2). Check if there are any duplicate alleles - same microhaplotypes with different IDs"
printf "\n  # If there are duplicates, output the duplicates and stop running the pipeline\n"
# https://stackoverflow.com/questions/947897/how-to-comment-a-block-of-code-in-bash
python3 "$SCRIPTS_DIR/step00_check_allele_uniqueness.py" "$REPORT_SNPID" "$FIRST_SAMPLE_COL"
REPORT_SNPID_DUP="${REPORT%????}_snpID_dup.csv"
if test -f "$REPORT_SNPID_DUP"; then
    printf "\n  # Check duplicated microhaplotypes\n"
    exit 0
else
    printf "\n  # No duplicates found in MADC file\n"
fi


printf "\n# 3). Filter alleles with missing data and generate RefMatch and AltMatch allele fasta and a temporary report\n"
if [ -n "${LUT_INSERT_SIZE_COL}" ]; then
  python3 "$SCRIPTS_DIR/step01_filter_missing_AND_ext_matchAlleles_from_madc_v1.py" "$REPORT_SNPID" --lut "$MARKERID_LUT" --min_samples "$MIN_SAMPLES" --min_reads "$MIN_READS" --split_by_colID "$LUT_INSERT_SIZE_COL"
else
  python3 "$SCRIPTS_DIR/step01_filter_missing_AND_ext_matchAlleles_from_madc_v1.py" "$REPORT_SNPID" --lut "$MARKERID_LUT" --min_samples "$MIN_SAMPLES" --min_reads "$MIN_READS" --split_by_colID ""
fi
# The output of this step is a temporary report file with renamed marker IDs and RefMatch and AltMatch allele fasta files
TMP_RENAME="${REPORT_SNPID%????}_tmp_rename.csv"
MATCH_ALLELES="${REPORT_SNPID%????}_match.fa"


printf "\n#  4). Remove adapters using cutadapt"
#  -n COUNT, --times COUNT
#                        Remove up to COUNT adapters from each read. Default: 1
#  -O MINLENGTH, --overlap MINLENGTH
#                        Require MINLENGTH overlap between read and adapter for an adapter to be found.
#                        Default: 3
# -a ACCGATCTCGTATGCCGTCTTCTGCTTG
if [ "$DESIGN_LEN" -ge "$SEQ_LEN" ]; then
    printf "\n  # Design length is greater than or equal to sequencing length. No need to run cntadapt.\n"
else
    printf "\n  # Panel with mixed insert sizes detected.\n"
    REPORT_DIR="$(dirname "$REPORT")"
    files_to_join=() # Keep all the fasta files to be concatenated
    for file in "$REPORT_DIR"/*_insert*bp.fa; do
        insert=$(echo "$file" | sed -n 's/.*insert\([0-9]*\)bp.*/\1/p') # insert size
        if [ "$insert" -ge "$SEQ_LEN" ]; then
            printf "\n  # Processing: $file"
            printf "\n  # Sequencing length ($SEQ_LEN) is less than or equal to panel insert size ($insert), skipping cutadapt.\n"
            files_to_join+=("$file")
        else
            printf "\n  # Processing: $file"
            printf "\n  # Sequencing length ($SEQ_LEN) is longer than panel insert size ($insert), run cutadapt\n"
            CUTADAPT="${file%.fa}_cutadapt.fa"
            CUT_LOG="${file%.fa}_cutadapt.log"
            "$CUTADAPT_PATH" -a ACCGATCTG -e 0.2 -n 1 --overlap 5 -m "$DESIGN_LEN" -o "$CUTADAPT" "$file" > "$CUT_LOG"

            files_to_join+=("$CUTADAPT")
            # Some Illumina instruments use a two-color chemistry to encode the four bases. This includes the NextSeq and the NovaSeq.
            # In those instruments, a ‘dark cycle’ (with no detected color) encodes a G.
            # However, dark cycles also occur when sequencing “falls off” the end of the fragment.
            # The read then contains a run of high-quality, but incorrect “G” calls at its 3’ end.

            #### CHANGE: `grep -c` exits 1 when the count is zero (even
            #### though it still prints "0") -- under set -e/pipefail that
            #### would wrongly abort the pipeline on a legitimately-empty
            #### fasta. `|| true` keeps the count but ignores that exit code.
            MATCH=$(grep -c ">" "$file" || true)
            MATCH_CUT=$(grep -c ">" "$CUTADAPT" || true)
            printf "\n    # Number of RefMatch and AltMatch extracted from MADC: $MATCH"
            printf "\n    # Number of RefMatch and AltMatch retained after cutadapt: $MATCH_CUT\n"
        fi
    done

    # Concatenate all files
    CUTADAPT="${REPORT_SNPID%.csv}_match_cutadapt.fa"
    if [ "${#files_to_join[@]}" -gt 0 ]; then
      cat "${files_to_join[@]}" > "$CUTADAPT"
      printf "\n  # Final combined microhaplotype fasta file: $CUTADAPT"
    else
      printf "\n  # [ERROR] No files were processed for concatenation.\n"
    fi
fi



printf "\n#  5). Check if there are duplicate alleles after removing adapters AND update allele sequences and read counts\n"
if test -f "$CUTADAPT"; then
    printf "  # CUTADAPT file: $CUTADAPT\n"
    printf "  # Temporary rename file: $TMP_RENAME\n"
    python3 "$SCRIPTS_DIR/step03_check_cutadapt_allele_uniqueness_AND_update_tmp_rename_report_v1.1.py" "$CUTADAPT" "$ALLELE_DB_DIR/$ALLELE_DB" "$TMP_RENAME"
    CUTADAPT_UNI="${CUTADAPT%???}_unique.fa"
    TMP_RENAME_UPDATED="${REPORT_SNPID%????}_tmp_rename_updatedSeq.csv"
fi



printf "\n#  6a). First, create the BLAST database"
# 1. First, create the BLAST database
makeblastdb -in "$ALLELE_DB_DIR/$ALLELE_DB" \
            -dbtype nucl


printf "\n#  6b). BLAST RefMatch and AltMatch against the allele db\n"
# If there are no duplicate alleles, there won't be the a '_match_cutadapt_unique.fa'
# There won't be a '_tmp_rename_updatedSeq.csv' either
# Here, use if else to execute different input files.
if test -f "$CUTADAPT"; then
    if test -f "$CUTADAPT_UNI"; then
      BLAST_DBBLAST="${REPORT_SNPID%????}_match_cutadapt_unique.fa.alleledb.bn"
      blastn -task blastn-short -dust no -soft_masking false -db "$ALLELE_DB_DIR/$ALLELE_DB" -query "$CUTADAPT_UNI" -out "$BLAST_DBBLAST" -evalue 1e-5 -num_threads 6 -max_target_seqs 15 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'
    else
      BLAST_DBBLAST="${REPORT_SNPID%????}_match_cutadapt.fa.alleledb.bn"
      blastn -task blastn-short -dust no -soft_masking false -db "$ALLELE_DB_DIR/$ALLELE_DB" -query "$CUTADAPT" -out "$BLAST_DBBLAST" -evalue 1e-5 -num_threads 6 -max_target_seqs 15 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'
    fi
else
    BLAST_DBBLAST="${REPORT_SNPID%????}_match.fa.alleledb.bn"
    blastn -task blastn-short -dust no -soft_masking false -db "$ALLELE_DB_DIR/$ALLELE_DB" -query "$MATCH_ALLELES" -out "$BLAST_DBBLAST" -evalue 1e-5 -num_threads 6 -max_target_seqs 15 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'
fi


printf "\n#  7). Determine status of RefMatch and AltMatch and Assign fixed IDs to them\n"
if test -f "$CUTADAPT_UNI"; then
    python3 "$SCRIPTS_DIR/step05_parse_madc_matchAlleles_blastn.py" "$ALLELE_DB_DIR/$MATCHCNT_LUT" "$ALLELE_DB_DIR/$ALLELE_DB" "$TMP_RENAME_UPDATED" "$BLAST_DBBLAST" "$SEQ_LEN" --cov_threshold "$COV" --iden_threshold "$IDEN" --allele-type-ground-truth "$ALLELE_TYPE"
    MADC_CLEANED="${REPORT_SNPID%????}_rename_updatedSeq.csv"
else
    python3 "$SCRIPTS_DIR/step05_parse_madc_matchAlleles_blastn.py" "$ALLELE_DB_DIR/$MATCHCNT_LUT" "$ALLELE_DB_DIR/$ALLELE_DB" "$TMP_RENAME" "$BLAST_DBBLAST" "$SEQ_LEN" --cov_threshold "$COV" --iden_threshold "$IDEN" --allele-type-ground-truth "$ALLELE_TYPE"
    MADC_CLEANED="${REPORT_SNPID%????}_rename.csv"
fi



printf "\n#  8). Check if there is a new version DB\n"
# If there are new alleles added to the db, a new version of db will be generated
# Otherwise, no new db
ALLELE_DB_NEW=$(bump_db_version "$ALLELE_DB" 1) || exit 1
if test -f "$ALLELE_DB_DIR/$ALLELE_DB_NEW"; then
    printf "  # New version of db created after adding noval alleles.\n"
    printf "  $ALLELE_DB_DIR/$ALLELE_DB_NEW\n"
    makeblastdb -in "$ALLELE_DB_DIR/$ALLELE_DB_NEW" -dbtype nucl

    printf "\n#  9). Check if there are duplicates."
    printf "  # If there are duplicates, retain only one allele of the duplicated ones\n"
    #cat $ALLELE_DB_NEW_BLAST | awk '$1!=$5 && $10==100 && $11==100.0' | more
    python3 "$SCRIPTS_DIR/step06_check_db_allele_uniqueness_v1.py" "$ALLELE_DB_DIR/$ALLELE_DB_NEW"

    printf "\n#  10). If there are duplicates in the new version DB, update MADC after removing duplicated alleles in db\n"
    DUP="$ALLELE_DB_DIR/$ALLELE_DB_NEW.dup.csv"
    if test -f "$DUP"; then
        ALLELE_DB_NEW_RMDUP=$(bump_db_version "$ALLELE_DB" 2) || exit 1
        makeblastdb -in "$ALLELE_DB_DIR/$ALLELE_DB_NEW_RMDUP" -dbtype nucl
        printf "  # There are duplicate alleles in db. Check if these duplicate alleles are in MADC file.\n"
        python3 "$SCRIPTS_DIR/step06_update_MADC_with_allele_uniqueness_v1.py" "$DUP" "$MADC_CLEANED"
        printf "\n#  11). Add version of the script to the MADC with fixed allele IDs\n"
        if test -f "$CUTADAPT_UNI"; then
            MADC_CLEANED_RMDUP="${REPORT_SNPID%????}_rename_updatedSeq_rmDup.csv"
            MADC_CLEANED_RMDUP_VER="${REPORT_SNPID%????}_rename_updatedSeq_rmDup_${CODE_VER}.csv"
        else
            MADC_CLEANED_RMDUP="${REPORT_SNPID%????}_rename_rmDup.csv"
            MADC_CLEANED_RMDUP_VER="${REPORT_SNPID%????}_rename_rmDup_${CODE_VER}.csv"
        fi
        awk -v val="$CODE_VER" 'NR==1{print "Code_version," $0} NR>1{print val "," $0}' "$MADC_CLEANED_RMDUP" > "$MADC_CLEANED_RMDUP_VER"
        MADC_CLEANED="$MADC_CLEANED_RMDUP"
        ALLELE_DB_NEW="$ALLELE_DB_NEW_RMDUP"
    else
        printf "  # No duplicate alleles found in microhap db\n"
        echo "$ALLELE_DB_DIR/$ALLELE_DB_NEW"
        printf "\n#  11). Add version of the script to the MADC with fixed allele IDs\n"
        if test -f "$CUTADAPT_UNI"; then
            MADC_CLEANED_VER="${REPORT_SNPID%????}_rename_updatedSeq_${CODE_VER}.csv"
        else
            MADC_CLEANED_VER="${REPORT_SNPID%????}_rename_${CODE_VER}.csv"
        fi
        awk -v val="$CODE_VER" 'NR==1{print "Code_version," $0} NR>1{print val "," $0}' "$MADC_CLEANED" > "$MADC_CLEANED_VER"
        printf "  # Version added to as the first column of the output file.\n"
    fi
else
    printf "  # No new alleles found, therefore, no new db generated.\n"
    mv "$PROCESS_README" "$NO_NEW_ALLELE_README"
fi



printf "\n#  12). Tidy up files.\n"
# Extract the directory
DIR_PATH=$(dirname "$REPORT")
HAP_DIR="${DIR_PATH}/matches"
mkdir -p "$HAP_DIR"
match_files=("${DIR_PATH}/"*match.fa*)
if [ "${#match_files[@]}" -gt 0 ]; then
    mv "${match_files[@]}" "$HAP_DIR/"
else
    printf "\n  # No *match.fa* files found to move into %s\n" "$HAP_DIR"
fi
rm -f "$TMP_RENAME" "$TMP_RENAME_UPDATED"
printf "\n###### Complete! #######\n"


####################################
# Generate microhaplotype presence-absence matrix for the MADC file with fixed allele IDs

if [[ -f "$MADC_CLEANED" ]]; then
    printf "\n#  13). Generate microhaplotype presence-absence matrix for the MADC file with fixed allele IDs\n"
    printf "\n  # Input MADC file: $MADC_CLEANED_VER\n"
    printf "\n  # Output haplotype presence-absence matrix will be generated in the same directory as the input MADC file.\n"
    python3 "$HAPPA_SCRIPT_DIR/util_haplotype_presence_absence_matrix.py" \
    --madc "$MADC_CLEANED"\
    --descriptor "$PROJECT_DESCRIPTOR" \
    --read_threshold "$MIN_READS_FOR_PRESENCE"
else
    printf "\n  # [ERROR] MADC file with fixed allele IDs not found. Cannot generate haplotype presence-absence matrix.\n"
fi
HAP_STATUS="${MADC_CLEANED%.csv}_hapStatus.csv"


printf "\n#  14). Create hapSearch, hapPA, and sample_hapStatus directories if they do not exist\n"
HAPSEARCH_DIR="$(dirname "$ALLELE_DB_DIR")/hapSearch"
HAPPA_MATRIX_DIR="$(dirname "$ALLELE_DB_DIR")/hapSearch/hapPA"
HAP_STATUS_DIR="$(dirname "$ALLELE_DB_DIR")/hapSearch/sample_hapStatus"
if [[ -d "${HAPSEARCH_DIR}" ]]; then
  echo "hapSearch directory exists."
  if [[ -d "${HAPPA_MATRIX_DIR}" ]]; then
    echo "HapPA directory exists."
  else
    echo "HapPA directory does not exist. Creating it..."
    mkdir -p "${HAPPA_MATRIX_DIR}"
  fi
  if [[ -d "${HAP_STATUS_DIR}" ]]; then
    echo "sample_hapStatus directory exists."
  else
    echo "sample_hapStatus directory does not exist. Creating it..."
    mkdir -p "${HAP_STATUS_DIR}"
  fi
else
  echo "Directory does not exist. Creating it along with hapPA and sample_hapStatus directories..."
  mkdir -p "${DIR_PATH}"
  mkdir -p "${HAPPA_MATRIX_DIR}"
  mkdir -p "${HAP_STATUS_DIR}"
fi


printf "\n#  15). Copy the haplotype presence-absence matrix to the sample_hapStatus directory\n"
if [[ -f "$HAP_STATUS" ]]; then
    printf "\n  # Copying microhaplotype presence/absence table to $HAP_STATUS_DIR\n"
    cp "$HAP_STATUS" "$HAP_STATUS_DIR/"
else
    printf "\n  # [ERROR] Haplotype presence-absence matrix not generated. Please check the previous steps for errors.\n"
fi

printf "\n#  16). Generate haplotype presence-absence matrix for the MADC file with fixed allele IDs in the hapPA directory\n"


HAPPA_MATRIX_NEW="${dirname, "$ALLELE_DB_DIR"}/hapSearch/hapPA/${ALLELE_DB_NEW%.fa}_hapPA.csv"
if [[ -f $LATEST_HAPPA_MATRIX ]]; then
    printf "\n  # Copying microhaplotype presence/absence table to $HAPPA_MATRIX_DIR\n"
    cp "$LATEST_HAPPA_MATRIX" "$HAPPA_MATRIX_DIR/"
else
    printf "\n  # [ERROR] Haplotype presence-absence matrix not found in hapPA directory. Please check the previous steps for errors.\n"
fi