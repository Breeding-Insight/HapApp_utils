# Output Files

The output folder must be new or empty. Both workflows write these files at its top level:

| File | Contents |
|------|----------|
| `resolved_config.yaml` | The settings the run used, after defaults and overrides |
| `run_manifest.json` | Status, counts, warnings, tool versions, input checksums and output paths |

Both also write a process readme in `versions/`. It opens with the input files, settings and tool versions,
and continues with everything each step printed.

## assign

### Main results

| File | Contents |
|------|----------|
| `report/*_rename_v2.csv` | The report with fixed allele IDs and a `Code_version` column |
| `report/*_rename.csv` | The same report without that column |
| `versions/*_v00N.fa` | The new allele database, with its BLAST index files |
| `versions/*_v00N_matchCnt_lut.txt` | The new count table |
| `versions/*_v00N_process.readme` | The run's log, named for the new version |
| `versions/*_noNewAllele.readme` | The run's log instead, if it found no new alleles, named for the database it kept |
| `hapSearch/hapPA/<database>_<project>_hapPA.csv` | The presence/absence matrix, merged if you gave an earlier one |
| `hapSearch/sample_hapStatus/*_hapStatus.csv` | One column per sample, with 1 for present and 0 for absent |

A report that adapter trimming changed has `updatedSeq` in its name, as in `*_rename_updatedSeq_v2.csv`.

### Supporting files

| File | Contents |
|------|----------|
| `report/*_snpID.csv` | The MADC with your standard marker IDs |
| `report/*_snpID_dup.csv` | Allele IDs that share a sequence, if the run stopped for that reason |
| `report/*_match.log` | The filter's settings and counts |
| `report/matches/` | The FASTA of kept alleles, the BLAST results and each allele's best hit |
| `versions/*_mismatchedType_alleles.csv` | New alleles whose closest match has the other Ref or Alt type |
| `versions/*_exact_match_type_corrections.csv` | Database alleles relabeled to match the report |
| `lut/` | The copy of the lookup table the run used |

The new database version sits next to a copy of the input version. If the new database holds duplicate sequences,
the run writes one more version without them, and the summary shows which version to use.

## build-db

### Main results

| File | Contents |
|------|----------|
| `versions/<name>_v001.fa` | Ref and Alt alleles in uppercase, one line per sequence, with BLAST index files |
| `versions/<name>_v001_matchCnt_lut.txt` | The count table, starting at zero |
| `versions/<name>_v001_process.readme` | The run's log |
| `probe/*_rmDupTag_snpID_lut.csv` | The lookup table without removed markers |
| `probe/*_snpID_lut.csv` | The lookup table before removal, which is also the one to use if the run removed none |
| `probe/*.botloci` | Markers on the minus strand of the genome |
| `report/*_markers_removed.tsv` | Removed markers, if any |

In `*_markers_removed.tsv`, the second column names the removed marker. The first names the marker it duplicates,
or says `no_alleles`.

### Supporting files

| File | Contents |
|------|----------|
| `probe/*.csv` | The probe file as plain CSV |
| `probe/*_dup.csv` | Markers that share a position in the genome |
| `probe/genome_lengths.txt` | Chromosome lengths read from the genome |
| `probe/*_sfetchKeys.txt`, `.fa` | Coordinates and sequences of the flanks |
| `probe/*_ref_alt.fa`, `*_ref_alt_rev.fa` | Ref and Alt flanks, before and after strand correction |
| `report/*_snpID.csv` | The MADC with your standard marker IDs |
| `report/*_ref_alt_amplicons.fa` | The MADC's Ref and Alt amplicons |
| `report/*_iupac.csv` | Ref alleles in the MADC that contain IUPAC codes |
| `report/*.bn` | BLAST results |
| `report/*_sfetch.fa` | All alleles before duplicate removal |
| `report/*_dupTag_cluster.tsv` | The MMseqs2 clusters |

## run_manifest.json

This file records each run. Its `counts` section holds the numbers the run printed:

```json
{
  "hapapp_utils_version": "2.0.0",
  "status": "finished",
  "first_sample_column": 17,
  "samples": 94,
  "tools": {"makeblastdb": "makeblastdb: 2.17.0+", "blastn": "blastn: 2.17.0+"},
  "inputs": {"madc": {"path": "...", "sha256": "..."}},
  "outputs": {"report": "report/...", "allele_db": "versions/..."},
  "counts": {
    "alleles_tested": 77491,
    "filter_samples_required": 2,
    "alleles_kept_by_filter": 9438,
    "alleles_existing": 708,
    "alleles_new": 8129,
    "alleles_removed_after_blast": 601,
    "database_before": 16400,
    "database_after": 24529
  },
  "warnings": ["1 MADC marker(s) are not in the lookup table and will be dropped: ['Chr06_021528607']"],
  "process_readme": "versions/species_allele_db_v004_process.readme"
}
```

`status` reads `finished`, `failed`, `interrupted` or `stopped after <stage>`.
