# Validation

We checked `hapapp-utils` against the bash workflows and an existing database, using a DArTag panel with
3,000 markers and 94 samples.

## build-db

| Check | Result |
|-------|--------|
| Against the panel's existing database (v003) | All 5,998 Ref and Alt alleles match, ignoring letter case, with the same 2,999 markers |
| Against the existing lookup table | Same markers, with matching IDs, positions and alleles |
| Duplicate detection | Finds `Chr06_021528607`, the one marker the panel had already excluded |
| Against the `build01` bash workflow | Same alleles |
| Python fetching against `esl-sfetch` | 3,000 flanks and 6,000 alleles match in name, order, sequence and letter case |
| Python filtering against `seqkit` | The same 3,000 Ref alleles |
| Reads cut to 81 bp, with 109 bp alleles | All 4,128 alleles match the `esl-sfetch` result and the full-length reads |
| NCBI chromosome names with a map, against a renamed genome | Byte-identical database and count table |
| Two runs in a row | Byte-identical database and count table |

v003 holds 1,294 alleles with lowercase bases from soft masking in the genome, and `build-db` writes uppercase.
BLAST ignores case here.

In the 81 bp test, 936 markers dropped out because cutting the reads removed their SNP. A real panel would not lose
these markers, since its reads cover the variant.

## assign

| Check | Result |
|-------|--------|
| Against the `build02` bash workflow | 11 output files match byte for byte, using the same code version and the same step fixes |
| A second threshold combination (1 sample, 30 reads) | Matches an earlier run by hand |

The comparison uses the fixed steps. An unfixed `build02` from before version 2 gives different results, because
its filter counted summary columns as samples.

## Test suite

The test suite runs each workflow on a small synthetic dataset. It covers cases that panel lacks: indels,
markers the tool cannot place, mixed insert sizes, adapter trimming and reports without summary columns.

To compare a real run with saved results:

```bash
HAPAPP_BASELINE_DIR=path/to/saved/run HAPAPP_BASELINE_CONFIG=path/to/run.yaml pytest -k real_data
```

## Known limits

- The real-data checks cover one SNP panel. Indels have only run on synthetic data.
- BLAST+, MMseqs2 and cutadapt can change results between versions. Each run records the versions it used in
  `run_manifest.json`.
