# How Allele IDs Work

## ID format

An allele ID has a marker, a type and a serial number:

```text
Chr01_000120735|RefMatch_0002
```

| Type | Meaning |
|------|---------|
| `Ref_0001` | The reference allele. Each marker has one |
| `Alt_0002` | The alternate allele. Each marker has one |
| `RefMatch_NNNN` | Another allele that carries the reference base at the target variant |
| `AltMatch_NNNN` | Another allele that carries the alternate base at the target variant |
| `Other_NNNN` | Another allele that carries neither |

The MADC reports each allele's type. The count table holds the last serial number issued for each type at each
marker, so the tool never issues a number twice.

## The read filter

The filter runs first, and decides which alleles go on to BLAST. It keeps an allele if enough samples each have
enough reads of it:

```text
samples with at least filter.min_reads reads >= min(filter.min_samples, 5% of samples, rounded up)
```

With the defaults and 94 samples, an allele needs at least 5 samples with at least 2 reads each. The tool
counts only sample columns, never the MADC's summary columns. Ref and Alt alleles always pass.

## Adapter trimming

If `panel.sequencing_length` exceeds `panel.design_length`, reads can run into the adapter. Cutadapt removes it,
and drops any allele shorter than the design length afterward. A warning tells you how many it dropped.

Two alleles that differed only past the adapter become identical at this point. The tool merges them, keeps the
first one's ID, and adds their read counts together.

For panels with mixed insert sizes, set `panel.insert_size_column`. The tool then trims only the markers whose
insert is shorter than the read.

## The BLAST decision

The tool aligns each kept allele to the database, and looks only at hits to the allele's own marker. It takes the
hit with the lowest E-value, then applies these rules in order:

| Best hit | Result |
|----------|--------|
| Identical, full length, and the same type | Reuse that ID |
| Identical to `Ref_0001` or `Alt_0002` | Remove the allele, since the database already has that sequence as Ref or Alt |
| Identical, but a different type | Depends on `blast.allele_type_truth` |
| Reaches `min_coverage` and `min_identity` | Issue a new ID |
| Anything else | Remove the allele |

Coverage here means the aligned length as a percent of the allele's length. BLAST stops an alignment at an indel
or a cluster of mismatches near the end, so a real allele can fall short of a high threshold. See
[Compare Settings](../tutorials/compare-settings.md).

An allele with no hit at its own marker also goes unassigned. The `blast` stage reports how many alleles hit their
own marker.

## When the MADC and the database disagree

Sometimes an allele's type in the MADC differs from the type of its match in the database.
`blast.allele_type_truth` decides what happens.

| Situation | `report` (default) | `db` |
|-----------|--------------------|------|
| Exact match with a different type | Issues an ID of the MADC's type, and retires the old ID | Keeps the database's ID |
| Near match with a different type | Issues a new ID of the MADC's type | Issues a new ID of the database's type |

The tool lists near matches with a different type in `versions/*_mismatchedType_alleles.csv`, and retired IDs in
`versions/*_exact_match_type_corrections.csv`.

## Consistency checks

The run stops if any of these fail:

- Two allele IDs in the MADC share a sequence.
- The report's sample count differs from the MADC's.
- A new ID has no sequence in the new database.
- Reordering the report during trimming would lose or invent a row.

## Duplicates in the new database

If the new database holds the same sequence under two IDs, the tool writes one more version with one copy, and
adds the read counts of the removed allele to the one it kept.
