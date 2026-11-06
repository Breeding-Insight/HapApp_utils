# Input Files

## MADC report

Used by `build-db` and `assign`. `hapapp-utils` accepts two layouts and finds the first sample column
from the header row:

- The raw DArT report, with rows of `*` above the header and 13 summary columns such as `CallRate`
  before the samples.
- A shorter report with only `AlleleID`, `CloneID`, `AlleleSequence` and the samples.

Each marker needs a `Ref` and an `Alt` row. Export the report as CSV, not Excel.

The MADC report can come straight from DArT genotyping, or you can generate one from FASTQ files with
[MADCATS](https://madcats.readthedocs.io). [What is a MADC file?](https://softwarehub.breedinginsight.org/learn/what-is-a-madc-file.html)
explains the format.

```text
AlleleID,CloneID,AlleleSequence,5236,5237,...
Chr01_000120735|Ref,Chr01_000120735,GGTGATGAGTGGAGTT...,0,12,...
Chr01_000120735|Alt,Chr01_000120735,GGTGATGAGTGGAGTT...,0,3,...
Chr01_000120735|RefMatch,Chr01_000120735,GGTGATGAGTGGAA...,0,0,...
```

## Probe design file

Used by `build-db`. This is the template you submit to DArT, as `.xlsx`, `.csv` or tab-delimited
`.txt`. The tool skips any notes above the row that starts with `Marker Name`, and reads these columns
by position:

| Column | Name | Example |
|--------|------|---------|
| 1 | Marker Name | `Chr01_000120735` |
| 4 | Chrom | `Chr01` |
| 5 | Chrom Pos | `120735` |
| 6 | VariantAllelesDef | `[T/A]` |
| 7 | Marker_Type | `SNP` or `Indel` |
| 8 | Essential Marker | `2` |
| 9 | Comments | `Marker` |

The marker's standard ID comes from columns 4 and 5, as `Chr01_000120735`.

## Reference genome

Used by `build-db`. A FASTA file, not gzipped. The chromosome names must match the `Chrom` column of the
probe file. If they differ, you have two choices:

- Use a copy of the genome renamed to match.
- Give a [chromosome map](#chromosome-map), and keep the original file.

Both give the same alleles. `hapapp-utils` never writes next to the genome.

## Chromosome map

Optional, for `build-db`. A CSV with no header, listing the probe file's name and then the genome's
name:

```text
Chr01,NC_012345.1
Chr02,NC_012346.1
```

## Marker lookup table

Used by `assign`, and written by `build-db`. A CSV whose first column holds the MADC's marker names and
whose second column holds your standard IDs. A header row is optional, and must start with
`Panel_markerID`. The tool ignores any other columns, unless you set
[`panel.insert_size_column`](../reference/settings.md#panel).

```text
Panel_markerID,Marker_ID,Chr,Pos,Ref,Alt,Type,Indel_pos,Priority,Note
Chr01_000120735,Chr01_000120735,Chr01,120735,T,A,SNP,,2,Marker
```

The tool drops any MADC marker missing from this table and warns you. Use this to exclude markers.

## Allele database

Used by `assign`, and written by `build-db` and `assign`. A FASTA file. Every marker the run uses needs
`|Ref_0001` and `|Alt_0002` alleles.

```text
>Chr01_000120735|Ref_0001
GGTGATGAGTGGAGTTGATGATTTGTCTAATTCC...
>Chr01_000120735|RefMatch_0001
GGTGATGAGTGGAGTTGATGATTTGTCTAATTCA...
```

Keep the version tag in the file name, such as `_v003`. The tool raises it to name the next version.

## Count table

Used by `assign`, and written by `build-db` and `assign`. A tab-delimited file with the number of IDs
issued for each allele class at each marker. The third column, when present, holds the count from the
previous version.

```text
Chr01_000120735|RefMatch	1	0
Chr01_000120735|AltMatch	0	0
```

## Duplicate tags file

Optional, for `assign`. A text file with one marker ID per line, using your standard IDs. The run removes
these markers before filtering.
