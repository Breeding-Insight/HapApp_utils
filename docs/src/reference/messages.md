# Warnings and Errors

## How the tool reports things

While a run goes on, each stage prints a line with its counts. Anything else you should know appears as a
`WARNING` line, and the run carries on. Warnings also go into the process readme and `run_manifest.json`.

An error stops the run and prints `[ERROR]`. If the tool caught the error before starting, it creates no output
folder. If the error came later, the output folder holds a partial run and a process readme with the details.

## Errors before a run starts

`check` reports the same problems without starting a run.

| Message | Cause | What to do |
|---------|-------|------------|
| `Settings file not found` | Wrong path to the settings file | Check the path |
| `Missing required setting(s)` | The file lacks a required setting | Add the settings it names |
| `Unknown setting(s) ... (did you mean ...?)` | A misspelled setting, or a setting that belongs to another command | Use the name it suggests, or run the `init` command it shows to see the valid settings |
| `holds settings for assign, not build-db` | The settings file belongs to the other command | Run the other command, or use `init --workflow` |
| `was read as False. Put its value in quotes.` | YAML treats words such as `NO` as true or false | Put the value in quotes |
| `is not installed or not on PATH` | A program is missing | Follow [Installation](../getting-started/installation.md) |
| `not found` or `is a folder, not a file` | Wrong path to an input | Check the path in the settings |
| `is not a text file` | The input is Excel, gzipped or another binary format | Export it as CSV, or decompress it |
| `no header row starting with AlleleID` | The file is not a MADC | Check `inputs.madc` |
| `are summary columns, not samples` | `panel.first_sample_column` points into the summary columns | Remove the setting. The tool finds the column itself |
| `None of the MADC's markers are in the lookup table` | The first column does not hold the MADC's marker names | Use the lookup table for this panel |
| `needs at least two comma-separated columns` | The lookup table has a blank or short line | Fix the line it names |
| `have no Ref_0001 and Alt_0002 allele` | The database lacks markers this run uses | Use this panel's database, or list the markers in `inputs.dup_tags` |
| `Could not read inputs.match_counts` | The count table has an unexpected format | Check the file |
| `Output folder is not empty` | The folder holds an earlier run | Choose a new folder with `-o` |
| `Cannot write to the output folder` | You lack permission there | Choose another folder |
| `does not look like a DArT probe design file` | The file lacks the expected columns | Check `inputs.probe_design`. A marker selection sheet will not work |
| `could not be read as an Excel file` | The file is corrupt or not an Excel file | Save it again as `.xlsx` or `.csv` |
| `Chromosome(s) not in the genome` | The probe file and genome use different names | Rename the genome, or set `inputs.chromosome_map` |

## Errors during a run

| Message | Cause | What to do |
|---------|-------|------------|
| `[STOPPED] The MADC has the same sequence under different allele IDs` | Two rows share a sequence, and the run exits with code 3 | Review the `_dup.csv` file it names |
| `No allele passed the read filter` | The filter removed every allele | Lower `filter.min_samples` or `filter.min_reads` |
| `Trimming left no alleles` | Every allele came out shorter than the design length | Check `panel.design_length` and `panel.adapter` |
| `sample columns but the MADC has` | The report's sample count looks wrong | Check `panel.first_sample_column` |
| `new allele ID(s) were issued with no sequence` | The new database would disagree with its count table | Report this as a bug, with the process readme |
| `No MADC amplicon aligns to its flanking sequence` | The genome does not fit the panel | Check the assembly and the chromosome names |
| `No marker's Ref or Alt base matches the genome` | The genome does not fit the panel | Check the assembly and the chromosome names |
| `<step> stopped with <error>` | A step failed unexpectedly | Read the process readme, which holds the full error |
| `Something unexpected went wrong` | The tool failed outside any step | Run again with `hapapp-utils --debug` to see the full error |

## Warnings

| Warning | Meaning | What to do |
|---------|---------|------------|
| `MADC marker(s) are not in the lookup table` | The run drops those markers | Nothing, if you excluded them on purpose |
| `the first N sample column(s) will be left out` | `panel.first_sample_column` is past the first sample | Remove the setting unless you meant this |
| `panel.sequencing_length is X, but the MADC's Ref alleles are Y bp` | The setting and the file disagree | Fix the setting |
| `Repeated sample names will get a numeric suffix` | Two samples share a name | Update your passport data to match |
| `inputs.dup_tags is empty` or `nothing removes them` | The tags file removes fewer markers than you may expect | Check that it uses your standard IDs |
| `has no version tag` | The database file name lacks a tag such as `_v003` | Rename the file |
| `allele ID(s) in the count table have no sequence` | An earlier run left the database and count table out of step | See below |
| `alleles were dropped for being shorter than` | Trimming shortened those alleles below the design length | Check the design length if the number looks high |
| `No kept allele was shortened by trimming` | No allele held the adapter | Check `panel.design_length` and `panel.adapter` |
| `No allele had a BLAST hit at its own marker` | The database does not fit this MADC | Check `inputs.allele_db` |
| `duplicate sequence(s) in the new database were collapsed` | The run wrote one more database version | Use the version the summary names |
| `flank_length ... is shorter than ... allele_length` | Some alleles may come out short | Raise `database.flank_length` |
| `allele_length ... is shorter than the MADC's` | The alleles will be shorter than the reads | Raise `database.allele_length` |
| `MADC marker(s) are not in the probe design file` | The run leaves those markers out | Check that the files belong to the same panel |
| `do not match the genome at their positions` | The run removes those markers | Check them if the number looks high |
| `SNP alleles are shorter than requested` | The flanks ran out | Raise `database.flank_length` |

### Count table and database out of step

Earlier versions of the pipeline could lose `Other` alleles when trimming found duplicates. The count table then
counted alleles that the database did not hold. A run on such a database still works, and issues new IDs above the
old numbers, which leaves gaps. The orders that ran with that version may also lack those rows in their reports,
so consider processing them again.

## Warnings from inside a step

The tool also shows any line a step prints with `[WARNING]`, up to five for each stage. The process readme holds
the rest.
