# FAQ

## Which thresholds should I use?

Start with the defaults, then test others on your data. [Compare Settings](../tutorials/compare-settings.md) shows
how, with results from a real panel.

## Why did my database grow so much?

A loose read filter lets rare alleles earn permanent IDs. A single sample with an unusual genome can also add many
alleles. The counts at each stage show where the growth came from, and a trial run lets you test a stricter filter
before you accept any new database.

## Why does the tool drop a marker?

A marker drops out when it is missing from the lookup table. The tool warns you and names the marker. Panels use this
to exclude markers, such as ones `build-db` found to be duplicates.

## Can reads be shorter than the alleles?

Yes. Set `database.allele_length` to the length you want. `build-db` extends each allele from the genome, so an 81 bp
read can give a 109 bp allele. This prepares the database for longer reads later.

## Do I need esl-sfetch and seqkit?

Not for the CLI. `build-db` does those two jobs in Python: cutting a range from a FASTA and picking the Ref alleles by
name. On the red clover panel, the results match both tools exactly. The Python version also writes no index file next
to your genome. The bash workflows still use both tools.

## My genome uses NCBI chromosome names. What do I do?

Either rename the genome to match the probe file, or give a
[chromosome map](../getting-started/input-files.md#chromosome-map). Both give the same result. `check` reports a mismatch
before a run starts.

## How do I accept a run's results?

Copy the new database and count table from the run's `versions` folder to where you keep your databases, and use them
for the next project. A run never changes its input database, so a trial has no effect until you do this.

## Can I run projects out of order?

You can, but IDs depend on the order of the runs, since a new allele takes the next free number. Choose one order,
and always start from the latest database.

## Why does the report say `v2`?

`Code_version` records the version of the code that made the report, and this tool is version 2 of the pipeline. Set
`output.code_version` to use another label.

## Can I still use the bash workflows?

Yes. They are in `workflow/`. Set `SCRIPTS_DIR` and `HAPPA_SCRIPT_DIR` to the folders under
`src/hapapp_utils/steps/`.

## Can I run one step on its own?

Yes. Each step still takes the same arguments as its original script, either as a file or as a module:

```bash
python src/hapapp_utils/steps/assign/step05_parse_madc_matchAlleles_blastn.py --help
python -m hapapp_utils.steps.assign.step05_parse_madc_matchAlleles_blastn --help
```

## Why can a run fail with a message about a step?

The tool catches any unexpected failure inside a step and names the step. The process readme holds the full Python
error. If you think you found a bug, please open an issue and attach that file.
