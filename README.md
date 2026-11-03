For building microhaplotype databases from targeted genotyping platforms

We developed a systematic workflow to process MADC reports and assign standardized microhaplotype IDs. The process begins with updating marker IDs to the standard “chromosomeN_000000000” format from the initial MADC report. As an initial quality control, we filter RefMatch and AltMatch sequences to require a minimum presence of either 10 samples or 5% of the total samples in a project (whichever is smaller), with each sample having at least two reads.

## Command-line tool

`hapapp-utils` runs this workflow from one settings file. Each step remains its own Python module.

```bash
pip install -e ".[test]"             # also needs BLAST+ and MMseqs2 on your PATH
hapapp-utils init run.yaml           # write a commented settings file
hapapp-utils check -c run.yaml       # check tools, inputs and settings
hapapp-utils build-db -c build.yaml  # a panel's first database, from the reference genome
hapapp-utils assign -c run.yaml      # fixed allele IDs and a new database version
```

See [docs/cli.md](docs/cli.md) for the settings, what a run writes, and how to work on the steps.
