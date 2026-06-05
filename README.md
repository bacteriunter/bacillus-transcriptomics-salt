# Bacillus transcriptomic responses to salt stress

This repository contains a reproducible workflow for analyzing public RNA-seq datasets from Bacillus species exposed to salt stress.

## Research question

Do different Bacillus species share a common transcriptomic response to salt stress, or do they deploy lineage-specific regulatory strategies?

## Project structure

- `metadata/`: dataset screening and sample annotation files.
- `raw_data/`: downloaded FASTQ files, reference genomes, and annotations.
- `scripts/`: standalone analysis scripts.
- `workflow/`: Snakemake rules, scripts, and environments.
- `results/`: quality control, alignments, counts, differential expression, and enrichment results.
- `figures/`: exploratory and final figures.
- `docs/`: manuscript drafts, methods notes, and documentation.
- `logs/`: execution logs.

## Reproducibility

The workflow will be managed with Snakemake. Raw sequencing data will not be stored in GitHub; only accession tables, scripts, configuration files, and processed summary outputs will be versioned.
