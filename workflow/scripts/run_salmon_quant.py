#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
import subprocess

samples = pd.read_csv("metadata/sample_annotation/samples.tsv", sep="\t")
refs = pd.read_csv("metadata/reference_genomes/sample_reference_map.tsv", sep="\t")

samples = samples.merge(
    refs[["species", "selected_accession"]],
    on="species",
    how="left"
)

Path("results/salmon_quant").mkdir(parents=True, exist_ok=True)
Path("logs/salmon_quant").mkdir(parents=True, exist_ok=True)

for _, row in samples.iterrows():
    srr = row["run_accession"]
    layout = str(row["library_layout"]).upper()
    acc = row["selected_accession"]

    outdir = Path(f"results/salmon_quant/{srr}")
    if (outdir / "quant.sf").exists():
        print(f"Skipping {srr}, already done")
        continue

    print(f"Quantifying {srr} ({layout}) against {acc}")

    index = f"results/salmon_index/{acc}"

    if layout == "PAIRED":
        cmd = [
            "salmon", "quant",
            "-i", index,
            "-l", "A",
            "-1", f"raw_data/fastq/{srr}_1.fastq.gz",
            "-2", f"raw_data/fastq/{srr}_2.fastq.gz",
            "-p", "2",
            "--validateMappings",
            "-o", str(outdir)
        ]
    else:
        cmd = [
            "salmon", "quant",
            "-i", index,
            "-l", "A",
            "-r", f"raw_data/fastq/{srr}.fastq.gz",
            "-p", "2",
            "--validateMappings",
            "-o", str(outdir)
        ]

    with open(f"logs/salmon_quant/{srr}.log", "w") as log:
        subprocess.run(cmd, stdout=log, stderr=log, check=True)

print("Done.")
