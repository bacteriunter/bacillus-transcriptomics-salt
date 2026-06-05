#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
import subprocess

samples = pd.read_csv("metadata/sample_annotation/samples.tsv", sep="\t")

Path("results/trimmed_reads").mkdir(parents=True, exist_ok=True)
Path("results/qc_reports/fastp").mkdir(parents=True, exist_ok=True)
Path("logs/fastp").mkdir(parents=True, exist_ok=True)

for _, row in samples.iterrows():
    srr = row["run_accession"]
    layout = str(row["library_layout"]).upper()

    print(f"Processing {srr} ({layout})")

    if layout == "PAIRED":
        cmd = [
            "fastp",
            "-i", f"raw_data/fastq/{srr}_1.fastq.gz",
            "-I", f"raw_data/fastq/{srr}_2.fastq.gz",
            "-o", f"results/trimmed_reads/{srr}_1.trimmed.fastq.gz",
            "-O", f"results/trimmed_reads/{srr}_2.trimmed.fastq.gz",
            "-h", f"results/qc_reports/fastp/{srr}.html",
            "-j", f"results/qc_reports/fastp/{srr}.json",
            "-w", "2"
        ]
    else:
        cmd = [
            "fastp",
            "-i", f"raw_data/fastq/{srr}.fastq.gz",
            "-o", f"results/trimmed_reads/{srr}.trimmed.fastq.gz",
            "-h", f"results/qc_reports/fastp/{srr}.html",
            "-j", f"results/qc_reports/fastp/{srr}.json",
            "-w", "2"
        ]

    with open(f"logs/fastp/{srr}.log", "w") as log:
        subprocess.run(cmd, stdout=log, stderr=log, check=True)
