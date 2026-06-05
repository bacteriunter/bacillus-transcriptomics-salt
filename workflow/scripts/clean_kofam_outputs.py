#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

indir = Path("results/enrichment/kofam")

for name in [
    "background_rbh_kofam_accepted",
    "conserved_up_accepted",
    "conserved_down_accepted",
    "opposite_accepted"
]:
    infile = indir / f"{name}.tsv"
    outfile = indir / f"{name}_clean.tsv"

    rows = []

    with open(infile) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 7:
                continue

            rows.append({
                "protein_id": parts[1],
                "KO": parts[2],
                "threshold": parts[3],
                "score": parts[4],
                "evalue": parts[5],
                "definition": " ".join(parts[6:])
            })

    df = pd.DataFrame(rows)
    df.to_csv(outfile, sep="\t", index=False)

    print(name, "rows:", len(df), "proteins:", df["protein_id"].nunique(), "KOs:", df["KO"].nunique())
