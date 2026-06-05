#!/usr/bin/env python3

import pandas as pd

cols = [
    "ve_protein",
    "pa_protein",
    "identity",
    "alignment_length",
    "mismatches",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore"
]

df = pd.read_csv(
    "results/orthology/mmseqs_rbh/velezensis_vs_paralicheniformis_rbh.tsv",
    sep="\t",
    header=None,
    names=cols
)

out = "results/tables/TS3_All_RBH_Orthologs.tsv"

df.to_csv(out, sep="\t", index=False)

print("Rows:", len(df))
print("Saved:", out)
