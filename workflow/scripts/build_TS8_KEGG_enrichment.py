#!/usr/bin/env python3

import pandas as pd

up = pd.read_csv(
    "results/enrichment/kegg_ko/conserved_up_kegg_ko_enrichment.tsv",
    sep="\t"
)

down = pd.read_csv(
    "results/enrichment/kegg_ko/conserved_down_kegg_ko_enrichment.tsv",
    sep="\t"
)

opp = pd.read_csv(
    "results/enrichment/kegg_ko/opposite_kegg_ko_enrichment.tsv",
    sep="\t"
)

combined = pd.concat(
    [up, down, opp],
    ignore_index=True
)

combined = combined.sort_values(
    ["set", "padj", "pvalue"],
    ascending=[True, True, True]
)

outfile = "results/tables/TS8_KEGG_KO_Enrichment.tsv"

combined.to_csv(
    outfile,
    sep="\t",
    index=False
)

print("\nRows:", len(combined))
print("Saved:", outfile)

print("\nCounts by set:")
print(combined["set"].value_counts())
