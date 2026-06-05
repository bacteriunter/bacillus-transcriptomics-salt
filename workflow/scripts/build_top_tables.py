#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

input_file = "results/tables/master_conserved_ortholog_deg_table.tsv"

outdir = Path("results/tables")
outdir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(input_file, sep="\t")

# One row per ortholog representative
df = df.drop_duplicates(subset=["ve_protein"], keep="first").copy()

df["score"] = df["ve_log2FC"].abs() + df["pa_log2FC"].abs()

# Table 1
up = (
    df[df["response_class"] == "conserved_up"]
    .sort_values("score", ascending=False)
    .head(50)
)

up.to_csv(
    "results/tables/Table1_Top50_ConservedUp.tsv",
    sep="\t",
    index=False
)

# Table 2
down = (
    df[df["response_class"] == "conserved_down"]
    .sort_values("score", ascending=False)
    .head(50)
)

down.to_csv(
    "results/tables/Table2_Top50_ConservedDown.tsv",
    sep="\t",
    index=False
)

# Table 3
osmocats = [
    "proline",
    "glycine_betaine_choline",
    "sodium_homeostasis",
    "poly_gamma_glutamate"
]

osmotic = df[
    df["functional_category"]
    .fillna("")
    .apply(lambda x: any(cat in x for cat in osmocats))
].copy()

osmotic = osmotic.sort_values("score", ascending=False)

osmotic.to_csv(
    "results/tables/Table3_Osmoadaptation_Candidates.tsv",
    sep="\t",
    index=False
)

print()
print("Top conserved_up:", len(up))
print("Top conserved_down:", len(down))
print("Osmoadaptation candidates:", len(osmotic))
print("Unique orthologs in source:", df["ve_protein"].nunique())
print()
