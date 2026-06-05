#!/usr/bin/env python3

import pandas as pd

df = pd.read_csv(
    "results/tables/Table3_Osmoadaptation_Candidates.tsv",
    sep="\t"
)

core_genes = {
    "proB",
    "proC",
    "putP",
    "betB",
    "gbsB",
    "opuAA",
    "opuAB",
    "nhaC",
    "pgsB",
    "pgsC"
}

def choose_symbol(row):
    ve = str(row["ve_symbol"]).strip()
    pa = str(row["pa_symbol"]).strip()

    if ve and ve != "nan":
        return ve
    if pa and pa != "nan":
        return pa

    return None

df["gene_symbol"] = df.apply(choose_symbol, axis=1)

df = df[df["gene_symbol"].isin(core_genes)].copy()

if "combined_score" not in df.columns:
    df["combined_score"] = (
        df["ve_log2FC"].abs() +
        df["pa_log2FC"].abs()
    )

df = (
    df.sort_values("combined_score", ascending=False)
      .drop_duplicates("gene_symbol")
)

table2 = df[
    [
        "gene_symbol",
        "KO",
        "KO_definition",
        "functional_category",
        "ve_log2FC",
        "pa_log2FC"
    ]
].copy()

table2.columns = [
    "Gene",
    "KO",
    "Function",
    "Category",
    "log2FC_B_velezensis",
    "log2FC_B_paralicheniformis"
]

table2 = table2.sort_values(
    "log2FC_B_velezensis",
    ascending=False
)

outfile = "results/tables/Table2_Core_Osmoadaptation_Genes.tsv"

table2.to_csv(
    outfile,
    sep="\t",
    index=False
)

print(table2)
print("\nSaved:", outfile)
