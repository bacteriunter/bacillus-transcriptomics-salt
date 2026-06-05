#!/usr/bin/env python3

import pandas as pd

ve_tested = pd.read_csv(
    "results/deseq2/velezensis_control_vs_salt/deseq2_summary.tsv",
    sep="\t"
)

pa_tested = pd.read_csv(
    "results/deseq2/paralicheniformis_control_vs_salt/deseq2_summary.tsv",
    sep="\t"
)

ve_deg = pd.read_csv(
    "results/deseq2/velezensis_control_vs_salt/filtered/filtered_deg_summary.tsv",
    sep="\t"
)

pa_deg = pd.read_csv(
    "results/deseq2/paralicheniformis_control_vs_salt/filtered/filtered_deg_summary.tsv",
    sep="\t"
)

table = pd.DataFrame([
    {
        "Species": "Bacillus velezensis",
        "Control_replicates": 3,
        "Salt_replicates": 3,
        "Genes_tested": int(ve_tested["n_genes_tested"].iloc[0]),
        "Upregulated_DEGs": int(ve_deg["upregulated"].iloc[0]),
        "Downregulated_DEGs": int(ve_deg["downregulated"].iloc[0]),
        "Total_DEGs": int(ve_deg["significant_abs_log2FC_1"].iloc[0])
    },
    {
        "Species": "Bacillus paralicheniformis",
        "Control_replicates": 3,
        "Salt_replicates": 3,
        "Genes_tested": int(pa_tested["n_genes_tested"].iloc[0]),
        "Upregulated_DEGs": int(pa_deg["upregulated"].iloc[0]),
        "Downregulated_DEGs": int(pa_deg["downregulated"].iloc[0]),
        "Total_DEGs": int(pa_deg["significant_abs_log2FC_1"].iloc[0])
    }
])

out = "results/tables/Table1_Transcriptome_DEG_Summary.tsv"

table.to_csv(
    out,
    sep="\t",
    index=False
)

print(table)
print(f"\nSaved: {out}")
