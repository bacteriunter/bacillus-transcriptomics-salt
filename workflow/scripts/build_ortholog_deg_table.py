#!/usr/bin/env python3

import pandas as pd
import re

# -------------------------
# GFF parser
# -------------------------

def parse_gff(gff_file):
    rows = []

    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue

            cols = line.rstrip().split("\t")

            if len(cols) < 9:
                continue

            if cols[2] != "CDS":
                continue

            attrs = cols[8]

            m1 = re.search(r"Parent=(gene-[^;]+)", attrs)
            m2 = re.search(r"protein_id=([^;]+)", attrs)

            if m1 and m2:
                rows.append(
                    {
                        "gene_id": m1.group(1),
                        "protein_id": m2.group(1)
                    }
                )

    return pd.DataFrame(rows)

# -------------------------
# mappings
# -------------------------

ve_map = parse_gff(
    "references/gff/GCF_057365115.1/genomic.gff"
)

pa_map = parse_gff(
    "references/gff/GCF_979930925.1/genomic.gff"
)

# -------------------------
# DEGs
# -------------------------

ve_deg = pd.read_csv(
    "results/deseq2/velezensis_control_vs_salt/filtered/DEGs_padj0.05_log2FC1.tsv",
    sep="\t"
)

pa_deg = pd.read_csv(
    "results/deseq2/paralicheniformis_control_vs_salt/filtered/DEGs_padj0.05_log2FC1.tsv",
    sep="\t"
)

ve_deg = ve_deg.merge(ve_map, on="gene_id")
pa_deg = pa_deg.merge(pa_map, on="gene_id")

# -------------------------
# RBH
# -------------------------

rbh = pd.read_csv(
    "results/orthology/mmseqs_rbh/velezensis_vs_paralicheniformis_rbh.tsv",
    sep="\t",
    header=None
)

rbh = rbh.iloc[:, [0,1]]
rbh.columns = [
    "ve_protein",
    "pa_protein"
]

# -------------------------
# merge
# -------------------------

merged = (
    rbh
    .merge(
        ve_deg[[
            "gene_id",
            "protein_id",
            "log2FoldChange"
        ]],
        left_on="ve_protein",
        right_on="protein_id"
    )
    .rename(
        columns={
            "gene_id":"ve_gene",
            "log2FoldChange":"ve_log2FC"
        }
    )
    .drop(columns=["protein_id"])
)

merged = (
    merged
    .merge(
        pa_deg[[
            "gene_id",
            "protein_id",
            "log2FoldChange"
        ]],
        left_on="pa_protein",
        right_on="protein_id"
    )
    .rename(
        columns={
            "gene_id":"pa_gene",
            "log2FoldChange":"pa_log2FC"
        }
    )
    .drop(columns=["protein_id"])
)

def classify(row):

    if row.ve_log2FC > 0 and row.pa_log2FC > 0:
        return "conserved_up"

    if row.ve_log2FC < 0 and row.pa_log2FC < 0:
        return "conserved_down"

    return "opposite"

merged["response_class"] = merged.apply(
    classify,
    axis=1
)

merged.to_csv(
    "results/orthology/conserved_deg_orthologs.tsv",
    sep="\t",
    index=False
)

print()
print(merged["response_class"].value_counts())
print()
print("total shared DEGs =", len(merged))
