#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
from Bio import SeqIO

outdir = Path("results/enrichment/input_sets")
outdir.mkdir(parents=True, exist_ok=True)

rbh = pd.read_csv(
    "results/orthology/mmseqs_rbh/velezensis_vs_paralicheniformis_rbh.tsv",
    sep="\t",
    header=None
)

rbh = rbh.iloc[:, [0,1]]
rbh.columns = ["ve_protein", "pa_protein"]

annot = pd.read_csv(
    "results/orthology/conserved_deg_orthologs_annotated.tsv",
    sep="\t"
)

# Background: all RBH velezensis proteins
background = pd.DataFrame({
    "protein_id": sorted(rbh["ve_protein"].unique())
})

background.to_csv(
    outdir / "background_rbh_velezensis_proteins.txt",
    index=False,
    header=False
)

# Class-specific sets
for cls in ["conserved_up", "conserved_down", "opposite"]:
    sub = annot[annot["response_class"] == cls]
    ids = pd.DataFrame({
        "protein_id": sorted(sub["ve_protein"].unique())
    })

    ids.to_csv(
        outdir / f"{cls}_velezensis_proteins.txt",
        index=False,
        header=False
    )

# Write FASTA files
proteins = SeqIO.to_dict(
    SeqIO.parse("results/orthology/mmseqs_rbh/velezensis.faa", "fasta")
)

for txt in outdir.glob("*_proteins.txt"):
    ids = [x.strip() for x in txt.read_text().splitlines() if x.strip()]
    fasta_out = txt.with_suffix(".faa")

    with open(fasta_out, "w") as handle:
        for pid in ids:
            if pid in proteins:
                SeqIO.write(proteins[pid], handle, "fasta")

print("Prepared enrichment input sets:")
for f in sorted(outdir.glob("*")):
    print(f)
