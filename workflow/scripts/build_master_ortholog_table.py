#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

outdir = Path("results/tables")
outdir.mkdir(parents=True, exist_ok=True)

orth = pd.read_csv(
    "results/orthology/conserved_deg_orthologs_annotated.tsv",
    sep="\t"
)

kofam = pd.read_csv(
    "results/enrichment/kofam/background_rbh_kofam_accepted_clean.tsv",
    sep="\t",
    header=None,
    names=["ve_protein","KO","threshold","score","evalue","KO_definition"]
)

osm = pd.read_csv(
    "results/orthology/osmoadaptation/osmoadaptation_candidate_summary.tsv",
    sep="\t"
)

# Merge KO annotations using velezensis representative protein
master = orth.merge(
    kofam[["ve_protein","KO","KO_definition"]],
    on="ve_protein",
    how="left"
)

# Functional candidate category assignment
patterns = {
    "glycine_betaine_choline": "opu|gbs|bet|glycine betaine|choline",
    "proline": "proB|proC|putP|proline|pyrroline",
    "sodium_homeostasis": "nha|mnh|Na\\+|cation/H|antiporter",
    "poly_gamma_glutamate": "pgs|poly-gamma-glutamate",
    "oxidative_stress": "oxidoreductase|hydroperoxide|superoxide|redox|peroxidase",
    "transport": "ABC transporter|MFS transporter|MATE|symporter|permease",
    "regulation": "regulator|histidine kinase|response regulator|MarR|GbsR",
    "biofilm_cell_envelope": "biofilm|TasA|matrix|glycosyltransferase|polysaccharide|cell wall|peptidoglycan"
}

def assign_categories(row):
    text = " ".join([
        str(row.get("ve_symbol","")),
        str(row.get("pa_symbol","")),
        str(row.get("ve_product","")),
        str(row.get("pa_product","")),
        str(row.get("KO_definition",""))
    ])

    hits = []
    for cat, pat in patterns.items():
        if pd.Series([text]).str.contains(pat, case=False, regex=True).iloc[0]:
            hits.append(cat)

    return ";".join(hits) if hits else "unclassified"

master["functional_category"] = master.apply(assign_categories, axis=1)

master.to_csv(
    outdir / "master_conserved_ortholog_deg_table.tsv",
    sep="\t",
    index=False
)

summary = (
    master.groupby(["response_class","functional_category"])
    .size()
    .reset_index(name="n")
    .sort_values(["response_class","n"], ascending=[True, False])
)

summary.to_csv(
    outdir / "master_functional_category_summary.tsv",
    sep="\t",
    index=False
)

print("Master rows:", len(master))
print()
print(summary.head(30))
