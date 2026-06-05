#!/usr/bin/env python3

import pandas as pd
import re
from pathlib import Path

infile = "results/orthology/conserved_deg_orthologs_annotated.tsv"
outdir = Path("results/orthology/osmoadaptation")
outdir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(infile, sep="\t")

patterns = {
    "glycine_betaine_choline": r"opu|gbs|bet|glycine betaine|choline",
    "proline": r"proB|proC|putP|proline|pyrroline",
    "sodium_homeostasis": r"nha|mnh|Na\+|cation/H|antiporter",
    "poly_gamma_glutamate": r"pgs|poly-gamma-glutamate",
    "oxidative_stress": r"oxidoreductase|hydroperoxide|superoxide|redox|peroxidase",
    "transport": r"ABC transporter|MFS transporter|MATE|symporter|permease",
    "regulation": r"regulator|histidine kinase|response regulator|MarR|GbsR",
    "biofilm_cell_envelope": r"biofilm|TasA|matrix|glycosyltransferase|polysaccharide|cell wall|peptidoglycan"
}

rows = []

for category, pattern in patterns.items():
    mask = (
        df["ve_symbol"].fillna("").str.contains(pattern, case=False, regex=True) |
        df["pa_symbol"].fillna("").str.contains(pattern, case=False, regex=True) |
        df["ve_product"].fillna("").str.contains(pattern, case=False, regex=True) |
        df["pa_product"].fillna("").str.contains(pattern, case=False, regex=True)
    )

    sub = df[mask].copy()
    sub["candidate_category"] = category

    sub.to_csv(
        outdir / f"{category}_candidates.tsv",
        sep="\t",
        index=False
    )

    counts = (
        sub["response_class"]
        .value_counts()
        .to_dict()
    )

    rows.append({
        "candidate_category": category,
        "total": len(sub),
        "conserved_up": counts.get("conserved_up", 0),
        "conserved_down": counts.get("conserved_down", 0),
        "opposite": counts.get("opposite", 0)
    })

summary = pd.DataFrame(rows)
summary.to_csv(
    outdir / "osmoadaptation_candidate_summary.tsv",
    sep="\t",
    index=False
)

print(summary)
