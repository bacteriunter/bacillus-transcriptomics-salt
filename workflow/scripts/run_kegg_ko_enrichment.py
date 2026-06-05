#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

outdir = Path("results/enrichment/kegg_ko")
outdir.mkdir(parents=True, exist_ok=True)

def read_kofam(path):
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["protein_id","KO","threshold","score","evalue","definition"]
    )
    df = df.drop_duplicates(subset=["protein_id","KO"])
    return df

background = read_kofam(
    "results/enrichment/kofam/background_rbh_kofam_accepted_clean.tsv"
)

sets = {
    "conserved_up": read_kofam(
        "results/enrichment/kofam/conserved_up_accepted_clean.tsv"
    ),
    "conserved_down": read_kofam(
        "results/enrichment/kofam/conserved_down_accepted_clean.tsv"
    ),
    "opposite": read_kofam(
        "results/enrichment/kofam/opposite_accepted_clean.tsv"
    )
}

bg_proteins = set(background["protein_id"])
bg_kos = set(background["KO"])

for set_name, fg in sets.items():

    fg_proteins = set(fg["protein_id"])
    fg_kos = set(fg["KO"])

    rows = []

    for ko in sorted(bg_kos):

        bg_with_ko = set(background.loc[background["KO"] == ko, "protein_id"])
        fg_with_ko = set(fg.loc[fg["KO"] == ko, "protein_id"])

        a = len(fg_with_ko)
        b = len(fg_proteins) - a
        c = len(bg_with_ko - fg_with_ko)
        d = len(bg_proteins) - len(fg_proteins) - c

        if a == 0:
            continue

        oddsratio, pvalue = fisher_exact(
            [[a, b], [c, d]],
            alternative="greater"
        )

        definition = (
            background.loc[background["KO"] == ko, "definition"]
            .dropna()
            .iloc[0]
        )

        rows.append({
            "set": set_name,
            "KO": ko,
            "definition": definition,
            "foreground_with_KO": a,
            "foreground_total": len(fg_proteins),
            "background_with_KO": len(bg_with_ko),
            "background_total": len(bg_proteins),
            "oddsratio": oddsratio,
            "pvalue": pvalue
        })

    res = pd.DataFrame(rows)

    if len(res) > 0:
        res["padj"] = multipletests(
            res["pvalue"],
            method="fdr_bh"
        )[1]

        res = res.sort_values(
            ["padj", "pvalue", "oddsratio"],
            ascending=[True, True, False]
        )

    outfile = outdir / f"{set_name}_kegg_ko_enrichment.tsv"
    res.to_csv(outfile, sep="\t", index=False)

    print()
    print(set_name)
    print("tested KOs:", len(res))
    print("significant FDR<0.05:", (res["padj"] < 0.05).sum() if len(res) else 0)
    print(res.head(10))
