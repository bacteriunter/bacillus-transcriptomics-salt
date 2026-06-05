#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

outdir = Path("results/tables/main_results")
outdir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(
    "results/tables/master_conserved_ortholog_deg_table.tsv",
    sep="\t"
)

df["combined_abs_log2FC"] = df["ve_log2FC"].abs() + df["pa_log2FC"].abs()

# Tabla 1: Top conserved_up
top_up = (
    df[df["response_class"] == "conserved_up"]
    .sort_values("combined_abs_log2FC", ascending=False)
    .head(50)
)

top_up.to_csv(
    outdir / "Table_1_top50_conserved_up.tsv",
    sep="\t",
    index=False
)

# Tabla 2: Top conserved_down
top_down = (
    df[df["response_class"] == "conserved_down"]
    .sort_values("combined_abs_log2FC", ascending=False)
    .head(50)
)

top_down.to_csv(
    outdir / "Table_2_top50_conserved_down.tsv",
    sep="\t",
    index=False
)

# Tabla 3: candidatos osmoadaptativos
osm_patterns = [
    "sodium_homeostasis",
    "proline",
    "glycine_betaine_choline",
    "poly_gamma_glutamate"
]

osm = df[
    df["functional_category"]
    .fillna("")
    .apply(lambda x: any(p in x for p in osm_patterns))
].copy()

osm = osm.sort_values(
    ["response_class", "functional_category", "combined_abs_log2FC"],
    ascending=[True, True, False]
)

osm.to_csv(
    outdir / "Table_3_osmoadaptation_candidates.tsv",
    sep="\t",
    index=False
)

# Resumen
summary = pd.DataFrame({
    "table": [
        "Table_1_top50_conserved_up",
        "Table_2_top50_conserved_down",
        "Table_3_osmoadaptation_candidates"
    ],
    "n_rows": [
        len(top_up),
        len(top_down),
        len(osm)
    ]
})

summary.to_csv(
    outdir / "main_result_tables_summary.tsv",
    sep="\t",
    index=False
)

print(summary)
