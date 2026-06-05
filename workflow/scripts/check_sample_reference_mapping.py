#!/usr/bin/env python3

import pandas as pd

samples = pd.read_csv(
    "metadata/sample_annotation/samples.tsv",
    sep="\t"
)

refs = pd.read_csv(
    "metadata/reference_genomes/sample_reference_map.tsv",
    sep="\t"
)

merged = samples.merge(
    refs,
    on="species",
    how="left"
)

print()
print("Samples:", len(merged))
print("Species:", merged["species"].nunique())
print()

print(
    merged[
        [
            "sample_id",
            "species",
            "selected_accession"
        ]
    ]
)

missing = merged["selected_accession"].isna().sum()

print()
print("Missing references:", missing)

