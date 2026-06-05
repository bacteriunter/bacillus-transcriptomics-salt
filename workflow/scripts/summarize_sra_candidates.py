#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

input_file = Path("metadata/sra_tables/bacillus_salt_rnaseq_runinfo.csv")
output_file = Path("metadata/sra_tables/bacillus_salt_rnaseq_summary.tsv")

df = pd.read_csv(input_file)

summary = (
    df.groupby(["BioProject", "SRAStudy", "ScientificName", "LibraryLayout", "Platform", "Model"], dropna=False)
      .agg(
          n_runs=("Run", "count"),
          total_size_MB=("size_MB", "sum"),
          first_release=("ReleaseDate", "min"),
          last_release=("ReleaseDate", "max"),
          runs=("Run", lambda x: ",".join(x.astype(str)))
      )
      .reset_index()
      .sort_values(["ScientificName", "BioProject"])
)

summary.to_csv(output_file, sep="\t", index=False)

print(f"Saved: {output_file}")
print(f"Number of summarized groups: {len(summary)}")
