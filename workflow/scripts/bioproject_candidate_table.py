#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

runinfo_file = Path("metadata/sra_tables/bacillus_salt_rnaseq_runinfo.csv")
output_file = Path("metadata/study_selection/bioproject_candidate_table.tsv")

keep_bioprojects = [
    "PRJNA1029469",
    "PRJNA1109437",
    "PRJNA1402157",
    "PRJNA194121",
    "PRJNA252753",
    "PRJNA321031",
    "PRJNA386487",
    "PRJNA431298",
    "PRJNA511580",
    "PRJNA675787",
    "PRJNA744162",
    "PRJNA854068",
    "PRJNA928739",
    "PRJNA1425010",
]

df = pd.read_csv(runinfo_file)
df = df[df["BioProject"].isin(keep_bioprojects)]

summary = (
    df.groupby(["BioProject", "SRAStudy", "ScientificName"], dropna=False)
      .agg(
          n_runs=("Run", "count"),
          runs=("Run", lambda x: ",".join(x.astype(str))),
          experiments=("Experiment", lambda x: ",".join(x.astype(str).unique())),
          biosamples=("BioSample", lambda x: ",".join(x.astype(str).unique())),
          platform=("Platform", lambda x: ",".join(x.astype(str).unique())),
          model=("Model", lambda x: ",".join(x.astype(str).unique())),
          layout=("LibraryLayout", lambda x: ",".join(x.astype(str).unique())),
          total_size_MB=("size_MB", "sum")
      )
      .reset_index()
      .sort_values(["ScientificName", "BioProject"])
)

summary["decision"] = "pending"
summary["reason"] = ""

output_file.parent.mkdir(parents=True, exist_ok=True)
summary.to_csv(output_file, sep="\t", index=False)

print(f"Saved: {output_file}")
print(f"BioProjects summarized: {summary['BioProject'].nunique()}")
print(f"Rows: {len(summary)}")
