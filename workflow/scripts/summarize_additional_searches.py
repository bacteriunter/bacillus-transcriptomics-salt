#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

input_dir = Path("metadata/study_selection/searches/results")
output_dir = Path("metadata/study_selection/searches")
output_dir.mkdir(parents=True, exist_ok=True)

all_rows = []

for file in sorted(input_dir.glob("*.csv")):
    try:
        df = pd.read_csv(file)
        if "Run" not in df.columns or df.empty:
            continue
        df["SearchFile"] = file.name
        all_rows.append(df)
    except Exception as e:
        print(f"Skipping {file}: {e}")

if not all_rows:
    raise SystemExit("No valid search results found.")

df = pd.concat(all_rows, ignore_index=True)

df = df.drop_duplicates(subset=["Run"])

summary = (
    df.groupby(["BioProject", "SRAStudy", "ScientificName"], dropna=False)
      .agg(
          n_runs=("Run", "count"),
          runs=("Run", lambda x: ",".join(x.astype(str))),
          library_layout=("LibraryLayout", lambda x: ",".join(sorted(set(x.astype(str))))),
          platform=("Platform", lambda x: ",".join(sorted(set(x.astype(str))))),
          model=("Model", lambda x: ",".join(sorted(set(x.astype(str))))),
          search_files=("SearchFile", lambda x: ",".join(sorted(set(x.astype(str)))))
      )
      .reset_index()
      .sort_values(["ScientificName", "BioProject"])
)

summary.to_csv(output_dir / "additional_search_summary.tsv", sep="\t", index=False)

print(f"Total unique runs: {df['Run'].nunique()}")
print(f"Total BioProjects: {summary['BioProject'].nunique()}")
print(f"Saved: {output_dir / 'additional_search_summary.tsv'}")
