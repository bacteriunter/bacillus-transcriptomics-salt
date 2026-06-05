#!/usr/bin/env python3

import pandas as pd

keep = [
    "PRJNA252753",
    "PRJNA1109437",
    "PRJNA1402157",
    "PRJNA321031",
    "PRJNA854068"
]

df = pd.read_csv(
    "metadata/sra_tables/bacillus_salt_rnaseq_runinfo.csv"
)

df = df[df["BioProject"].isin(keep)]

cols = [
    "BioProject",
    "ScientificName",
    "Run",
    "Experiment",
    "Sample",
    "BioSample",
    "LibraryLayout",
    "Model"
]

df[cols].to_csv(
    "metadata/study_selection/manual_curation/candidate_runs.tsv",
    sep="\t",
    index=False
)

print(df[cols].shape)
