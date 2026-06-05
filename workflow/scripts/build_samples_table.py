#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

runinfo_main = pd.read_csv("metadata/sra_tables/bacillus_salt_rnaseq_runinfo.csv")
runinfo_velezensis = pd.read_csv("metadata/sra_tables/PRJNA1402157_full_runinfo.csv")

runinfo = (
    pd.concat([runinfo_main, runinfo_velezensis], ignore_index=True)
      .drop_duplicates(subset=["Run"])
)

manual = [
    # PRJNA1402157 - B. velezensis
    ("PRJNA1402157", "SRR36805042", "Bacillus_velezensis", "HR6-1", "control", "Control1", "core"),
    ("PRJNA1402157", "SRR36805041", "Bacillus_velezensis", "HR6-1", "control", "Control2", "core"),
    ("PRJNA1402157", "SRR36805040", "Bacillus_velezensis", "HR6-1", "control", "Control3", "core"),
    ("PRJNA1402157", "SRR36805039", "Bacillus_velezensis", "HR6-1", "salt", "Salt1", "core"),
    ("PRJNA1402157", "SRR36805038", "Bacillus_velezensis", "HR6-1", "salt", "Salt2", "core"),
    ("PRJNA1402157", "SRR36805037", "Bacillus_velezensis", "HR6-1", "salt", "Salt3", "core"),

    # PRJNA1109437 - B. paralicheniformis
    ("PRJNA1109437", "SRR28973049", "Bacillus_paralicheniformis", "BL-1", "salt", "S1", "core"),
    ("PRJNA1109437", "SRR28973048", "Bacillus_paralicheniformis", "BL-1", "salt", "S2", "core"),
    ("PRJNA1109437", "SRR28973047", "Bacillus_paralicheniformis", "BL-1", "salt", "S3", "core"),
    ("PRJNA1109437", "SRR28973052", "Bacillus_paralicheniformis", "BL-1", "control", "NS1", "core"),
    ("PRJNA1109437", "SRR28973051", "Bacillus_paralicheniformis", "BL-1", "control", "NS2", "core"),
    ("PRJNA1109437", "SRR28973050", "Bacillus_paralicheniformis", "BL-1", "control", "NS3", "core"),

    # PRJNA252753 - B. licheniformis
    ("PRJNA252753", "SRR1565402", "Bacillus_licheniformis", "WX-02", "control", "normal", "core_limited"),
    ("PRJNA252753", "SRR1565404", "Bacillus_licheniformis", "WX-02", "salt", "NaCl_6_percent", "core_limited"),

    # PRJNA321031 - B. subtilis
    ("PRJNA321031", "SRR3488624", "Bacillus_subtilis", "168", "salt", "salt_T30_rep1", "complementary"),
    ("PRJNA321031", "SRR3488629", "Bacillus_subtilis", "168", "salt", "salt_T30_rep2", "complementary"),
    ("PRJNA321031", "SRR3488625", "Bacillus_subtilis", "168", "salt", "salt_T60_rep1", "complementary"),
    ("PRJNA321031", "SRR3488634", "Bacillus_subtilis", "168", "salt", "salt_T60_rep2", "complementary"),
    ("PRJNA321031", "SRR3488626", "Bacillus_subtilis", "168", "salt", "salt_T90_rep1", "complementary"),
    ("PRJNA321031", "SRR3488635", "Bacillus_subtilis", "168", "salt", "salt_T90_rep2", "complementary"),
    ("PRJNA321031", "SRR3488628", "Bacillus_subtilis", "168", "control", "nosalt_T30_rep1", "complementary"),
    ("PRJNA321031", "SRR3488631", "Bacillus_subtilis", "168", "control", "nosalt_T30_rep2", "complementary"),
    ("PRJNA321031", "SRR3488630", "Bacillus_subtilis", "168", "control", "nosalt_T60_rep1", "complementary"),
    ("PRJNA321031", "SRR3488632", "Bacillus_subtilis", "168", "control", "nosalt_T60_rep2", "complementary"),
    ("PRJNA321031", "SRR3488633", "Bacillus_subtilis", "168", "control", "nosalt_T90_rep1", "complementary"),
]

samples = pd.DataFrame(
    manual,
    columns=["bioproject", "run_accession", "species", "strain", "condition", "sample_id", "role"]
)

samples = samples.merge(
    runinfo[["Run", "Experiment", "BioSample", "LibraryLayout", "Platform", "Model"]],
    left_on="run_accession",
    right_on="Run",
    how="left"
)

samples = samples.drop(columns=["Run"])
samples = samples.rename(columns={
    "Experiment": "experiment_accession",
    "BioSample": "biosample",
    "LibraryLayout": "library_layout",
    "Platform": "platform",
    "Model": "sequencer_model",
})


def assign_comparison_group(row):
    if row["bioproject"] == "PRJNA1402157":
        return "velezensis_control_vs_salt"
    if row["bioproject"] == "PRJNA1109437":
        return "paralicheniformis_control_vs_salt"
    if row["bioproject"] == "PRJNA252753":
        return "licheniformis_control_vs_salt"
    if row["bioproject"] == "PRJNA321031":
        return "subtilis_timecourse"
    return "unknown"

samples["comparison_group"] = samples.apply(assign_comparison_group, axis=1)

samples["fastq_1"] = "raw_data/fastq/" + samples["run_accession"] + "_1.fastq.gz"
samples["fastq_2"] = "raw_data/fastq/" + samples["run_accession"] + "_2.fastq.gz"

out = Path("metadata/sample_annotation/samples.tsv")
out.parent.mkdir(parents=True, exist_ok=True)
samples.to_csv(out, sep="\t", index=False)

print(f"Saved: {out}")
print("Total samples:", len(samples))
print(samples[["bioproject", "species", "condition", "sample_id", "run_accession", "role"]])
