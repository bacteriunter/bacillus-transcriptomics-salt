shell.executable("/bin/bash")

configfile: "config.yaml"

MAIN_TABLES = [
    "manuscript/tables/Table1_Transcriptome_DEG_Summary.tsv",
    "manuscript/tables/Table2_Core_Osmoadaptation_Genes.tsv"
]

SUPPLEMENTARY_TABLES = [
    "manuscript/supplementary_tables/TS1_Bvelezensis_All_DEGs.tsv",
    "manuscript/supplementary_tables/TS2_Bparalicheniformis_All_DEGs.tsv",
    "manuscript/supplementary_tables/TS3_All_RBH_Orthologs.tsv",
    "manuscript/supplementary_tables/TS4_Master_Conserved_Orthologs.tsv",
    "manuscript/supplementary_tables/TS5_Top50_ConservedUp.tsv",
    "manuscript/supplementary_tables/TS6_Top50_ConservedDown.tsv",
    "manuscript/supplementary_tables/TS7_Osmoadaptation_Candidates.tsv",
    "manuscript/supplementary_tables/TS8_KEGG_KO_Enrichment.tsv",
    "manuscript/supplementary_tables/TS9_Opposite_Category_Enrichment.tsv"
]

MAIN_FIGURES = [
    "manuscript/figures/Figure1_VolcanoPlots_Combined.png",
    "manuscript/figures/Figure1_VolcanoPlots_Combined.pdf",
    "manuscript/figures/Figure2_ResponseClassSummary.png",
    "manuscript/figures/Figure2_ResponseClassSummary.pdf",
    "manuscript/figures/Figure3_FunctionalCategories.png",
    "manuscript/figures/Figure3_FunctionalCategories.pdf",
    "manuscript/figures/Figure4_Top50ConservedOrthologsHeatmap.png",
    "manuscript/figures/Figure4_Top50ConservedOrthologsHeatmap.pdf",
    "manuscript/figures/Figure5_OsmoadaptationHeatmapTop30.png",
    "manuscript/figures/Figure5_OsmoadaptationHeatmapTop30.pdf",
    "manuscript/figures/Figure6_KeyOsmoadaptationGenes.png",
    "manuscript/figures/Figure6_KeyOsmoadaptationGenes.pdf",
    "manuscript/figures/Figure7_OrthologLog2FCCorrelation.png",
    "manuscript/figures/Figure7_OrthologLog2FCCorrelation.pdf",
    "manuscript/figures/Figure8_Conceptual_Model_Osmoadaptation.png",
    "manuscript/figures/Figure8_Conceptual_Model_Osmoadaptation.pdf"
]

SUPPLEMENTARY_FIGURES = [
    "manuscript/supplementary_figures/FigureS1_PCA_Combined.png"
]

rule all:
    input:
        MAIN_TABLES,
        SUPPLEMENTARY_TABLES,
        MAIN_FIGURES,
        SUPPLEMENTARY_FIGURES


rule table1_transcriptome_summary:
    input:
        ve_summary = "results/deseq2/velezensis_control_vs_salt/deseq2_summary.tsv",
        ve_filtered = "results/deseq2/velezensis_control_vs_salt/filtered/filtered_deg_summary.tsv",
        pa_summary = "results/deseq2/paralicheniformis_control_vs_salt/deseq2_summary.tsv",
        pa_filtered = "results/deseq2/paralicheniformis_control_vs_salt/filtered/filtered_deg_summary.tsv"
    output:
        "results/tables/Table1_Transcriptome_DEG_Summary.tsv"
    shell:
        "python workflow/scripts/build_table1_dataset_summary.py"


rule table2_core_osmoadaptation:
    input:
        "results/tables/Table3_Osmoadaptation_Candidates.tsv"
    output:
        "results/tables/Table2_Core_Osmoadaptation_Genes.tsv"
    shell:
        "python workflow/scripts/build_table2_core_osmoadaptation.py"


rule ts3_all_rbh:
    input:
        "results/orthology/mmseqs_rbh/velezensis_vs_paralicheniformis_rbh.tsv"
    output:
        "results/tables/TS3_All_RBH_Orthologs.tsv"
    shell:
        "python workflow/scripts/build_TS3_RBH_table.py"


rule ts8_kegg_enrichment:
    input:
        up = "results/enrichment/kegg_ko/conserved_up_kegg_ko_enrichment.tsv",
        down = "results/enrichment/kegg_ko/conserved_down_kegg_ko_enrichment.tsv",
        opposite = "results/enrichment/kegg_ko/opposite_kegg_ko_enrichment.tsv"
    output:
        "results/tables/TS8_KEGG_KO_Enrichment.tsv"
    shell:
        "python workflow/scripts/build_TS8_KEGG_enrichment.py"


rule ts9_opposite_category_enrichment:
    input:
        "results/tables/master_conserved_ortholog_deg_table.tsv"
    output:
        ts9 = "results/tables/TS9_Opposite_Category_Enrichment.tsv",
        enrichment = "results/tables/opposite_category_enrichment.tsv",
        summary = "results/tables/opposite_functional_category_summary.tsv",
        topko = "results/tables/opposite_top_KO_definitions.tsv"
    shell:
        """
        /home/bacte/miniconda3/envs/docking/bin/python - << 'PY'
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

df = pd.read_csv("results/tables/master_conserved_ortholog_deg_table.tsv", sep="\\t")
opp = df[df["response_class"] == "opposite"].copy()

cats = (
    opp["functional_category"]
    .fillna("unclassified")
    .str.split(";")
    .explode()
    .value_counts()
    .rename_axis("functional_category")
    .reset_index(name="count")
)
cats.to_csv("results/tables/opposite_functional_category_summary.tsv", sep="\\t", index=False)

topko = (
    opp["KO_definition"]
    .fillna("unannotated")
    .replace("", "unannotated")
    .value_counts()
    .rename_axis("KO_definition")
    .reset_index(name="count")
)
topko.to_csv("results/tables/opposite_top_KO_definitions.tsv", sep="\\t", index=False)

categories = set()
for x in df["functional_category"].fillna("unclassified"):
    for c in str(x).split(";"):
        categories.add(c)

rows = []
for cat in sorted(categories):
    in_opp = df["response_class"] == "opposite"
    has_cat = df["functional_category"].fillna("").str.contains(cat, regex=False)

    a = ((in_opp) & (has_cat)).sum()
    b = ((in_opp) & (~has_cat)).sum()
    c = ((~in_opp) & (has_cat)).sum()
    d = ((~in_opp) & (~has_cat)).sum()

    odds, p = fisher_exact([[a, b], [c, d]])
    rows.append([cat, a, c, odds, p])

res = pd.DataFrame(
    rows,
    columns=[
        "functional_category",
        "opposite_orthologs",
        "non_opposite_orthologs",
        "odds_ratio",
        "pvalue"
    ]
)

res["padj"] = multipletests(res["pvalue"], method="fdr_bh")[1]
res = res.sort_values(["padj", "pvalue"])

res.to_csv("results/tables/opposite_category_enrichment.tsv", sep="\\t", index=False)
res.to_csv("results/tables/TS9_Opposite_Category_Enrichment.tsv", sep="\\t", index=False)
PY
        """


rule copy_main_tables:
    input:
        table1 = "results/tables/Table1_Transcriptome_DEG_Summary.tsv",
        table2 = "results/tables/Table2_Core_Osmoadaptation_Genes.tsv"
    output:
        table1 = "manuscript/tables/Table1_Transcriptome_DEG_Summary.tsv",
        table2 = "manuscript/tables/Table2_Core_Osmoadaptation_Genes.tsv"
    shell:
        """
        mkdir -p manuscript/tables
        cp {input.table1} {output.table1}
        cp {input.table2} {output.table2}
        """


rule copy_supplementary_tables:
    input:
        ts1 = "results/deseq2/velezensis_control_vs_salt/filtered/DEGs_padj0.05_log2FC1.tsv",
        ts2 = "results/deseq2/paralicheniformis_control_vs_salt/filtered/DEGs_padj0.05_log2FC1.tsv",
        ts3 = "results/tables/TS3_All_RBH_Orthologs.tsv",
        ts4 = "results/tables/master_conserved_ortholog_deg_table.tsv",
        ts5 = "results/tables/Table1_Top50_ConservedUp.tsv",
        ts6 = "results/tables/Table2_Top50_ConservedDown.tsv",
        ts7 = "results/tables/Table3_Osmoadaptation_Candidates.tsv",
        ts8 = "results/tables/TS8_KEGG_KO_Enrichment.tsv",
        ts9 = "results/tables/TS9_Opposite_Category_Enrichment.tsv"
    output:
        ts1 = "manuscript/supplementary_tables/TS1_Bvelezensis_All_DEGs.tsv",
        ts2 = "manuscript/supplementary_tables/TS2_Bparalicheniformis_All_DEGs.tsv",
        ts3 = "manuscript/supplementary_tables/TS3_All_RBH_Orthologs.tsv",
        ts4 = "manuscript/supplementary_tables/TS4_Master_Conserved_Orthologs.tsv",
        ts5 = "manuscript/supplementary_tables/TS5_Top50_ConservedUp.tsv",
        ts6 = "manuscript/supplementary_tables/TS6_Top50_ConservedDown.tsv",
        ts7 = "manuscript/supplementary_tables/TS7_Osmoadaptation_Candidates.tsv",
        ts8 = "manuscript/supplementary_tables/TS8_KEGG_KO_Enrichment.tsv",
        ts9 = "manuscript/supplementary_tables/TS9_Opposite_Category_Enrichment.tsv"
    shell:
        """
        mkdir -p manuscript/supplementary_tables
        cp {input.ts1} {output.ts1}
        cp {input.ts2} {output.ts2}
        cp {input.ts3} {output.ts3}
        cp {input.ts4} {output.ts4}
        cp {input.ts5} {output.ts5}
        cp {input.ts6} {output.ts6}
        cp {input.ts7} {output.ts7}
        cp {input.ts8} {output.ts8}
        cp {input.ts9} {output.ts9}
        """


rule figure1_volcano:
    input:
        ve = "results/volcano/velezensis_control_vs_salt_volcano_table.tsv",
        pa = "results/volcano/paralicheniformis_control_vs_salt_volcano_table.tsv"
    output:
        png = "manuscript/figures/Figure1_VolcanoPlots_Combined.png",
        pdf = "manuscript/figures/Figure1_VolcanoPlots_Combined.pdf"
    shell:
        """
        source ~/miniconda3/etc/profile.d/conda.sh
        conda activate deseq2_env
        Rscript workflow/scripts/build_figure1_volcano_publication.R
        """


rule figure2_response_classes:
    input:
        "results/orthology/response_classes/response_class_summary.tsv"
    output:
        png = "manuscript/figures/Figure2_ResponseClassSummary.png",
        pdf = "manuscript/figures/Figure2_ResponseClassSummary.pdf"
    shell:
        """
        source ~/miniconda3/etc/profile.d/conda.sh
        conda activate deseq2_env
        Rscript workflow/scripts/build_figure2_response_class_publication.R
        """


rule figure3_functional_categories:
    input:
        "results/tables/functional_category_response_matrix.tsv"
    output:
        png = "manuscript/figures/Figure3_FunctionalCategories.png",
        pdf = "manuscript/figures/Figure3_FunctionalCategories.pdf"
    shell:
        """
        source ~/miniconda3/etc/profile.d/conda.sh
        conda activate deseq2_env
        Rscript workflow/scripts/build_figure3_functional_categories_publication.R
        """


rule figure4_top50_heatmap:
    input:
        "results/tables/master_conserved_ortholog_deg_table.tsv"
    output:
        png = "manuscript/figures/Figure4_Top50ConservedOrthologsHeatmap.png",
        pdf = "manuscript/figures/Figure4_Top50ConservedOrthologsHeatmap.pdf"
    shell:
        """
        source ~/miniconda3/etc/profile.d/conda.sh
        conda activate deseq2_env
        Rscript workflow/scripts/plot_top50_conserved_heatmap.R
        """


rule figure5_osmoadaptation_heatmap:
    input:
        "results/tables/Table3_Osmoadaptation_Candidates.tsv"
    output:
        png = "manuscript/figures/Figure5_OsmoadaptationHeatmapTop30.png",
        pdf = "manuscript/figures/Figure5_OsmoadaptationHeatmapTop30.pdf"
    shell:
        """
        source ~/miniconda3/etc/profile.d/conda.sh
        conda activate deseq2_env
        Rscript workflow/scripts/plot_osmoadaptation_candidates_heatmap_publication.R
        """


rule figure6_keygenes_barplot:
    input:
        "results/tables/Table2_Core_Osmoadaptation_Genes.tsv"
    output:
        png = "manuscript/figures/Figure6_KeyOsmoadaptationGenes.png",
        pdf = "manuscript/figures/Figure6_KeyOsmoadaptationGenes.pdf"
    shell:
        """
        source ~/miniconda3/etc/profile.d/conda.sh
        conda activate deseq2_env
        Rscript workflow/scripts/plot_osmoadaptation_keygenes_barplot_publication.R
        """


rule figure7_ortholog_correlation:
    input:
        "results/tables/master_conserved_ortholog_deg_table.tsv"
    output:
        png = "manuscript/figures/Figure7_OrthologLog2FCCorrelation.png",
        pdf = "manuscript/figures/Figure7_OrthologLog2FCCorrelation.pdf",
        summary = "results/correlation/ortholog_log2FC_correlation_summary.tsv",
        scatter = "results/correlation/ortholog_log2FC_scatter_data.tsv"
    shell:
        """
        source ~/miniconda3/etc/profile.d/conda.sh
        conda activate deseq2_env
        Rscript workflow/scripts/plot_ortholog_log2fc_correlation.R
        cp figures/correlation/Figure_Ortholog_log2FC_Correlation.png {output.png}
        cp figures/correlation/Figure_Ortholog_log2FC_Correlation.pdf {output.pdf}
        """


rule figure8_conceptual_model:
    input:
        "results/correlation/ortholog_log2FC_correlation_summary.tsv",
        "results/tables/TS9_Opposite_Category_Enrichment.tsv"
    output:
        png = "manuscript/figures/Figure8_Conceptual_Model_Osmoadaptation.png",
        pdf = "manuscript/figures/Figure8_Conceptual_Model_Osmoadaptation.pdf"
    shell:
        """
        /home/bacte/miniconda3/envs/docking/bin/python workflow/scripts/build_conceptual_model_figure.py
        cp figures/conceptual_model/Figure8_Conceptual_Model_Osmoadaptation.png {output.png}
        cp figures/conceptual_model/Figure8_Conceptual_Model_Osmoadaptation.pdf {output.pdf}
        """


rule copy_supplementary_figures:
    input:
        "figures/supplementary/FigureS1_PCA_Combined.png"
    output:
        "manuscript/supplementary_figures/FigureS1_PCA_Combined.png"
    shell:
        """
        mkdir -p manuscript/supplementary_figures
        cp {input} {output}
        """
