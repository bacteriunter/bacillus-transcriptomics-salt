# Snakemake workflow for:
# Conserved transcriptomic responses to salt stress in Bacillus
#
# This workflow assumes FASTQ files and reference files are already available.
# Heavy intermediate files such as raw FASTQ, Salmon indices, and Salmon quantification
# outputs are excluded from GitHub through .gitignore.

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
    "manuscript/supplementary_tables/TS8_KEGG_KO_Enrichment.tsv"
]

MAIN_FIGURES = [
    "manuscript/figures/Figure1_VolcanoPlots_Combined.png",
    "manuscript/figures/Figure2_ResponseClassSummary.png",
    "manuscript/figures/Figure3_FunctionalCategories.png",
    "manuscript/figures/Figure4_Top50ConservedOrthologsHeatmap.png",
    "manuscript/figures/Figure5_OsmoadaptationHeatmapTop30.png",
    "manuscript/figures/Figure6_KeyOsmoadaptationGenes.png"
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
        ts8 = "results/tables/TS8_KEGG_KO_Enrichment.tsv"
    output:
        ts1 = "manuscript/supplementary_tables/TS1_Bvelezensis_All_DEGs.tsv",
        ts2 = "manuscript/supplementary_tables/TS2_Bparalicheniformis_All_DEGs.tsv",
        ts3 = "manuscript/supplementary_tables/TS3_All_RBH_Orthologs.tsv",
        ts4 = "manuscript/supplementary_tables/TS4_Master_Conserved_Orthologs.tsv",
        ts5 = "manuscript/supplementary_tables/TS5_Top50_ConservedUp.tsv",
        ts6 = "manuscript/supplementary_tables/TS6_Top50_ConservedDown.tsv",
        ts7 = "manuscript/supplementary_tables/TS7_Osmoadaptation_Candidates.tsv",
        ts8 = "manuscript/supplementary_tables/TS8_KEGG_KO_Enrichment.tsv"
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
        """


rule copy_main_figures:
    input:
        fig1 = "figures/manuscript/Figure1_VolcanoPlots_Combined.png",
        fig2 = "figures/response_summary/response_class_summary_barplot.png",
        fig3 = "figures/functional_categories_conserved_response.png",
        fig4 = "figures/heatmap/top50_conserved_orthologs_heatmap.png",
        fig5 = "figures/osmoadaptation/Osmoadaptation_Heatmap_Top30.png",
        fig6 = "figures/osmoadaptation/Osmoadaptation_KeyGenes_Barplot.png"
    output:
        fig1 = "manuscript/figures/Figure1_VolcanoPlots_Combined.png",
        fig2 = "manuscript/figures/Figure2_ResponseClassSummary.png",
        fig3 = "manuscript/figures/Figure3_FunctionalCategories.png",
        fig4 = "manuscript/figures/Figure4_Top50ConservedOrthologsHeatmap.png",
        fig5 = "manuscript/figures/Figure5_OsmoadaptationHeatmapTop30.png",
        fig6 = "manuscript/figures/Figure6_KeyOsmoadaptationGenes.png"
    shell:
        """
        mkdir -p manuscript/figures
        cp {input.fig1} {output.fig1}
        cp {input.fig2} {output.fig2}
        cp {input.fig3} {output.fig3}
        cp {input.fig4} {output.fig4}
        cp {input.fig5} {output.fig5}
        cp {input.fig6} {output.fig6}
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
