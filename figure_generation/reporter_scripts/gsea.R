#The defined comparisons were generated in the edgeR.R script. This script was also ran for the 50_v_0 comparison.

#Preparing a ranked list for Gene Set Enrichment Analysis (GSEA) based on logFC and F statistic
res_20_v_0 <- results[["20_vs_0"]]
gene_list_20_v_0 <- sign(res_20_v_0$logFC) * res_20_v_0$F
names(gene_list_20_v_0) <- rownamesres_20_v_0)
gene_list_20_v_0 <- gene_list_20_v_0[!is.na(gene_list_20_v_0)]
gene_list_20_v_0 <- sort(gene_list_20_v_0, decreasing = TRUE)

library(clusterProfiler) #(v4.18.4)
library(org.Mm.eg.db) #(v3.22.0)

#GSEA of gene ontology (biological processes, mouse)
gsea_go_20_v_0 <- gseGO(
  geneList = gene_list_20_v_0,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pvalueCutoff = 1,
  pAdjustMethod = "BH"
)

#Extract pathways table
gsea_go_20_v_0_df <- as.data.frame(gsea_go_20_v_0)

#Create a master pathway list
master_pathways_20_v_0 <- gsea_go_20_v_0_df |>
  dplyr::select(
    GO_ID = ID,
    Pathway = Description,
    setSize,
    enrichmentScore,
    NES,
    pvalue,
    p.adjust,
    qvalue,
    rank,
    leading_edge,
    core_enrichment
  )

#Write to csv
write.csv(
  master_pathways_20_v_0,
  file = "gsea_go_20_v_0_master_pathways.csv",
  row.names = FALSE
)

#Plot net enrichment scores (NES) and adjusted p-value of select enriched pathways
#Subsetting select pathways
selected_pathways <- c(
  "regulation of oligodendrocyte differentiation",
  "urogenital system development",
  "smoothened signaling pathway"
)

plot_df <- gsea_go_20_v_0_df %>%
  select(Description, NES, p.adjust) %>%
  filter(Description %in% selected_pathways) 

library(ggplot2) #(v4.0.2)
library(dplyr) #(v1.2.0)
library(stringr) #(v1.6.0)

# Wrap long pathway names
plot_df$Description <- str_wrap(plot_df$Description, width = 30)

ggplot(plot_df, aes(x = NES,
                    y = reorder(Description, NES))) +
  
  geom_point(aes(size = -log10(p.adjust),
                 color = NES)) +
  
  scale_color_gradient(low = "#9e9ac8",
                       high = "#54278f",
                       name = "NES") +
  
  scale_size_continuous(name = "-log10(adj p)") +
  
  labs(
    title = "Generation 20 Relative to Generation 0",
    x = "Normalized Enrichment Score (NES)",
    y = "GO Pathway"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank()
  )
