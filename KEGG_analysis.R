# Load necessary libraries
library(KEGGprofile)
library(limma)

# Load differential expression results
source("DE_analysis.R")

#---- KEGG Pathway Analysis for ETR1 vs ETR1E

# Prepare gene list for KEGG analysis
ni_genes_ETR <- rownames(significant_limma_ETR)

# Perform KEGG pathway enrichment analysis
ni_KEGG_ETR <- find_enriched_pathway(ni_genes_ETR, species = "ath", returned_genenumber = 2, download_latest = FALSE, returned_pvalue = 0.05)
write.csv(ni_KEGG_ETR$stastic, "./results/KEGG/KEGG_thaliana_statistic_ETR.csv")

# Prepare normalized expression data for plotting
expr_ETR <- normalized_counts_ETR
colnames(expr_ETR) <- c("ETR1_1", "ETR1_2", "ETR1_3", "ETR1_4", "ETR1_5",
                        "ETR1E_1", "ETR1E_2", "ETR1E_3", "ETR1E_4", "ETR1E_5")

# Plot correlation of enriched pathways
plot_pathway_cor(gene_expr = expr_ETR, kegg_enriched_pathway = ni_KEGG_ETR)

# Prepare data for pathway map plotting
E_expr_ETR <- expr_ETR[, 6:10]
temp_ETR <- apply(E_expr_ETR, 1, function(x) length(which(is.na(x))))
E_expr_ETR <- E_expr_ETR[which(temp_ETR == 0), ]
col <- col_by_value(E_expr_ETR, col = colorRampPalette(c("blue", "violet", "red"))(1024), range = c(-6, 6))

# Plot specific KEGG pathways
KEGGdf_ETR <- ni_KEGG_ETR$stastic
for (i in 1:nrow(KEGGdf_ETR)) {
  out <- tryCatch({
    KEGG_n <- rownames(KEGGdf_ETR)[i]
    KEGG_id <- KEGGdf_ETR[i, 1]
    download_KEGGfile(pathway_id = KEGG_n, species = "ath", target_dir = getwd())
    KEGG_database <- parse_XMLfile(KEGG_n, species = "ath", database_dir = getwd())
    plot_profile(gene_expr = E_expr_ETR, pathway_name = paste0("ath", KEGG_n), result_name = paste0("./plots/KEGG/ETR/", KEGG_id, ".png"),
                 KEGG_database = KEGG_database, text_col = "white", type = "bg", bg_col = col,
                 species = "ath", magnify = 1.2, pathway_min = 0)
  }, error = function(cond) {
    print("Problem with XML or too global perspective")
  }, warning = function(cond) {
    print("Problem with XML or too global perspective")
  }, finally = {
    print("Problem with XML or too global perspective")
  })
}

#---- KEGG Pathway Analysis for WT vs WTE

# Prepare gene list for KEGG analysis
ni_genes_WT <- rownames(significant_limma_WT)

# Perform KEGG pathway enrichment analysis
ni_KEGG_WT <- find_enriched_pathway(ni_genes_WT, species = "ath", returned_genenumber = 1, download_latest = FALSE, returned_pvalue = 0.05)
write.csv(ni_KEGG_WT$stastic, "./results/KEGG/KEGG_thaliana_statistic_WT.csv")

# Prepare normalized expression data for plotting
expr_WT <- normalized_counts_WT
colnames(expr_WT) <- c("WT_1", "WT_2", "WT_3", "WT_4", "WT_5",
                       "WTE_1", "WTE_2", "WTE_3", "WTE_4", "WTE_5")

# Plot correlation of enriched pathways
plot_pathway_cor(gene_expr = expr_WT, kegg_enriched_pathway = ni_KEGG_WT)

# Prepare data for pathway map plotting
E_expr_WT <- expr_WT[, 6:10]
temp_WT <- apply(E_expr_WT, 1, function(x) length(which(is.na(x))))
E_expr_WT <- E_expr_WT[which(temp_WT == 0), ]
col <- col_by_value(E_expr_WT, col = colorRampPalette(c("blue", "violet", "red"))(1024), range = c(-6, 6))

# Plot specific KEGG pathways
KEGGdf_WT <- ni_KEGG_WT$stastic
for (i in 1:nrow(KEGGdf_WT)) {
  out <- tryCatch({
    KEGG_n <- rownames(KEGGdf_WT)[i]
    KEGG_id <- KEGGdf_WT[i, 1]
    download_KEGGfile(pathway_id = KEGG_n, species = "ath", target_dir = getwd())
    KEGG_database <- parse_XMLfile(KEGG_n, species = "ath", database_dir = getwd())
    plot_profile(gene_expr = E_expr_WT, pathway_name = paste0("ath", KEGG_n), result_name = paste0("./plots/KEGG/WT/", KEGG_id, ".png"),
                 KEGG_database = KEGG_database, text_col = "white", type = "bg", bg_col = col,
                 species = "ath", magnify = 1.2, pathway_min = 0)
  }, error = function(cond) {
    print("Problem with XML or too global perspective")
  }, warning = function(cond) {
    print("Problem with XML or too global perspective")
  }, finally = {
    print("Problem with XML or too global perspective")
  })
}
