# Load necessary libraries
library(pheatmap)
library(RColorBrewer)

# Convert GO terms to human-readable terms
glukozinolates <- c("GO:0019762", "GO:0019761")
glukozinolates <- TermConversion(glukozinolates)
sulfur <- c("GO:0044272", "GO:0072348")
sulfur <- TermConversion(sulfur)

# Load the normalized expression data
ETR_voom <- TMM_voom_counts_ETR
WT_voom <- TMM_voom_counts_WT

ETR_voom_counts <- voom_limma_ETR$E
WT_voom_counts <- voom_limma_WT$E

# Load the normalized mean data
ETR_mean <- normalized_mean_ETR
WT_mean <- normalized_mean_WT

# Combine data for both conditions
unified_genes <- union(rownames(WT_mean), rownames(ETR_mean))
counts_merged <- cbind(ETR_mean[unified_genes, ], WT_mean[unified_genes, ])
colnames(counts_merged) <- c("ETR1", "ETR1E", "WT", "WTE")

# Combine GO term lists
ETR_GO <- genes_GO_limma_ETR
WT_GO <- genes_GO_limma_WT
GO_union <- union(names(WT_GO), names(ETR_GO))
merged_GO <- append(ETR_GO, WT_GO)

# Additional genes for a specific GO term
merged_GO[["glucosinolate catabolic process"]] <- c(merged_GO[["glucosinolate catabolic process"]], "AT3G16400", "AT3G16410")

# Function to create and save heatmaps for each GO term
CreateSetOfPlots <- function(normcounts, golist, deresult_arenosa, deresult_thaliana, gogenes, savedir, root = TRUE) {
  for (i in 1:length(golist)) {
    print(golist[i])
    all_genes <- as.vector(unlist(gogenes[golist[i]]))
    sig_genes_arenosa <- intersect(rownames(deresult_arenosa[deresult_arenosa["adj.P.Val"] <= 0.05, ]), all_genes)
    sig_genes_thaliana <- intersect(rownames(deresult_thaliana[deresult_thaliana["adj.P.Val"] <= 0.05, ]), all_genes)
    sig_genes <- union(sig_genes_arenosa, sig_genes_thaliana)
    sig_counts <- normcounts[sig_genes, ]
    goname <- gsub(":", "", golist[i])
    
    subdir <- savedir
    fname <- paste0("./heatmaps/2024/", subdir, "/", goname, ".png")
    print("All genes:\n")
    print(all_genes)
    print(fname)
    
    if (nrow(sig_counts) > 0) {
      pheatmap(sig_counts, display_numbers = FALSE, main = goname, filename = fname,
               fontsize = 5, number_color = "black", border_color = "black", cellwidth = 10, cellheight = 10,
               color = colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Set2"))(100)) # Set2 is good for colorblind users
    } else {
      print(paste("Not enough significant genes for", goname))
    }
    
    if (length(dev.list()) > 0) {
      for (dev_id in dev.list()) {
        dev.off(dev_id)
      }
    }
  }
}

# Create and save heatmaps for the specified GO terms
CreateSetOfPlots(counts_merged, glukozinolates, ETR_voom, WT_voom, merged_GO, "glukozinolates")
CreateSetOfPlots(counts_merged, sulfur, ETR_voom, WT_voom, merged_GO, "sulfur")
