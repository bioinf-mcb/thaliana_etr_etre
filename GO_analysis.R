# Load necessary libraries
library(topGO)
library(org.At.tair.db)
library(limma)

# Load differential expression results
source("DE_analysis.R")

# Function to convert GO IDs to GO terms
IDtoTerm <- function(GOlist) {
  GO_names <- names(GOlist)
  for (i in 1:length(GO_names)) {
    GO_names[i] <- AnnotationDbi::Term(GOTERM[GO_names[i]])
  }
  names(GOlist) <- GO_names
  return(GOlist)
}

# Function to convert a vector of GO IDs to GO terms
TermConversion <- function(GOvec) {
  for (i in 1:length(GOvec)) {
    GOvec[i] <- AnnotationDbi::Term(GOTERM[GOvec[i]])
  }
  return(GOvec)
}

#---- Gene Ontology Analysis for ETR1 vs ETR1E

# Prepare gene list for topGO analysis
genelist_limma_ETR <- TMM_voom_counts_ETR$adj.P.Val
names(genelist_limma_ETR) <- rownames(TMM_voom_counts_ETR)

# Function to select significant genes
topDiffGenes <- function(allGene, thr = 0.05) {
  return(allGene <= thr)
}

# Create topGOdata object for Biological Process (BP) ontology
topGO_limma_ETR <- new("topGOdata",
                       description = "Limma BP", ontology = "BP",
                       allGenes = genelist_limma_ETR, geneSel = topDiffGenes,
                       nodeSize = 1,
                       annot = annFUN.org, mapping = "org.At.tair.db")

# Run enrichment test using Fisher's exact test with the "elim" algorithm
limma_Fisher_ETR <- runTest(topGO_limma_ETR, algorithm = "elim", statistic = "fisher")

# Generate results table
GO_res_ETR <- GenTable(topGO_limma_ETR, elimFisher = limma_Fisher_ETR,
                       orderBy = "elimFisher", ranksOf = "elimFisher",
                       topNodes = 500)

# Save results to CSV file
write.csv(GO_res_ETR, "ETR_BP.csv")

# Identify differentially expressed genes associated with significant GO terms
changed_genes_ETR <- rownames(TMM_voom_counts_ETR[TMM_voom_counts_ETR$adj.P.Val <= 0.05, ])
genes_GO_limma_ETR <- genesInTerm(topGO_limma_ETR)
genes_GO_limma_ETR <- IDtoTerm(genes_GO_limma_ETR)

for (i in 1:length(genes_GO_limma_ETR)) {
  tmp_vec <- genes_GO_limma_ETR[[i]]
  tmp_vec <- tmp_vec[tmp_vec %in% changed_genes_ETR]
  genes_GO_limma_ETR[[i]] <- tmp_vec
}

# Save gene-GO term associations
save(genes_GO_limma_ETR, file ="genes_GO_limma_ETR.RData")
sink("genes_limma_ETR.txt")
print(genes_GO_limma_ETR)
sink()

#---- Gene Ontology Analysis for WT vs WTE

# Prepare gene list for topGO analysis
genelist_limma_WT <- TMM_voom_counts_WT$adj.P.Val
names(genelist_limma_WT) <- rownames(TMM_voom_counts_WT)

# Create topGOdata object for Biological Process (BP) ontology
topGO_limma_WT <- new("topGOdata",
                      description = "Limma BP", ontology = "BP",
                      allGenes = genelist_limma_WT, geneSel = topDiffGenes,
                      nodeSize = 1,
                      annot = annFUN.org, mapping = "org.At.tair.db")

# Run enrichment test using Fisher's exact test with the "elim" algorithm
limma_Fisher_WT <- runTest(topGO_limma_WT, algorithm = "elim", statistic = "fisher")

# Generate results table
GO_res_WT <- GenTable(topGO_limma_WT, elimFisher = limma_Fisher_WT,
                      orderBy = "elimFisher", ranksOf = "elimFisher",
                      topNodes = 500)

# Save results to CSV file
write.csv(GO_res_WT, "WT_BP.csv")

# Identify differentially expressed genes associated with significant GO terms
changed_genes_WT <- rownames(TMM_voom_counts_WT[TMM_voom_counts_WT$adj.P.Val <= 0.05, ])
genes_GO_limma_WT <- genesInTerm(topGO_limma_WT)
genes_GO_limma_WT <- IDtoTerm(genes_GO_limma_WT)

for (i in 1:length(genes_GO_limma_WT)) {
  tmp_vec <- genes_GO_limma_WT[[i]]
  tmp_vec <- tmp_vec[tmp_vec %in% changed_genes_WT]
  genes_GO_limma_WT[[i]] <- tmp_vec
}

# Save gene-GO term associations
save(genes_GO_limma_WT, file ="genes_GO_limma_WT.RData")
sink("genes_limma_WT.txt")
print(genes_GO_limma_WT)
sink()

#---- Functional Group Analysis

# Define functional groups
auxin <- TermConversion(c("GO:0009733", "GO:0009850"))
root <- TermConversion(c("GO:0048765", "GO:0010054", "GO:0048767", "GO:0048364"))
other <- TermConversion(c("GO:0006833", "GO:0009828", "GO:0071555", "GO:0009913", "GO:0048527"))
nutrients <- TermConversion(c("GO:0016036", "GO:0006995", "GO:0006817", "GO:0015706"))
metals <- TermConversion(c("GO:0006826", "GO:0000041", "GO:0046686", "GO:0010043", "GO:0006865"))
ethylene <- TermConversion(c("GO:0009723", "GO:0009693", "GO:0071369"))

# Function to create and save plots for each functional group
CreateSetOfPlots <- function(normcounts, golist, deresult, gogenes, savedir) {
  if (!dir.exists(savedir)) {dir.create(savedir)}
  for (i in 1:length(golist)) {
    goname <- gsub(":", "", golist[i])
    fname <- paste0(savedir, "/", goname, ".png")
    csv_name <- paste0(savedir, "/", goname, ".csv")
    pheatmap(normcounts[gogenes[[i]], ], display_numbers = FALSE, angle_col = 0,
             main = goname, filename = fname,
             fontsize = 5, number_color = "black", color = colorRampPalette(rev(brewer.pal(n = 7, name = "YlGnBu")))(100))
    write.csv(normcounts[gogenes[[i]], ], csv_name)
  }
}

# Normalized counts for ETR1 vs ETR1E
normalized_counts_ETR <- voom_limma_ETR$E
MeanE_ETR <- rowMeans(normalized_counts_ETR[, 1:5])
MeanR_ETR <- rowMeans(normalized_counts_ETR[, 6:10])
normalized_mean_ETR <- data.frame(MeanE_ETR, MeanR_ETR)
colnames(normalized_mean_ETR) <- c("Mean normalized expression ETR1", "Mean normalized expression ETR1E")
rownames(normalized_mean_ETR) <- rownames(normalized_counts_ETR)

# Create and save plots for ETR1 vs ETR1E
CreateSetOfPlots(normalized_mean_ETR, auxin, TMM_voom_counts_ETR, genes_GO_limma_ETR, "./plots/GO/ETR_X_WT/auxin")
CreateSetOfPlots(normalized_mean_ETR, root, TMM_voom_counts_ETR, genes_GO_limma_ETR, "./plots/GO/ETR_X_WT/root")
CreateSetOfPlots(normalized_mean_ETR, other, TMM_voom_counts_ETR, genes_GO_limma_ETR, "./plots/GO/ETR_X_WT/other")
CreateSetOfPlots(normalized_mean_ETR, nutrients, TMM_voom_counts_ETR, genes_GO_limma_ETR, "./plots/GO/ETR_X_WT/nutrients")
CreateSetOfPlots(normalized_mean_ETR, metals, TMM_voom_counts_ETR, genes_GO_limma_ETR, "./plots/GO/ETR_X_WT/metals")
CreateSetOfPlots(normalized_mean_ETR, ethylene, TMM_voom_counts_ETR, genes_GO_limma_ETR, "./plots/GO/ETR_X_WT/ethylene")

# Normalized counts for WT vs WTE
normalized_counts_WT <- voom_limma_WT$E
MeanE_WT <- rowMeans(normalized_counts_WT[, 1:5])
MeanR_WT <- rowMeans(normalized_counts_WT[, 6:10])
normalized_mean_WT <- data.frame(MeanE_WT, MeanR_WT)
colnames(normalized_mean_WT) <- c("Mean normalized expression WT", "Mean normalized expression WTE")
rownames(normalized_mean_WT) <- rownames(normalized_counts_WT)

# Create and save plots for WT vs WTE
CreateSetOfPlots(normalized_mean_WT, auxin, TMM_voom_counts_WT, genes_GO_limma_WT, "./plots/GO/WT_X_WTE/auxin")
CreateSetOfPlots(normalized_mean_WT, root, TMM_voom_counts_WT, genes_GO_limma_WT, "./plots/GO/WT_X_WTE/root")
CreateSetOfPlots(normalized_mean_WT, other, TMM_voom_counts_WT, genes_GO_limma_WT, "./plots/GO/WT_X_WTE/other")
CreateSetOfPlots(normalized_mean_WT, nutrients, TMM_voom_counts_WT, genes_GO_limma_WT, "./plots/GO/WT_X_WTE/nutrients")
CreateSetOfPlots(normalized_mean_WT, metals, TMM_voom_counts_WT, genes_GO_limma_WT, "./plots/GO/WT_X_WTE/metals")
CreateSetOfPlots(normalized_mean_WT, ethylene, TMM_voom_counts_WT, genes_GO_limma_WT, "./plots/GO/WT_X_WTE/ethylene")
