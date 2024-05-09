# Load necessary libraries
library(limma)
library(DESeq2)
library(edgeR)
library(topGO)
library(dplyr)
library(xlsx)
library(org.At.tair.db)
library(jsonlite)
library(glmnetUtils)
library(ggplot2)

# Function to calculate quartiles
quart <- function(x) {
  x <- sort(x)
  n <- length(x)
  m <- (n + 1) / 2
  if (floor(m) != m) {
    l <- m - 1 / 2
    u <- m + 1 / 2
  } else {
    l <- m - 1
    u <- m + 1
  }
  c(Q1 = median(x[1:l]), Q3 = median(x[u:n]))
}

#---- Data Preparation

# Load count data for ETR1 and WT conditions
ETR_counts <- read.csv("./counts/ETR_ETRE.csv")
WT_counts <- read.csv("./counts/WT_WTE.csv")

# Set row names and remove first column
rownames(ETR_counts) <- ETR_counts$X
rownames(WT_counts) <- WT_counts$X
ETR_counts <- ETR_counts[, -1]
WT_counts <- WT_counts[, -1]

# Extract sample identifiers
samples_ETR <- substr(colnames(ETR_counts), 9, 13)
samples_WT <- substr(colnames(WT_counts), 9, 11)

# Create design matrices for limma analysis
design_limma_ETR <- data.frame(ETR1s = ifelse(samples_ETR != "ETR1E", 1, 0),
                               ETRE1s = ifelse(samples_ETR == "ETR1E", 1, 0))
design_limma_WT <- data.frame(WTs = ifelse(samples_WT != "WTE", 1, 0),
                              WTEs = ifelse(samples_WT == "WTE", 1, 0))

#---- Limma Analysis for ETR1 vs ETR1E

# Prepare DGEList object and normalize
DGElimma_ETR <- DGEList(counts = ETR_counts)
DGElimma_ETR <- calcNormFactors(DGElimma_ETR)
voom_limma_ETR <- voom(DGElimma_ETR, design_limma_ETR, plot = TRUE)

# Fit linear model and apply contrasts
contrasts_ETR <- makeContrasts(ETR1EsvsETR1s = ETRE1s - ETR1s,
                               levels = design_limma_ETR)
voom_fit_ETR <- lmFit(voom_limma_ETR, contrast = contrasts_ETR)
cf_limma_ETR <- contrasts.fit(voom_fit_ETR, contrasts_ETR)
eBayes_limma_ETR <- eBayes(cf_limma_ETR, proportion = 0.01)

# Extract significant genes and filter based on logFC and average expression
TMM_voom_counts_ETR <- topTable(eBayes_limma_ETR, number = Inf, adjust.method = "BH", sort.by = "none", confint = TRUE)
significant_limma_ETR <- TMM_voom_counts_ETR[TMM_voom_counts_ETR$adj.P.Val <= 0.05, ]
LogFC_boolvec <- abs(significant_limma_ETR$logFC) >= 1
significant_limma_ETR <- significant_limma_ETR[LogFC_boolvec, ]
AveExp_Q1 <- quart(significant_limma_ETR$AveExpr)[1]
AveExp_boolvec <- significant_limma_ETR$AveExpr >= AveExp_Q1
significant_limma_ETR <- significant_limma_ETR[AveExp_boolvec, ]
write.csv(file = "significant_ETR1.csv", x = significant_limma_ETR)

# Plot MA plot for ETR1 vs ETR1E
pchs <- rep('.', nrow(TMM_voom_counts_ETR))
colors <- rep('black', nrow(TMM_voom_counts_ETR))
names(pchs) <- names(colors) <- rownames(TMM_voom_counts_ETR)
selected_transcripts <- rownames(significant_limma_ETR)
pchs[selected_transcripts] <- "x"
colors[selected_transcripts] <- "purple"

plot(TMM_voom_counts_ETR$AveExpr, TMM_voom_counts_ETR$logFC,
     xlab = 'Average Expression', ylab = "logFC",
     main = "Limma MA plot", pch = pchs, col = colors)

# Save normalized counts for ETR1 vs ETR1E
write.csv(voom_limma_ETR$E, "normalized_counts_ETR1.csv")

# Identify low copy number root genes and save
potential <- TMM_voom_counts_ETR[TMM_voom_counts_ETR$adj.P.Val <= 0.05, ]
LogFC_boolvec <- abs(potential$logFC) >= 1
potential <- potential[LogFC_boolvec, ]
AveExp_Q1 <- quart(potential$AveExpr)[1]
AveExp_boolvec <- potential$AveExpr < AveExp_Q1
potential <- potential[AveExp_boolvec, ]
write.csv(file = "ETR_significant_low_copy_root_limma.csv", x = potential)

#---- Limma Analysis for WT vs WTE

# Prepare DGEList object and normalize
DGElimma_WT <- DGEList(counts = WT_counts)
DGElimma_WT <- calcNormFactors(DGElimma_WT)
voom_limma_WT <- voom(DGElimma_WT, design_limma_WT, plot = TRUE)

# Fit linear model and apply contrasts
contrasts_WT <- makeContrasts(WTsvsWTEs = WTEs - WTs, levels = design_limma_WT)
voom_fit_WT <- lmFit(voom_limma_WT, contrast = contrasts_WT)
cf_limma_WT <- contrasts.fit(voom_fit_WT, contrasts_WT)
eBayes_limma_WT <- eBayes(cf_limma_WT, proportion = 0.01)

# Extract significant genes and filter based on logFC and average expression
TMM_voom_counts_WT <- topTable(eBayes_limma_WT, number = Inf, adjust.method = "BH", sort.by = "none", confint = TRUE)
significant_limma_WT <- TMM_voom_counts_WT[TMM_voom_counts_WT$adj.P.Val <= 0.05, ]
LogFC_boolvec <- abs(significant_limma_WT$logFC) >= 1
significant_limma_WT <- significant_limma_WT[LogFC_boolvec, ]
AveExp_Q1 <- quart(significant_limma_WT$AveExpr)[1]
AveExp_boolvec <- significant_limma_WT$AveExpr >= AveExp_Q1
significant_limma_WT <- significant_limma_WT[AveExp_boolvec, ]
write.csv(file = "significant_WT.csv", x = significant_limma_WT)

# Plot MA plot for WT vs WTE
pchs <- rep('.', nrow(TMM_voom_counts_WT))
colors <- rep('black', nrow(TMM_voom_counts_WT))
names(pchs) <- names(colors) <- rownames(TMM_voom_counts_WT)
selected_transcripts <- rownames(significant_limma_WT)
pchs[selected_transcripts] <- "x"
colors[selected_transcripts] <- "purple"

plot(TMM_voom_counts_WT$AveExpr, TMM_voom_counts_WT$logFC,
     xlab = 'Average Expression', ylab = "logFC",
     main = "Limma MA plot", pch = pchs, col = colors)

# Save normalized counts for WT vs WTE
write.csv(voom_limma_WT$E, "normalized_counts_WT.csv")

#---- DESeq2 Analysis for ETR1 vs ETR1E

# Prepare metadata and DESeq2 dataset
coldata_ETR <- data.frame(sample = colnames(ETR_counts),
                          condition = ifelse(samples_ETR != "ETR1E", "ETR1", "ETR1E"))
coldata_ETR$condition <- as.factor(coldata_ETR$condition)

dds_ETR <- DESeqDataSetFromMatrix(countData = ETR_counts,
                                  colData = coldata_ETR,
                                  design = ~condition)

# Run DESeq2 analysis
dds_ETR <- DESeq(dds_ETR)
results_ETR <- results(dds_ETR)
results_ETR <- results_ETR[!is.na(results_ETR$padj), ]

# Extract significant genes
significant_ETR_DESeq2 <- results_ETR[results_ETR$padj <= 0.05, ]

# Plot MA plot for DESeq2 results
pchs <- rep('.', nrow(ETR_counts))
colors <- rep('black', nrow(ETR_counts))
names(pchs) <- names(colors) <- rownames(ETR_counts)
selected_transcripts <- rownames(significant_ETR_DESeq2)
pchs[selected_transcripts] <- "x"
colors[selected_transcripts] <- "purple"

plot(results_ETR$baseMean, results_ETR$log2FoldChange,
     xlab = 'Average Expression', ylab = "logFC",
     main = "DESeq2 MA plot", pch = pchs, col = colors,
     xlim = c(0, 2000))

#---- DESeq2 Analysis for WT vs WTE

# Prepare metadata and DESeq2 dataset
coldata_WT <- data.frame(sample = colnames(WT_counts),
                         condition = ifelse(samples_WT != "WTE", "WT", "WTE"))
coldata_WT$condition <- as.factor(coldata_WT$condition)

dds_WT <- DESeqDataSetFromMatrix(countData = WT_counts,
                                 colData = coldata_WT,
                                 design = ~condition)

# Run DESeq2 analysis
dds_WT <- DESeq(dds_WT)
results_WT <- results(dds_WT)
results_WT <- results_WT[!is.na(results_WT$padj), ]

# Extract significant genes
significant_WT_DESeq2 <- results_WT[results_WT$padj <= 0.05, ]

# Plot MA plot for DESeq2 results
pchs <- rep('.', nrow(WT_counts))
colors <- rep('black', nrow(WT_counts))
names(pchs) <- names(colors) <- rownames(WT_counts)
selected_transcripts <- rownames(significant_WT_DESeq2)
pchs[selected_transcripts] <- "x"
colors[selected_transcripts] <- "purple"

plot(results_WT$baseMean, results_WT$log2FoldChange,
     xlab = 'Average Expression', ylab = "logFC",
     main = "DESeq2 MA plot", pch = pchs, col = colors,
     xlim = c(0, 2000))
