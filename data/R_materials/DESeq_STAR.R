# Manually tell RStudio Server where to find the pre-installed R Packages
.libPaths(c("/lustre/isaac24/proj/UTK0386/R"))

# Load up the packages
library(tidyverse)
library(tximport)
library(GenomicFeatures)
library(txdbmaker)
library(pheatmap)
library(DESeq2)
library(apeglm)
library(EnhancedVolcano)
library(ggpubr)

# create a workding directory
setwd(paste0("/lustre/isaac24/proj/UTK0386/analysis/", Sys.getenv("USER"), "/05_count"))


################################################################################
# load counts, metadata, and create deseq objects
################################################################################
# load the featureCounts table
fc <- read.table("combined.counts.txt", header = TRUE, row.names = 1, comment.char = "#")

# get the sample columns only (row names are genes)
count_matrix <- fc[, 6:ncol(fc)]

# clean up the column names to be more tidy
colnames(count_matrix) <- gsub("\\.Aligned.sortedByCoord.out.bam$", "", colnames(count_matrix))

# load up the metadata
samples <- read_csv("/lustre/isaac24/proj/UTK0386/data/star_data.csv")

# narrow down to those samples in our featurecounts matrix
samples <- samples[samples$sample_id %in% colnames(count_matrix), ]

# convert fields to factors
samples$treatment <- factor(samples$treatment)
samples$time      <- factor(samples$time)
samples$replicate <- factor(samples$replicate)


# make the order of samples match the metadata
count_matrix <- count_matrix[, samples$sample_id]

# build the deseq object (usually called dds) including a simple design
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = samples, design = ~ time)

# run deseq2 (this step does the actual normalization)
dds <- DESeq(dds)

# plot dispersion
plotDispEsts(dds)


################################################################################
# summarize results
################################################################################
# our simple model lets us make a very simple results object
res <- results(dds)

head(res)
summary(res)

# order the results by smallest padj, then largest absolute log change
res_ord <- res[order(res$padj, -abs(res$log2FoldChange)),]

# save the file
write.csv(res_ord, file="star_featurecounts_results.csv")

# get a list of the groupings (to use in contrasts if > 2 groups are in your design)
resultsNames(dds)


# perform variance stabilization transformation (move data to a scale less influenced by highly expressed genes)
vsd <- vst(dds)

# plot a PCA using the vst object, color groups by time (0 vs. 3 hours)
plotPCA(vsd, intgroup = c("time"))


# create a contrast with alpha cutoff (first list item is condition from samples object)
res_time <- results(dds, alpha = 0.05, contrast = c("time", "3h", "0h"))

# make a MA plot of differential expression
plotMA(res_time, ylim=c(-12,12))

# now create a contrast in the opposite direction
res_time <- results(dds, alpha = 0.05, contrast = c("time", "0h", "3h"))

# does the MA plot look identical?
plotMA(res_time, ylim=c(-12,12))

# run this, then click on a point on the plot, be sure to hit "Finish" in the upper right of the plot
idx <- identify(res$baseMean, res$log2FoldChange)

# now see the gene name of that point
rownames(res)[idx]

# manually filter the file down to desired padj level...
res_sig <- as.data.frame(res[ which(res$padj < 0.05),])

# ...and sort based on log2fold values ()
res_sig <- res_sig[order(res_sig$padj, -abs(res_sig$log2FoldChange)),]

# write the significant results to file
write.csv(res_sig, file="star_featurecounts_sig_results.csv")

# make a detailed volcano plot
EnhancedVolcano(res, x = 'log2FoldChange',
                lab = row.names(res),
                pCutoff = 1e-100,
                FCcutoff = 2,
                y = 'padj',)

# create a plot for a single gene using the normalized counts
plotCounts(dds, gene="gene-AT4G25480", intgroup="time")

# create a heatmap using the vst object from above (top 20 genes based on log2foldChange)
mat <- assay(vsd)[ head(order(-abs(res$log2FoldChange)), 20), ]
mat <- mat - rowMeans(mat) 
mat <- data.frame(mat)

pheatmap(mat, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE)
