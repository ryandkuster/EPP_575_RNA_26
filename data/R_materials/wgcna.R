# Manually tell RStudio Server where to find the pre-installed R Packages
.libPaths(c("/lustre/isaac24/proj/UTK0386/R"))

library(tidyverse)
library(magrittr)
library(WGCNA) 
library(tibble)
library(readr)
library(DESeq2)
library(tximport)
library(GenomicFeatures)
library(genefilter)
library(R.utils)
library(txdbmaker)

allowWGCNAThreads()

setwd(paste0("/lustre/isaac24/proj/UTK0386/analysis/", Sys.getenv("USER"), "/05_count"))

# give the full path to the "salmon_results" (it should end in "salmon_results")
salmon_dir <- "/lustre/isaac24/proj/UTK0386/data/EPP_575_RNA_26/data/salmon_output/"

# load the gff3 file, then create a transcript database/dataframe for use with deseq
download.file("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/735/GCA_000001735.2_TAIR10.1/GCA_000001735.2_TAIR10.1_genomic.gff.gz", destfile = "GCA_000001735.2_TAIR10.1_genomic.gff.gz", method = "wget")
gunzip("GCA_000001735.2_TAIR10.1_genomic.gff.gz", remove = TRUE)

txdb <- makeTxDbFromGFF("GCA_000001735.2_TAIR10.1_genomic.gff")
keytypes(txdb)
k <- keys(txdb, keytype = "CDSNAME")
str(k)

txdf = AnnotationDbi::select(txdb, k, "GENEID", "CDSNAME")

samples <- read_csv(paste0(salmon_dir, "salmon_data.csv"))
Qfiles <- file.path(salmon_dir, samples$quant_file)


# this step imports the count data from salmon
txi <- tximport(files = Qfiles, type = "salmon", txOut = TRUE)
head(txi$counts)
colnames(txi$counts) <- samples$sample_id

# convert fields to factors
samples$treatment = factor(samples$treatment)
samples$time = factor(samples$time)
samples$replicate = factor(samples$replicate)

# now we convert the txi object into a deseq-formatted object
dds <- DESeqDataSetFromTximport(txi = txi, colData = samples, design = ~ treatment + time + treatment:time)
dds <- DESeq(dds, test="LRT", reduced = ~ treatment + time)

# if a model fails to converge (usually due to many zero count genes) remake model dropping those genes
dds <- dds[which(mcols(dds)$fullBetaConv),]

# transform data
vsd <- varianceStabilizingTransformation(dds)

# get vst data as a matrix for wgcna
wpn_vsd <- getVarianceStabilizedData(dds)
rv_wpn <- rowVars(wpn_vsd)
summary(rv_wpn)

# limit the count by quantiles
sapply(c(0.50, 0.75, 0.90), function(q) sum(rv_wpn > quantile(rv_wpn, q)))

# (this will remove MANY genes, try a few settings and see how it affects your dataset)
# .95 will keep only the top 10% most variable genes
q_cutoff = .9

quantile_wpn <- quantile( rowVars(wpn_vsd), q_cutoff)  # <= changed to 90 quantile to reduce dataset
expr_normalized <- wpn_vsd[ rv_wpn > quantile_wpn, ]
dim(expr_normalized)

# make data "long" for ggplot
expr_normalized_df <- data.frame(expr_normalized) %>%
  mutate(
    Gene_id = row.names(expr_normalized)
  ) %>%
  pivot_longer(-Gene_id)

# make a violin plot of normalized expression values
expr_normalized_df %>% ggplot(., aes(x = name, y = value)) +
  geom_violin() +
  geom_point() +
  theme_bw() +
  theme(
    axis.text.x = element_text( angle = 90)
  ) +
  ylim(0, NA) +
  labs(
    title = "Normalized and 90 quantile Expression",
    x = "treatment",
    y = "normalized expression"
  )

input_mat = t(expr_normalized)
input_mat[1:5,1:6]

# visualize the matrix of normalized expression values as hierarchical clustering of samples
sampleTree = hclust(dist(input_mat), method = "average")
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# now we will visualize our topology to determine a target cutoff
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# WGCNA converts a correlation matrix into a weighted network by raising
# correlations to a power β. This step tests a range of powers to find the one
# where the network best approximates scale-free topology — a property of most
# real biological networks where a few "hub" genes have many connections and
# most genes have few.

sft = pickSoftThreshold(
  input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)

# make a plot; this will help determine the threshold to use
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")


# STOP!
# visually inspect the plot to see if there is a point where we see an "elbow"
# this elbow tells us where we get diminishing returns from increasing the soft
# power threshold (pick a value just at or below the .90 R squared line)
# these are simulated networks based on different power thresholds
picked_power = 7
temp_cor <- cor       
cor <- WGCNA::cor

# this step takes a minute...
netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

cor <- temp_cor

# convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)

# plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

# wow! so here we see the dendogram of the topological overlap matrix (TOM)
# each leaf is a gene, and they have been clustered based on the network of all
# gene correlations (clusters are impacted by soft-threshold value above) 

netwk$colors[netwk$blockGenes[[1]]]
table(netwk$colors)

module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

module_df[1:5,]

# save the network
write_delim(module_df,
            file = "gene_modules.txt",
            delim = "\t")

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$treatment = row.names(MEs0)

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")

# pick out a few modules of interest here
modules_of_interest = c("magenta")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id

subexpr = expr_normalized[submod$gene_id,]

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "treatment",
       y = "normalized expression")


genes_of_interest = module_df %>%
  subset(colors %in% modules_of_interest)

expr_of_interest = expr_normalized[genes_of_interest$gene_id,]

# TOMs are "topical overlay matrices
# how similar are gene A and gene B's entire connection patterns across all other genes?
TOM = TOMsimilarityFromExpr(t(expr_of_interest),
                            power = picked_power)

row.names(TOM) = row.names(expr_of_interest)
colnames(TOM) = row.names(expr_of_interest)

edge_list = data.frame(TOM) %>%
  mutate(
    gene1 = row.names(.)
  ) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, correlation = value) %>%
  unique() %>%
  subset(!(gene1==gene2)) %>%
  mutate(
    module1 = module_df[gene1,]$colors,
    module2 = module_df[gene2,]$colors
  )

# inspect the edge list; it contains all gene to gene correlations
head(edge_list)

# save edge list
write_delim(edge_list,
            file = "edgelist.tsv",
            delim = "\t")
