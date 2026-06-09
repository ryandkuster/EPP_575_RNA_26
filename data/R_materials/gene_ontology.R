# Manually tell RStudio Server where to find the pre-installed R Packages
.libPaths(c("/lustre/isaac24/proj/UTK0386/R"))

# Load up the packages
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(pathview)
library(DOSE)
library(gage)
library(pathview)
library("org.At.tair.db", character.only = TRUE)

# create a workding directory
setwd(paste0("/lustre/isaac24/proj/UTK0386/analysis/", Sys.getenv("USER"), "/05_count"))

# confirm you've loaded up the Arabidopsis Database for GO
keytypes(org.At.tair.db)

############################################################
# import DGE data
############################################################
# open the full results file from previous lab
df = read.csv("star_featurecounts_results.csv", header=TRUE)

# rename the former row variable (called X by default) to "gene"
names(df)[names(df) == "X"] <- "gene"

# make a new gene names column ("Name") by stripping off the "gene-" prefix using gsub function
df$Name <- sub("^gene-", "", df$gene)

# manually make a list of the log2fold change values only
l2f_gene_list <- df$log2FoldChange
length(l2f_gene_list)

# make this a named list by adding the Name column (in same order as list we just made)
names(l2f_gene_list) <- df$Name

# remove any missing genes from our list
gene_list<-na.omit(l2f_gene_list)
length(l2f_gene_list)

# sort the list from greatest to least log2fold change values (this is critical!)
gene_list = sort(gene_list, decreasing = TRUE)


############################################################
# GO approach
############################################################
# we want to see if GO terms found in our genes tend to cluster near the top or bottom of our list (recall we just sorted it)
# ont = "ALL" means we want biological, molecular, and cellular terms
# min and max GSS (gene set size) determine how small and large clusters of terms to 
# the minGSS value here is rather low, but we have a stringent pvalueCutoff and use multiple testing correction (pAdjustMethod)

gsea <- gseGO(geneList=gene_list, 
              ont ="ALL", 
              keyType = "TAIR",  
              minGSSize = 10, 
              maxGSSize = 1000, 
              pvalueCutoff = 0.01, 
              verbose = TRUE, 
              OrgDb = org.At.tair.db,
              pAdjustMethod = "BH")

dotplot(gsea, showCategory=10, split=".sign") + facet_grid(.~.sign)

# optional, but useful, plots
ridgeplot(gsea) + labs(x = "enrichment distribution")

# for (i in 1:5) {
#   p <- gseaplot(gsea, by = "all", title = gsea$Description[i], geneSetID = i)
#   print(p)
# }
# 
# # grab a few terms and pull any relevant literature to show trends in findings
# terms <- gsea$Description[1:4]
# pmcplot(terms, 2010:2026, proportion=FALSE)

############################################################
# gage approach
############################################################
# for many organisms, we won't be lucky enough to have a TAIR dedicated GO database
# the gage database works for some species (supply the shorthand name below)

# kg.ath <- kegg.gsets(species="ath")
# kegg_df = kg.ath[["kg.sets"]]
# unique_entries <- unique(unlist(kegg_df))
# binary_matrix <- sapply(kegg_df, function(x) as.integer(unique_entries %in% x))
# result_df <- as.data.frame(binary_matrix, row.names = unique_entries)
# colnames(result_df) <- names(kegg_df)
# result_df <- tibble::rownames_to_column(result_df, var = "Name")
# merged_df <- merge(df, result_df, by="Name", all.x=TRUE)


############################################################
# KEGG approach
############################################################
kegg_species <- clusterProfiler:::kegg_species_data()
ath <- search_kegg_organism("ath", by='kegg_code')

# over-representation test (only uses gene names)
res_sig = subset(df, abs(df$log2FoldChange) > 2)
res_sig = subset(res_sig, res_sig$padj < 0.05)
gene <- res_sig$Name
length(gene)
"AT4G25480" %in% gene

kk <- enrichKEGG(gene,
                 organism="ath",
                 pvalueCutoff=0.05,
                 pAdjustMethod="BH")
head(kk)

# let's look at the top two hits
# green means the gene is in "ath" set, and red means it is in our enriched set
browseKEGG(kk, 'ath04075')
browseKEGG(kk, 'ath04626')

# gseKEGG gene set enrichment analysis (gsea) (takes in log2FC with gene names)
kk2 <- gseKEGG(geneList      = gene_list,
               organism      = "ath",
               pvalueCutoff  = 0.05,
               verbose       = FALSE,
               pAdjustMethod = "fdr")

summary(as.data.frame(kk2))
dotplot(kk2, showCategory=10, split=".sign") + facet_grid(.~.sign)


############################################################
# work-around for kegg ath database not importing by default
############################################################
korg <- cbind("ktax.id" = "T00041", "tax.id" = "3702", "kegg.code" = "ath",
              "scientific.name" = "Arabidopsis thaliana", "common.name" = "thale cress",
              "entrez.gnodes" = "1", "kegg.geneid" = "816394", "ncbi.geneid" = "816394",
              "ncbi.proteinid" = "NP_001325249.1 ", "uniprot" = "A0A178VXX1")

pathview(gene.data  = gene_list,
         pathway.id = "ath04075",
         species    = "ath")
