##LOADING IN DATA
## load SN503 meta data
## NOTE: you will need to change the path so that it matches the mounted directory path on your laptop/computer
filepath  <- "Z:/PD_GSC_pstat/Results/scRNA_seq/10X_chromium/SN503_pstat_2019.04.29/R_code/"
filename  <- "SN503_metaData.RData"

load(paste0(filepath, filename))

## load SN520 meta data
filepath  <- "Z:/PD_GSC_pstat/Results/scRNA_seq/10X_chromium/SN520_pstat_2019.01.30/R_code/"
filename  <- "SN520_metaData.RData"

load(paste0(filepath, filename))

## from SN520_pstat_2019.01.30_normalization.Rmd
load("/mnt/omics4tb2/jpark/Projects/PD_GSC_pstat/Results/scRNA_seq/10X_chromium/SN520_pstat_2019.01.30/R_code/SN520_seuratObjs_scran_norm_2022.01.24.RData")

filepath  <- "Z:/PD_GSC_pstat/Results/scRNA_seq/10X_chromium/SN520_pstat_2019.01.30/R_code/"
filename  <- "SN520_SN503_miner3_dataPrep_objs_2022_04_28.RData"

load(paste0(filepath, filename))
# common_genes
# genes2use 
# ensembl
# ensemblID.df
# sn503.data.scaled.4miner
# sn503.data.scaled.4miner.zscore
# sn520.data.scaled.4miner
# sn520.data.scaled.4miner.zscore

dim(sn503.data.scaled.4miner)
dim(sn503.data.scaled.4miner.zscore)

sn503.data.scaled.4miner[1:5, 1:5]

## convert ensemblID to gene symbol
gene_symb  <- ensemblID.df$hgnc_symbol[match(rownames(sn503.data.scaled.4miner), ensemblID.df$ensembl_gene_id)]

## assigne gene symbols to rownames of data matrices
rownames(sn503.data.scaled.4miner)  <- gene_symb
rownames(sn520.data.scaled.4miner)  <- gene_symb

rownames(sn503.data.scaled.4miner.zscore)  <- gene_symb
rownames(sn520.data.scaled.4miner.zscore)  <- gene_symb

##INSTALL PACKAGES
install.packages("cluster")
BiocManager::install("SingleCellExperiment")
BiocManager::install("scater")
BiocManager::install("scran")
install.packages("hdf5r")
BiocManager::install("limma")
BiocManager::install("GO.db")
BiocManager::install("org.Hs.eg.db")

##RUN PACKAGES##
library(Seurat)
library(SingleCellExperiment)
library(scater)
library(scran)
library(cluster)
library(PCAtools)
library(hdf5r)
library(GO.db)
library(limma)
library(org.Hs.eg.db)

##CREATING NEW DMSO DATASETS##
#SN503 DMSO DATA#
#creating new sn503 dataset with only DMSO data. Subset the dmso values from the 503 metadata. 
sn503_dmso_metadata <- subset(sn503.meta.data, condition_pstat_d4_union == "DMSO_d4" | condition_pstat_d4_union == "DMSO_d2" | condition_pstat_d4_union == "DMSO_d3")

#get dmso row names as an array
sn503_dmso_row_names <- c(row.names(sn503_dmso_metadata))

#get only the sn503 rna seq values that are dmso treated
sn503_dmso_only_subset <- sn503.data.scaled.4miner.zscore[, colnames(sn503.data.scaled.4miner.zscore) %in% sn503_dmso_row_names]

#SN520 DMSO DATA#
#do the same thing for sn520 - subset the dmso values from sn520 metadata
sn520_dmso_metadata <- subset(sn520.meta.data, condition_pstat_d4_union == "DMSO_d4" | condition_pstat_d4_union == "DMSO_d2" | condition_pstat_d4_union == "DMSO_d3")

#get dmso row names as an array
sn520_dmso_row_names <- c(row.names(sn520_dmso_metadata))

#get only the sn520 rna seq values that are dmso treated
sn520_dmso_only_subset <- sn520.data.scaled.4miner.zscore[, colnames(sn520.data.scaled.4miner.zscore) %in% sn520_dmso_row_names]

##SN520 DMSO TREATED ONLY: DIMENSIONALITY REDUCTION + CLUSTERING##
#turn data into seurat object
seu <- CreateSeuratObject(counts=sn520_dmso_only_subset, assay="RNA", meta.data=sn520_dmso_metadata) 

#turn data into single cell experiment object
sce <- as.SingleCellExperiment(seu)
sce

#run pca and reduce dimensionality
#set seed
set.seed(1234)
sce <- runPCA(sce)
percent.var <- attr(reducedDim(sce),"percentVar")
#find optimal number of PCs using an elbow plot, plot elbow plot
chosen.elbow <- findElbowPoint(percent.var)
chosen.elbow
plot(percent.var, xlab="PC", ylab="Variance explained (%)")
abline(v=chosen.elbow, col="red")
#reduce number of PCs and plot data based on day
reducedDim(sce, "PCA") <- reducedDim(sce, "PCA")[,1:50]
plotPCA(sce, colour_by = "condition_pstat_d4_union")

#GRAPH BASED CLUSTERING#
#overlay clusters on PCA graph, k value can be changed
g <- buildSNNGraph(sce, k=20, use.dimred = "PCA")
clust <- igraph::cluster_walktrap(g)$membership
table(clust)
colLabels(sce) <- factor(clust)
plotReducedDim(sce, "PCA", colour_by="label") #label = cluster

##FIND DIFFERENTIAL EXPRESSION PATTERNS ACROSS CLUSTERS##

marker.info <- scoreMarkers(sce, colLabels(sce))
marker.info
#take marker information for x cluster
chosen <- marker.info[["4"]]
#order the clusters based on differential expression patterns
ordered <- chosen[order(chosen$mean.AUC, decreasing=TRUE),]
head(ordered[,1:4])
#plot clusters based on differential expression patterns
plotExpression(sce, features=head(rownames(ordered)), 
               x="label", colour_by="label")
#prints the area underneath the curve, showing us upregulation/downregulation patterns (1 is up, 0 is none, -1 is down)
auc.only <- chosen[,grepl("AUC", colnames(chosen))]
auc.only[order(auc.only$mean.AUC,decreasing=TRUE),]

##SN503 DMSO TREATED ONLY: DIMENSIONALITY REDUCTION + CLUSTERING##
#turn data into seurat object
seu <- CreateSeuratObject(counts=sn503_dmso_only_subset, assay="RNA", meta.data=sn503_dmso_metadata) 

#turn data into single cell experiment object
sce <- as.SingleCellExperiment(seu)
sce

#set seed
set.seed(1234)

#run pca and reduce dimensionality
sce <- runPCA(sce)
percent.var <- attr(reducedDim(sce),"percentVar")
chosen.elbow <- findElbowPoint(percent.var)
chosen.elbow
plot(percent.var, xlab="PC", ylab="Variance explained (%)")
abline(v=chosen.elbow, col="red")
reducedDim(sce, "PCA") <- reducedDim(sce, "PCA")[,1:2]
plotPCA(sce, colour_by = "condition_pstat_d4_union")

#GRAPH BASED CLUSTERING#
g <- buildSNNGraph(sce, k=10, use.dimred = "PCA")
clust <- igraph::cluster_walktrap(g)$membership
table(clust)
colLabels(sce) <- factor(clust)
plotReducedDim(sce, "PCA", colour_by="label", theme_size = 5)

##FIND DIFFERENTIAL EXPRESSION PATTERNS ACROSS CLUSTERS##
marker.info <- scoreMarkers(sce, colLabels(sce), lfc = 1)
marker.info
#take marker information for x cluster
chosen <- marker.info[["1"]]
#order the clusters based on differential expression patterns
ordered <- order(chosen$mean.AUC, decreasing=TRUE)
chosen[ordered,1:4]#prints out 6 potential marker genes w/ highest mean AUC values (most upregulated)
#plot clusters based on differential expression patterns
plotExpression(sce, features=head(rownames(chosen[ordered,1:4])), 
               x="label", colour_by="label", theme_size = 5)
#prints the area underneath the curve, showing us upregulation/downregulation patterns (1 is up, 0 is none, -1 is down)
auc.only <- chosen[,grepl("AUC", colnames(chosen))]
auc.only[order(auc.only$mean.AUC,decreasing=TRUE)[1:6],]

##GENE ENRICHMENT ANALYSIS##
#finds GO terms (actions of gene products/pathways) overrepresented in our dataset
lfc.markers <- scoreMarkers(sce, lfc = 1)
#only choose the gene markers in x cluster
cur.markers <- lfc.markers[["4"]]
#order the genes based on mean AUC and then subset the first 100
is.de <- order(cur.markers$mean.AUC, decreasing=TRUE)[1:100]
cur.markers[is.de,1:4]

#map genes to pathways
entrez.ids <- mapIds(org.Hs.eg.db, keys=rownames(cur.markers), 
                     column="ENTREZID", keytype="SYMBOL")
go.out <- goana(unique(entrez.ids[is.de]), species="Hs", 
                universe=unique(entrez.ids))
#subset reference genes (only looking at biological pathways that are not overly general)
go.out <- go.out[order(go.out$P.DE),]
go.useful <- go.out[go.out$Ont=="BP" & go.out$N <= 200,]
head(go.useful[,c(1,3,4)], 30)