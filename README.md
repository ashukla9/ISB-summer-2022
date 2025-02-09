# Analyzing Glioblastoma Data Using PCA And K-Means Clustering

I used this Bioconductor tutorial (http://bioconductor.org/books/3.15/OSCA.basic/cell-type-annotation.html#assigning-cluster-labels-from-markers) and this clustering tutorial (https://www.gaininginsight.org/cmwg/2020/coronavirus/scrna-seq-analysis) to learn more about analyzing RNA-seq data. 

First, I subset the sn503.data.scaled.4miner.zscore and sn520.data.scaled.4miner.zscore data to only include the cells treated with DMSO. Because my RStudio wouldn’t let me merge the two datasets together (I kept getting a memory allocation error), I analyzed the two datasets separately. I used the elbow method to find the optimal number of PCs, then ran PCA on both “DMSO only” datasets and labeled the cells based on day. I didn’t see a large difference between the day 2, day 3, and day 4 cells. 

For this reason, I focused on just looking at the sn503 and sn520 DMSO data across all time points rather than breaking up the datasets further into individual days. For future analysis, it might be interesting to look at clustering within each day to see whether the same genes/pathways are highly differentially expressed across all three days.

I then ran graph-based clustering on top of my PCA graph for both the sn503 and sn520 “DMSO only” datasets. I saw much more variation within the sn503 data (as defined by number of clusters) than within the sn520 data, which might speak to the heterogeneity of the non-responder cells. 

I used a k of 10 for the sn503 data and a k of 20 for the sn520 data  (a k of 10 for the sn520 data created 6 clusters, but not all of them were well-defined). These “k”s (the number of nearest neighbors) were determined through an ad hoc, “what makes my clusters the most well-defined without overfitting the data” process, so a future step would be to use mathematical processes to find the optimal k for each dataset. I chose to focus on graph-based clustering, but it might be interesting to run k-means clustering and compare the results to my graph-based clusters, as well as run subclustering to detect subtle variations within major clusters.

I then looked within each cluster to determine the genes that were upregulated in that cluster compared to other clusters (essentially, genes that are highly differentially expressed across clusters). I measured differential expression based on the mean AUC, or area under the curve, which “represents the probability that a randomly chosen observation from our cluster of interest is greater than a randomly chosen observation from the other cluster” (as per the Bioconductor tutorial). A value closer to 1 means the gene is upregulated, and a value closer to 0 means the gene is downregulated. It is also possible to measure up/downregulation by looking at the minimum or median AUC—I chose the mean as it is a middle-of-the-road default for this analysis. I applied a log-fold threshold to reduce the likelihood of seeing genes with large effect sizes but low differential expression across clusters.

Within each cluster, I found the 6 genes that had the greatest mean AUC and plotted their expression against other clusters. In several of the clusters, the genes with the highest mean AUCs had a mean AUC value in the 0.6 - 0.7 range, which doesn’t seem very significant in comparison to mean AUC values from other clusters. A future step might be to analyze clusters with genes with lower mean AUCs to see if there is significant downregulation of genes instead of upregulation. 

Finally, per cluster, I ran enrichment analysis on 100 genes with the highest mean AUC. I used Bioconductor’s genome-wide human annotation package to determine the 30 most active pathways within each cluster. sn503 had a diverse array of pathways (e.g. DNA/RNA regulation and repair, ion response, cell differentiation, nervous system regulation), while sn520’s pathways mostly related to the cell cycle/division.  

sn503 Data - Summary of Pathways:
Cluster #1: cardiovascular and ion-related
Cluster #2: ribosomal activity, RNA transportation, chromatin/nucleosome assembly
Cluster #3: O-glycan processing, glycosylation
Cluster #4: neuron differentiation and nervous system development
Cluster #5: membrane proteins, ion/cation response
Cluster #6: ribosomal assembly, DNA and RNA regulation
Cluster #7: myeloid dendritic cells, transmembrane transportation
Cluster #8: intracellular transportation, metabolic process, muscle contractions
Cluster #9: DNA repair and regulation, extracellular organization
Cluster #10: protein localization regulation, regulation of protein catabolic processes
Cluster #11: lipid biosynthetic processes, Schwann cells
Cluster #12: cell growth and development, negative regulation of nervous system
Cluster #13: axonogenesis, axon and nervous system development regulation
Cluster #14: ion stress response/detoxification
Cluster #15: ribosomal assembly, endothelial/epithelial cell fate commitment
Cluster #16: import across plasma membrane/blood-brain barrier, D-aspertate import
Cluster #17: myeloid cell differentiation, epithelial cell development

sn520 Data - Summary of Pathways:
Cluster #1: cell division
Cluster #2: neuron/axon regeneration and death, glial cell development, interleukin-1 beta production, ion transport/response
Cluster #3: protein catabolic and biosynthetic processes, telomere maintenance
Cluster #4: histone phosphorylation, cell division

As it doesn’t make sense that, post-processing, sn520 would have so many pathways related to cell division, a future step would be to extract relevant genes from each pathway to cross-validate the results, as well as statistically analyze the significance of the pathways determined.
