setwd("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_integration_dna_rna/scrna_vs_scdna")

source("complete_workflow_functions.R")

###############
## scDNA data

cellids <- list.files("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/sample1.genemodel/tmp.parallel",pattern = "intcellids_",full.names = T)

ListCells <- lapply(1:length(cellids), loadCells, cellids)

genes <- unique(unlist( lapply(1:length(ListCells), getUniqueGenes, ListCells) ))

DTgenes <- data.table(gene.id = genes,stringsAsFactors = F)
rownames(DTgenes) <- DTgenes$gene.id

columns <- mclapply(1:length(ListCells), fillGenesCellsMatrix, ListCells, DTgenes, genes, mc.cores = 4)

scDNAmat <- do.call(cbind, columns)
rownames(scDNAmat) <- genes
colnames(scDNAmat) <- gsub(basename(cellids),pattern = "intcellids_|-genemodel.sort.bed",replacement = "")

rm(ListCells)

# clean dna data
load("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/check_gained_genes/masterDFgenemodel_diploid_annotation_sample1.Rdata")

keep <- as.character(master$id[which(master$flag10x_is.noisy == 0 & master$loss < .5)])

scDNAmat <- scDNAmat[, keep ]

scDNAmat <- scDNAmat[rowSums(is.na(scDNAmat)) == 0, ]

dim(scDNAmat)

# PCA
pca_dna <- compute_pca(matdata = scDNAmat)
plot(x = pca_dna$pca.data$PCA_1,pca_dna$pca.data$PCA_2,pch=16,col=grey(.3,.6))

# tSNE
set.seed(9)
tsne_dna <- compute_tsne(matdata = scDNAmat)

par(pty='s')
plot(x = tsne_dna$TSNE_1,y = tsne_dna$TSNE_2,pch=16,col=grey(.3,.6),main = paste('n. cells: ',nrow(tsne_dna)))

# scDNA clustering 

# hierarchical cluster model
fit_cluster_hierarchical <- hclust(dist(scale(tsne_dna[,2:3])))
tsne_dna$cl_hierarchical = factor(cutree(fit_cluster_hierarchical, k=3)) 

par(pty='s',mfrow=c(1,2))

colDNA <- c("black","red","green")
names(colDNA) <- 1:3
plot(x = tsne_dna$TSNE_1,y = tsne_dna$TSNE_2,pch=16,col=colDNA[tsne_dna$cl_hierarchical])

# assign cell-barcode to cell-id
per_cell_summary_metrics <- '/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/data/10xCNV/sample1-CNV/outs/per_cell_summary_metrics.csv'

barcodes_tab <- read.delim(per_cell_summary_metrics,header = T,as.is = T,stringsAsFactors = F,sep = ',')

scDNA_clustering_info <- merge(x = tsne_dna,y = barcodes_tab,by = 'cell_id',all.x = T)

write.table(scDNA_clustering_info,file = 'outs/Sample1_scDNA_clustering_info.tsv',col.names = T,row.names = F,quote = F,sep = '\t')

# get median cn per gene by cluster

clusters <- as.character(unique(tsne_dna$cl_hierarchical))

agg_scDNAmat <- matrix(nrow = nrow(scDNAmat),
                       ncol = length(clusters),
                       data = NA,
                       dimnames = list(rownames(scDNAmat), clusters))
for(i in clusters){
  cluster_barcodes <- tsne_dna$cell_id[grep(tsne_dna$cl_hierarchical,pattern = i)]
  agg_scDNAmat[,i] <- apply(scDNAmat[, cluster_barcodes ],MARGIN = 1,FUN = median)
}

boxplot(agg_scDNAmat[names(pca_dna$loading_scores[1:1000]),],outline=F,ylab="Median Gene Copy Number")

agg_scDNA_table <- as.data.frame(agg_scDNAmat,stringsAsFactors = F)
agg_scDNA_table <- cbind(rownames(agg_scDNA_table),agg_scDNA_table)
colnames(agg_scDNA_table)[1] <- 'gene_id'

write.table(agg_scDNA_table,file = 'outs/Sample1_scDNA_medianCN_clusters.tsv',col.names = T,row.names = F,quote = F,sep = '\t')

######################################################################################################################################################
## scRNA data

filtered_feature_bc_matrix <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/data/10XRNA5P/sample1PrimaryNuclei/outs/filtered_feature_bc_matrix/"

MatCounts <- read_filt_feature_bc_matrix(filtered_feature_bc_matrix = filtered_feature_bc_matrix,
                                         min.detected.genes = 1000,
                                         genes_umi_threshold = .2)

# load QC-pass clustering [ keep normal cells ] 
load("Seurat.objects.all.four.samples.RData")

head( seu.objects.list$sample1PrimaryNuclei@reductions$tsne@cell.embeddings )
head( seu.objects.list$sample1PrimaryNuclei@meta.data )

meta.data <- merge(x = seu.objects.list$sample1PrimaryNuclei@meta.data,
                   y = seu.objects.list$sample1PrimaryNuclei@reductions$tsne@cell.embeddings,
                   by = "row.names",all.x = T)

rownames(meta.data) <- paste0(meta.data$Row.names,"-1")
meta.data <- meta.data[,-1]
meta.data$cell_id <- rownames(meta.data)

head( meta.data )

nrm_counts <- MatCounts$nrm_counts
nrm_counts <- nrm_counts[, intersect(rownames(meta.data), colnames(nrm_counts)) ]

counts <- MatCounts$counts
counts <- counts[, colnames(nrm_counts) ]

col <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6')
names(col) <- 0:8
plot(x = meta.data$tSNE_1,y = meta.data$tSNE_2,pch=16,col=grey(.3,.1))
plot(x = meta.data$tSNE_1,y = meta.data$tSNE_2,pch=21,bg=col[meta.data$RNA_snn_res.0.6],col="black",lwd=.5)

# different approach to check within cluster expression and cn

head( agg_scDNAmat )

nrm_counts <- nrm_counts + 1
nrm_counts <- t(t(nrm_counts) / apply(nrm_counts,MARGIN = 1,FUN = median ))

gene.name_list <- rownames( agg_scDNAmat[names(pca_dna$loading_scores)[1:50],] )

par(pty="s",mfrow=c(1,2))
for(gene.name in gene.name_list){
  brPal <- colorRampPalette(c('grey','red'))
  if( gene.name %in% rownames(nrm_counts) ){
    message(gene.name)
    col <- data.frame(value=nrm_counts[gene.name,],
                      color=brPal(10)[as.numeric(cut(nrm_counts[gene.name,],breaks = 10))],
                      stringsAsFactors = F)
    df <- merge(x = meta.data,y = col,by = "row.names",all = F)
    boxplot(value~RNA_snn_res.0.6,data = df,varwidth=T,outline=T)
    plot(x = df$tSNE_1,y = df$tSNE_2,pch=19,col=paste(df$color, "80", sep=""),main = gene.name)
  }
}

# # tSNE
# set.seed(9)
# tsne_rna <- compute_tsne(matdata = nrm_counts)
# 
# plot(x = tsne_rna$TSNE_1,y = tsne_rna$TSNE_2,pch=16,col=grey(.3,.6))
# 
# # scRNA clustering
# 
# n.clusters <- 10
# 
# # hierarchical cluster model
# fit_cluster_hierarchical <- hclust(dist(scale(tsne_rna[,2:3])))
# tsne_rna$cl_hierarchical = cutree(fit_cluster_hierarchical, k=n.clusters)
# plot(x = tsne_rna$TSNE_1,y = tsne_rna$TSNE_2,pch=16,col=adjustcolor(tsne_rna$cl_hierarchical,.6))

# aggregate UMI counts by clusters

clusters <- as.character(unique(meta.data$RNA_snn_res.0.6))

agg_nrm_counts <- matrix(nrow = nrow(counts),
                         ncol = length(clusters),
                         data = NA,
                         dimnames = list(rownames(counts), clusters))

for(i in clusters){
  cluster_barcodes <- meta.data$cell_id[grep(meta.data$RNA_snn_res.0.6,pattern = i)]
  agg_nrm_counts[,i] <- rowSums(counts[, intersect(cluster_barcodes,colnames(counts)) ])
}

agg_nrm_counts <- t( t(agg_nrm_counts) / colSums(agg_nrm_counts) )

#save(pca_dna, agg_nrm_counts, agg_scDNAmat,file = "agg_rna_dna_matrices.RData",compress = T)

###########################
# compare scDNA vs scRNA

load("agg_rna_dna_matrices.RData")

dim(agg_nrm_counts)
dim(agg_scDNAmat)

loading_scores_genes <- names(pca_dna$loading_scores)[1:1000]

agg_scDNAmat <- agg_scDNAmat[ loading_scores_genes ,]
agg_nrm_counts <- agg_nrm_counts[intersect(rownames(agg_nrm_counts), loading_scores_genes ),]

common.genes <- intersect(rownames(agg_nrm_counts), rownames(agg_scDNAmat) )
agg_scDNAmat <- agg_scDNAmat[common.genes,]
agg_nrm_counts <- agg_nrm_counts[common.genes,]

agg_nrm_counts <- as.matrix( agg_nrm_counts + 1 )
agg_nrm_counts <- t(t(agg_nrm_counts) / apply(agg_nrm_counts,MARGIN = 1,FUN = mean ))

dim(agg_nrm_counts)
dim(agg_scDNAmat)

matR <- matrix(nrow = ncol(agg_nrm_counts), ncol = ncol(agg_scDNAmat), data = NA)

for(i in 1:ncol(agg_nrm_counts)){
  cat(paste("[",Sys.time(),"]\tsc-rna:",i,"\n"))
  exprData <- agg_nrm_counts[,i]
  out <- mclapply(1:ncol(agg_scDNAmat),FUN = getreg, exprData, agg_scDNAmat, mc.cores = 2 )
  matR[i,] <- unlist(out)
}

rownames(matR) <- colnames(agg_nrm_counts)
colnames(matR) <- colnames(agg_scDNAmat)

my_palette <- colorRampPalette(c("white", "forestgreen"))(n = 50)

RowSideColors <- col[rownames(matR)]
ColSideColors <- colDNA[colnames(matR)]
library(gplots)
heatmap.2(matR,
          ColSideColors = ColSideColors,
          RowSideColors = RowSideColors,
          trace="none",
          scale="none",
          col=my_palette,
          dendrogram = "both",
          labCol=FALSE,labRow = FALSE)

# assign rna to cluster
mat_data_max <- matrix(nrow = nrow(matR),ncol = ncol(matR),data = 0)
mat_data_max[ cbind(seq(nrow(mat_data_max)), apply(matR,MARGIN = 1,FUN = which.max)) ] <- 1

my_palette <- colorRampPalette(c("white", "black"))(n = 2)
heatmap.2(mat_data_max,
          scale="none",
          col=my_palette,trace="none",
          dendrogram = "column",
          labCol=FALSE,labRow = FALSE)

# plot to tSNE
meta.data$cl_based_on_cn <- NA
for(i in 1:nrow(meta.data)){
  tsnecluster <- as.character( meta.data$RNA_snn_res.0.6[i] )
  meta.data$cl_based_on_cn[i] <- colnames(matR)[which.max(matR[tsnecluster,])]
}

table(meta.data$cl_based_on_cn)
plot(x = meta.data$tSNE_1,y = meta.data$tSNE_2,pch=19,col=meta.data$cl_based_on_cn)



