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
load("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/check_gained_genes/masterDFgenemodel_diploid_annotation.Rdata")

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

plot(x = tsne_dna$TSNE_1,y = tsne_dna$TSNE_2,pch=16,col=grey(.3,.6))

# scDNA clustering 

# hierarchical cluster model
fit_cluster_hierarchical <- hclust(dist(scale(tsne_dna[,2:3])))
tsne_dna$cl_hierarchical = factor(cutree(fit_cluster_hierarchical, k=3)) 
plot(x = tsne_dna$TSNE_1,y = tsne_dna$TSNE_2,pch=16,col=tsne_dna$cl_hierarchical)

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

boxplot( agg_scDNAmat[names(pca_dna$loading_scores[1:100]),] )

###############
## scRNA data

filtered_feature_bc_matrix <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/data/10XRNA5P/sample1PrimaryNuclei/outs/filtered_feature_bc_matrix/"

MatCounts <- read_filt_feature_bc_matrix(filtered_feature_bc_matrix = filtered_feature_bc_matrix,
                                         min.detected.genes = 1000,
                                         genes_umi_threshold = .2)

# load QC-pass clustering [ keep normal cells ] 
load("Seurat.objects.all.four.samples.RData")

meta.data <- seu.objects.list$sample1PrimaryNuclei@meta.data

keep_barcodes <- paste0(rownames(meta.data),"-1")

nrm_counts <- MatCounts$nrm_counts
nrm_counts <- nrm_counts[, intersect(keep_barcodes, colnames(nrm_counts)) ]

counts <- MatCounts$counts
counts <- counts[,colnames(nrm_counts)]

# PCA
pca_rna <- compute_pca(matdata = nrm_counts)

plot(x = pca_rna$pca.data$PCA_1,pca_rna$pca.data$PCA_2,pch=16,col=grey(.3,.6))

# tSNE
set.seed(9)
tsne_rna <- compute_tsne(matdata = nrm_counts)

plot(x = tsne_rna$TSNE_1,y = tsne_rna$TSNE_2,pch=16,col=grey(.3,.6))

# scRNA clustering

n.clusters <- 10

# hierarchical cluster model
fit_cluster_hierarchical <- hclust(dist(scale(tsne_rna[,2:3])))
tsne_rna$cl_hierarchical = cutree(fit_cluster_hierarchical, k=n.clusters)
plot(x = tsne_rna$TSNE_1,y = tsne_rna$TSNE_2,pch=16,col=adjustcolor(tsne_rna$cl_hierarchical,.6))

# aggregate UMI counts by clusters

clusters <- as.character(unique(tsne_rna$cl_hierarchical))

agg_nrm_counts <- matrix(nrow = nrow(nrm_counts),
                         ncol = length(clusters),
                         data = NA,
                         dimnames = list(rownames(nrm_counts), clusters))

for(i in clusters){
  cluster_barcodes <- tsne_rna$cell_id[grep(tsne_rna$cl_hierarchical,pattern = i)]
  agg_nrm_counts[,i] <- rowSums(counts[, cluster_barcodes ])
}

agg_nrm_counts <- t( t(agg_nrm_counts) / colSums(agg_nrm_counts) )

save(agg_nrm_counts, agg_scDNAmat,file = "agg_rna_dna_matrices.RData",compress = T)

###########################
# compare scDNA vs scRNA

load("agg_rna_dna_matrices.RData")

dim(agg_nrm_counts)
dim(agg_scDNAmat)

loading_scores_genes <- names(pca_dna$loading_scores)

agg_scDNAmat <- agg_scDNAmat[loading_scores_genes,]

keep <- intersect(rownames(agg_nrm_counts), loading_scores_genes[1:100] )
agg_nrm_counts <- agg_nrm_counts[keep,]

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

my_palette <- colorRampPalette(c("white", "forestgreen"))(n = 50)

library(gplots)
heatmap.2(matR,
          scale="none",
          col=my_palette,trace="none",
          dendrogram = "both",
          labCol=FALSE,labRow = FALSE)




