setwd("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_integration_dna_rna/scrna_vs_scdna")

library( Matrix )
library( Rtsne )
library( tidyverse )

read_filt_feature_bc_matrix <- function(filtered_feature_bc_matrix,min.detected.genes,genes_umi_threshold){
  
  genes <- read_tsv(file = file.path(filtered_feature_bc_matrix,"features.tsv.gz"), col_names = c( "gene_id", "gene_name", "info") )
  
  duplicate_gene_names <- genes$gene_name[ duplicated(genes$gene_name) ]
  
  genes$gene_name <- ifelse(
    genes$gene_name %in% duplicate_gene_names,
    str_c( genes$gene_name, "_", genes$gene_id ),
    genes$gene_name )
  
  barcodes <- readLines(file.path(file.path(filtered_feature_bc_matrix,"barcodes.tsv.gz")))
  
  matrixTable <- read_delim( file.path(file.path(filtered_feature_bc_matrix,"matrix.mtx.gz")),
                             delim=" ", comment = "%", col_names = c( "gene", "barcode", "count" ) )
  
  counts <-
    sparseMatrix( 
      i = matrixTable$gene[-1], j = matrixTable$barcode[-1], x = matrixTable$count[-1],
      dims = c( matrixTable$gene[1], matrixTable$barcode[1] ),
      dimnames = list( gene = genes$gene_id, barcode = barcodes ) )
  
  rm( barcodes, genes, matrixTable )
  
  # call everything with at least 1000 detected genes (above the gray line) a cell.
  cell_barcodes <- names(which( colSums( counts>0 ) >= min.detected.genes ))
  counts <- counts[ , cell_barcodes ]
  
  # remove the genes that appear in none of the remaining cells
  counts <- counts[ rowSums(counts) > 0,]
  
  # Normalisation
  nrm <- t( t(counts) / colSums(counts) )
  
  # keep only genes having UMI count not null in at leat x fraction of cells
  genes_umi <- names(which( (rowSums( nrm > 0 )/ncol(nrm)) >= genes_umi_threshold))
  
  nrm <- nrm[ genes_umi ,]
  counts <- counts[ genes_umi ,]
  
  return(list(as.matrix(counts),as.matrix(nrm)))
  
}

filtered_feature_bc_matrix <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/data/10XRNA5P/sample1PrimaryNuclei/outs/filtered_feature_bc_matrix/"

MatCounts <- read_filt_feature_bc_matrix(filtered_feature_bc_matrix = filtered_feature_bc_matrix,
                                         min.detected.genes = 1000,
                                         genes_umi_threshold = .4)

counts <- MatCounts[[1]]
nrm_counts <- MatCounts[[2]]
   
par(pty="s",mfrow=c(1,2))

# PCA
pca <- prcomp(x = t(nrm_counts))

pca.data <- data.frame(cell_id=rownames(t(nrm_counts)), PCA_1=pca$x[,1], PCA_2=pca$x[,2], stringsAsFactors = F)

plot(x = pca.data$PCA_1,pca.data$PCA_2,pch=16,col=grey(.3,.6))

# tSNE
set.seed(9) 

tsne = Rtsne(t(nrm_counts), check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=.5, dims=2)

tsne.data = data.frame(cell_id=rownames(t(nrm_counts)), TSNE_1=tsne$Y[,1], TSNE_2=tsne$Y[,2], stringsAsFactors = F)  

plot(x = tsne.data$TSNE_1,y = tsne.data$TSNE_2,pch=16,col=grey(.3,.6))

# scRNA clustering 

n.clusters <- 10

# k-means clustering model
fit_cluster_kmeans <- kmeans(scale(tsne$Y),n.clusters)  
tsne.data$cl_kmeans <- fit_cluster_kmeans$cluster
plot(x = tsne.data$TSNE_1,y = tsne.data$TSNE_2,pch=16,col=adjustcolor(tsne.data$cl_kmeans,.6))

# hierarchical cluster model
fit_cluster_hierarchical <- hclust(dist(scale(tsne$Y)))
tsne.data$cl_hierarchical = cutree(fit_cluster_hierarchical, k=n.clusters)
plot(x = tsne.data$TSNE_1,y = tsne.data$TSNE_2,pch=16,col=adjustcolor(tsne.data$cl_hierarchical,.6))

# aggregate UMI counts by clusters

clusters <- as.character(unique(tsne.data$cl_hierarchical))

agg_nrm_counts <- matrix(nrow = nrow(nrm_counts),
                         ncol = length(clusters),
                         data = NA,
                         dimnames = list(rownames(nrm_counts), clusters))

for(i in clusters){
  cluster_barcodes <- tsne.data$cell_id[grep(tsne.data$cl_hierarchical,pattern = i)]
  agg_nrm_counts[,i] <- rowSums(counts[, cluster_barcodes ])
}

agg_nrm_counts <- t( t(agg_nrm_counts) / colSums(agg_nrm_counts) )

save(agg_nrm_counts,file = "agg_nrm_counts.RData",compress = T)




