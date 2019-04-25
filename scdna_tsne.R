setwd("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_integration_dna_rna/scrna_vs_scdna")

library( data.table )
library( Rtsne )
library( parallel )

loadCells <- function(j,cellids, max.cn = 10){
  cellRaw <- fread(cellids[j],header = F,stringsAsFactors = F,data.table = F,select = c(5,9),colClasses=list(character=5,numeric=9))
  cell <- as.numeric(cellRaw[,2])
  names(cell) <- cellRaw[,1]
  cell <- cell[ cell <= max.cn ]
  return( cell )
}

cellids <- list.files("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/sample1.genemodel/tmp.parallel",pattern = "intcellids_",full.names = T)

ListCells <- lapply(1:length(cellids), loadCells, cellids)

getUniqueGenes <- function(j,ListCells){
  return( unique(names(ListCells[[j]])) )
}

genes <- unique(unlist( lapply(1:length(ListCells), getUniqueGenes, ListCells) ))

DTgenes <- data.table(gene.id = genes,stringsAsFactors = F)
rownames(DTgenes) <- DTgenes$gene.id

fillGenesCellsMatrix <- function(j, ListCells, DTgenes, genes){
  cn <- ListCells[[j]]
  dna <- data.table(gene.id = names(cn), copy.number = as.numeric(cn),stringsAsFactors = F)
  m <- merge(x = DTgenes,y = dna,by = "gene.id",all.x = T)
  m <- m[genes,]
  return( as.numeric(m$copy.number) )
}

columns <- mclapply(1:length(ListCells), fillGenesCellsMatrix, ListCells, DTgenes, genes, mc.cores = 4)

scDNAmat <- do.call(cbind, columns)
rownames(scDNAmat) <- genes
colnames(scDNAmat) <- gsub(basename(cellids),pattern = "intcellids_|-genemodel.sort.bed",replacement = "")

save(ListCells, genes, scDNAmat, file = "scDNAmat.RData",compress = T)

# clean data
load("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/check_gained_genes/masterDFgenemodel_diploid_annotation.Rdata")

keep <- as.character(master$id[which(master$flag10x_is.noisy == 0 & master$loss < .5)])

scDNAmat <- scDNAmat[, keep ]
dim(scDNAmat)

scDNAmat <- scDNAmat[rowSums(is.na(scDNAmat)) == 0, ]

# PCA
pca <- prcomp(x = t(scDNAmat))

par(pty="s")
plot(pca$x[,1],pca$x[,2],pch=16,col=grey(0.3,0.2))

pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)
barplot(pca.var.per[1:10],ylim=c(0,100))

library(ggplot2)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
pca.data

ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("My PCA Graph")

## get the name of the top 10 measurements (genes) that contribute
## most to pc1.
loading_scores <- pca$rotation[,1]
gene_scores <- abs(loading_scores) ## get the magnitudes
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
top_10_genes <- names(gene_score_ranked[1:10])

# tSNE

# UMAP


