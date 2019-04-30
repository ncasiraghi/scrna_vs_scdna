setwd("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_integration_dna_rna/scrna_vs_scdna")

library( data.table )
library( Rtsne )
library( parallel )

if(FALSE){
  
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

}

# upload
load("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_integration_dna_rna/scrna_vs_scdna/scDNAmat.RData")

# clean dna data
load("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/check_gained_genes/masterDFgenemodel_diploid_annotation.Rdata")

keep <- as.character(master$id[which(master$flag10x_is.noisy == 0 & master$loss < .5)])

scDNAmat <- scDNAmat[, keep ]
dim(scDNAmat)

scDNAmat <- scDNAmat[rowSums(is.na(scDNAmat)) == 0, ]

# PCA
pca <- prcomp(x = t(scDNAmat))

pca.data <- data.frame(cell_id=rownames(t(scDNAmat)), PCA_1=pca$x[,1], PCA_2=pca$x[,2], stringsAsFactors = F)

plot(x = pca.data$PCA_1,pca.data$PCA_2,pch=16,col=grey(.3,.6))

# tSNE
set.seed(9) 

tsne = Rtsne(t(scDNAmat), check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=.5, dims=2)

tsne.data = data.frame(cell_id=rownames(t(scDNAmat)), TSNE_1=tsne$Y[,1], TSNE_2=tsne$Y[,2], stringsAsFactors = F)  

plot(x = tsne.data$TSNE_1,y = tsne.data$TSNE_2,pch=16,col=grey(.3,.6))


