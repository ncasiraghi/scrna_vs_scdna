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

getUniqueGenes <- function(j,ListCells){
  return( unique(names(ListCells[[j]])) )
}

fillGenesCellsMatrix <- function(j, ListCells, DTgenes, genes){
  cn <- ListCells[[j]]
  dna <- data.table(gene.id = names(cn), copy.number = as.numeric(cn),stringsAsFactors = F)
  m <- merge(x = DTgenes,y = dna,by = "gene.id",all.x = T)
  m <- m[genes,]
  return( as.numeric(m$copy.number) )
}

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
  
  return(list(counts=as.matrix(counts),
              nrm_counts=as.matrix(nrm)))
  
}

compute_pca <- function(matdata,get_loadscore_from_pc = 1){
  pca <- prcomp(x = t(matdata))
  pca.var <- (pca$sdev^2)/sum(pca$sdev^2)
  head( pca.var )
  loading_scores <- sort(abs(pca$rotation[,get_loadscore_from_pc]),decreasing = T)
  pca.data <- data.frame(cell_id=rownames(t(matdata)), PCA_1=pca$x[,1], PCA_2=pca$x[,2], stringsAsFactors = F)
  return(list(pca=pca,
              pca.var=pca.var,
              loading_scores=loading_scores,
              pca.data=pca.data))
}

compute_tsne <- function(matdata){
  tsne = Rtsne(t(matdata), check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=.5, dims=2)
  tsne.data = data.frame(cell_id=rownames(t(matdata)), TSNE_1=tsne$Y[,1], TSNE_2=tsne$Y[,2], stringsAsFactors = F)
  return(tsne.data)
}

getreg <- function(j,exprData,agg_scDNAmat){
  cell <- agg_scDNAmat[,j]
  common.genes <- sort(intersect(names(exprData),names(cell)))
  exprData <- exprData[common.genes]
  cell <- cell[common.genes]
  m.regression <- summary(lm(exprData ~ cell))
  if(F){
    par(pty="s")
    plot(x = cell ,y = exprData,cex.main=0.6,pch = 16, col = gray(.3, .6),ylab="norm UMI counts",xlab="copy number")
    abline(h = mean(exprData,na.rm=T),col="gold")
    abline(m.regression,col="blue")
    mtext(paste0("R2 = ",m.regression$r.squared),cex = 0.5)
  }
  return(m.regression$r.squared)
}
