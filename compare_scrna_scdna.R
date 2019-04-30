setwd("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_integration_dna_rna/scrna_vs_scdna")

load("agg_nrm_counts.RData")
load("agg_scDNAmat.RData")

dim(agg_nrm_counts)
dim(agg_scDNAmat)

common.genes <- intersect(rownames(agg_nrm_counts), rownames(agg_scDNAmat) )

agg_scDNAmat <- agg_scDNAmat[common.genes,]
agg_nrm_counts <- agg_nrm_counts[common.genes,]

agg_nrm_counts <- as.matrix( agg_nrm_counts + 1 )
agg_nrm_counts <- t(t(agg_nrm_counts) / apply(agg_nrm_counts,MARGIN = 1,FUN = mean ))

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
