load('/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/frozen_scRNA_seurat_obj/12062019/Seurat.objects.all.four.samples.RData')

labels <- read.delim('inputs/labels.tsv',as.is = T,header = T,stringsAsFactors = F)
labels$barcode <- gsub(labels$barcode,pattern = '-1',replacement = '')

meta.data <- seu.object.list.all$sample1PrimaryNuclei@meta.data
meta.data$barcode <- rownames(meta.data)

tsne <- data.frame(seu.object.list.all$sample1PrimaryNuclei@reductions$tsne@cell.embeddings,stringsAsFactors = F)
tsne$barcode <- rownames(tsne)

seurat_out <- merge(x = meta.data[,c('barcode','seurat_clusters')],y = tsne,by = 'barcode',all.x = T)
seurat_out$barcode <- paste0(seurat_out$barcode,'-1')

# labels Seurat clustering and tSNE coordinates
write.table(x = seurat_out,'inputs/sample1_scRNA_Seurat_clusters.tsv',col.names = T,row.names = F,quote = F,sep = '\t')

tsne <- merge(x = tsne,y = labels,by = 'barcode',all.x = T)
# plot to tSNE color by Seurat clustering
plot(x = tsne$tSNE_1,y = tsne$tSNE_2,pch=19,col='grey')

# plot to tSNE color by MLE 

cols <- c('#fdae61','#abd9e9','#66bd63')
barplot(table(tsne$cluster),col = cols)

names(cols) <- sort(unique(labels$cluster))
plot(x = tsne$tSNE_1,y = tsne$tSNE_2,pch=21,lwd=0.2,bg=cols[as.character(tsne$cluster)])
legend('bottomleft',legend = sort(unique(labels$cluster)),fill=cols)
