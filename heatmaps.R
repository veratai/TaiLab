#sample code for making heat maps from microbiome data
#from Tai lab

#heat map, using base heatmap plotting
library(dendextend)
library(gplots)

#cluster based on unifrac distances
qiimez.r = rarefy_even_depth(qiimez)
qiimez.r.uwuf <- distance(qiimez.r, "unifrac")
qiimez.r.uwuf.hclust <- hclust(qiimez.r.uwuf, method="average")
plot(qiimez.r.uwuf.hclust)
summary(qiimez.r.uwuf.hclust)
#relabel samples, for labeling
qiimez.r.uwuf.hclust.beachsite <- qiimez.r.uwuf.hclust
qiimez.r.uwuf.hclust.beachsite$labels <- sample_data(rank2glommed.relA)$beach_site
#rotate left most branch, and high interstitial samples
plot(rotate(qiimez.r.uwuf.hclust.beachsite, c(2:15,22:85,16:21,1)))
dendro <- as.dendrogram(rotate(qiimez.r.uwuf.hclust.beachsite, c(2:15,22:85,16:21,1)))

#convert relative abundance info to dataframe
df <- as.data.frame(otu_table(rank2glommed.relA), row.names=tax_table(rank2glommed.relA)[,2], stringsAsFactors=FALSE)
#heatmap(as.matrix(df), Colv=dendro)

#rename column names, for labeling
colnames(df) <- sample_data(rank2glommed.relA)$beach_site

#make distance matrix for rel abundance data, to cluster taxa with similar rel abundances
df_dist <- dist(as.matrix(df))
df_hclust <- hclust(df_dist, method="average")
plot(df_hclust)
#rotate right most branches
plot(rotate(df_hclust, c(1:2,10:11,3:9)))
df_dendro_rotate <- as.dendrogram(rotate(df_hclust, c(1:2,10:11,3:9)))

#make 0 values into NAs, then can assign colour to NAs in heatmap
df.na <- df
df.na[df.na == 0] <- NA

#make colour bar matrix
#c("#009E73", "#E69F00", "#56B4E9")  #sand, interstitial, seawater
sample_col <- sample_data(rank2glommed.relA)$sample_type
sample_col <- gsub("beach sand only", "#009E73", sample_col)
sample_col <- gsub("beach interstitial", "#E69F00", sample_col)
sample_col <- gsub("sea water", "#56B4E9", sample_col)

#plot heatmap (relative abundance, hclustered by unifrac distances)heat_palette <- colorRampPalette(c(brewer.pal(9, "YlGnBu")[2:9]), bias=3)(n=500)
pdf("18SV4_protists_relA_rank2_clust_heatmap.pdf")
heatmap.2(as.matrix(df.na), Rowv=df_dendro_rotate, Colv=dendro, col=heat_palette, na.color="white", 
          trace="none", cexRow=0.8, cexCol=0.6, lhei=c(1,5),margins=c(5,8),keysize=1,
          key.title="relative abundance", key.xlab=NA, density.info="none",
          ColSideColors=sample_col
)
dev.off()
