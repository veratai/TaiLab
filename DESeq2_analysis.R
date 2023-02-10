#DESeq analysis

# function to calculate geometric mean over the normalized counts for each ASV
# 0 values are ignored
# as a work around to estimateSizeFactors 

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


#function to order results based on adj p-value, and return significant results (adjusted p-value < 0.01)

order_results <- function(results, phyloseq) {
  results = results[order(results$padj, na.last=NA), ]   #order by the adjusted p-value, remove entries with an NA value
  alpha = 0.01  #set adjusted p-value alpha
  sigtab = results[(results$padj < alpha), ]      #get only results with adj p-value less than alpha
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyloseq)[rownames(sigtab), ], "matrix"))       #organize data into nice format
  return(sigtab) 
}

#function to plot results

plot_results <- function(sigtab, tax1, tax2, title) {
  # tax1 order, to set factors for colours, can also adjust code to use for plotting order
  col_index1 = which(names(sigtab)==tax1)
  x = tapply(sigtab$log2FoldChange, sigtab[,col_index1], function(x) max(x))  #gets max OTU fold change for each tax1 group?
  x = sort(x, TRUE)       #sort, highest fold change first
  sigtab[,col_index1] = factor(as.character(sigtab[,col_index1]), levels=names(x))  #therefore, factor order is named by highest fold change first
  # tax2 order, to plot taxa in the order of having the ASV with highest fold change first
  col_index2 = which(names(sigtab)==tax2)
  x = tapply(sigtab$log2FoldChange, sigtab[,col_index2], function(x) max(x))
  x = sort(x, TRUE)
  sigtab[,col_index2] = factor(as.character(sigtab[,col_index2]), levels=names(x))
  #make plot
  deseq2_plot <- ggplot(sigtab, aes(x=sigtab[,col_index2], y=log2FoldChange, color=sigtab[,col_index1])) + 
    geom_point(size=3) + 
    theme_set(theme_bw()) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5), legend.title=element_blank(), legend.position="bottom", axis.title.x=element_blank()) +
    ggtitle(title) +
    xlab(tax2) + 
    #scale_colour_manual(name=tax1,values = myPalette)  #calls colour by specific taxon name
    guides(col=guide_legend(ncol=8))
  
  return(deseq2_plot)
}



library(phyloseq)
library(ggplot2)
library(scales)
library(grid)
library("DESeq2")
library(plyr)

#normalize data
#rarefying counts produces decent results, can also normalize counts with cumulative sum scaling
#see Weiss et al. 2017

##rarefy counts with phyloseq object
qiimez.rare = rarefy_even_depth(qiimez)


#deseq2 analysis based on sample_type groups
diagdds = phyloseq_to_deseq2(qiimez.rare, ~ sample_type)
#calculating geometric mean to be used in estimating size factors (for each ASV across samples)
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="local")

#examine results - contrasting groups in pairs (eg. sand (numerator) to seawater (denominator))
res_ssw = results(diagdds, contrast=c("sample_type", "sand", "sea water"))
#get results that are significant with adj p < 0.01
sigtab_ssw = order_results(res_ssw, qiimez.rare)

#plot results, choosing which taxonomic ranks for colouring points (eg. Rank2) and grouping, ordering ASVs in the plot (eg. Rank6)  
wald_ssw_sig_plot <- plot_results(sigtab_ssw, "Rank2", "Rank6", "18SV4 protists CSSnorm, sand vs seawater")
wald_ssw_sig_plot
ggsave("18SV4_protists_CSS_ddseqWald_sand_seawater.pdf")


#Can create a custom color scale for plot, eg.

#to list rank_names used in taxonomy
rank_names(qiimez)
#list unique taxa names in Rank2, use same rank as in plot_results function
Rank2tax <- get_taxa_unique(qiimez, "Rank2")
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
myPalette=getPalette(length(Rank2tax))
names(myPalette) <- Rank2tax
#then uncomment scale_colour_manual setting in plot_results function (and re-define function)


#output results in tabular format
#get results with positive fold change
posigtab_ssw = sigtab_ssw[sigtab_ssw[, "log2FoldChange"] > 0, ]
posigtab_ssw = posigtab_ssw[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6", "Rank7")]
head(posigtab_ssw, n=5)
#negative fold change
negsigtab_ssw = sigtab_ssw[sigtab_ssw[, "log2FoldChange"] < 0, ]
negsigtab_ssw = negsigtab_ssw[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6", "Rank7")]
head(negsigtab_ssw, n=5)
#write out
header = as.matrix(t(c("OTU_ID", colnames(posigtab_ssw))))
write.table(header, "18SV4_protists_CSS_ddseqWald_sand_seawater_results.txt", row.names=F, append=F, quote= FALSE, sep="\t", col.names=F)
write.table(posigtab_ssw,"18SV4_protists_CSS_ddseqWald_sand_seawater_results.txt", row.names=T, append=T, quote= FALSE, sep="\t", col.names=F)
write.table(negsigtab_ssw,"18SV4_protists_CSS_ddseqWald_sand_seawater_results.txt", row.names=T, append=T, quote= FALSE, sep="\t", col.names=F)

