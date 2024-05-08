<<<<<<< HEAD
=======
library(ALDEx2)
library(phyloseq)
library(ggplot2)
library(scales)
library(grid)
library(RColorBrewer)

>>>>>>> main
#Create a custom color scale
Rank2tax <- c("Alveolata", "Amoebozoa", "Apusozoa", "Archaeplastida", "Excavata", "Hacrobia", "Opisthokonta", "Rhizaria", "Stramenopiles", "Eukaryota_X", "Eukaryota", "Unassigned")
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
myPalette=getPalette(length(Rank2tax))
names(myPalette) <- Rank2tax

#function to plot significantly different ASVs
plot_results <- function(sigtab, tax1, tax2, title) {
  # Phylum order
  col_index1 = which(names(sigtab)==tax1)
  x = tapply(sigtab$diff.btw, sigtab[,col_index1], function(x) max(x))  #gets max OTU fold change for each phylum?
  x = sort(x, TRUE)       #sort, highest fold change first
  sigtab[,col_index1] = factor(as.character(sigtab[,col_index1]), levels=names(x))  #therefore, factor order is named by highest fold change first
  # Order order
  col_index2 = which(names(sigtab)==tax2)
  x = tapply(sigtab$diff.btw, sigtab[,col_index2], function(x) max(x))
  x = sort(x, TRUE)
  sigtab[,col_index2] = factor(as.character(sigtab[,col_index2]), levels=names(x))
  aldex2_plot <- ggplot(sigtab, aes(x=sigtab[,col_index2], y=diff.btw, color=sigtab[,col_index1])) + 
    geom_point(size=3) + 
    theme_set(theme_bw()) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5), legend.title=element_blank(), axis.title.x=element_blank()) +
    ggtitle(title) +
    xlab(tax2) + 
    scale_colour_manual(name=tax1,values = myPalette)
  #scale_fill_discrete(name=tax1)
  
  return(aldex2_plot)
}


<<<<<<< HEAD
library(ALDEx2)
library(phyloseq)
library(ggplot2)
library(scales)
library(grid)
library(RColorBrewer)

=======
>>>>>>> main

#from phyloseq object qiimez
#make data frame object with just otu counts
# = otu_table(qiimez)
otu_table_df <- as.data.frame(otu_table(qiimez))
#need vector of conditions in same order as samples in otu_table
sample_type_vec <- sample_data(qiimez)$sample_type

#BUT, will want samples only for pairwise comparison of 2 groups (here, interstitial and sand, so all except sea water samples)
qiimez_interstitial_sand <- subset_samples(qiimez, sample_type!="sea water")
#make data frame object with just otu counts
otu_table_is_df <- as.data.frame(otu_table(qiimez_interstitial_sand))
#need vector of groups(conditions) in same order as samples in otu_table
sample_type_is_vec <- sample_data(qiimez_interstitial_sand)$sample_type

#generate instances of centered-log ratio transformed values using monte carlo samples of the dirichlet distribution
#denom = "iqlr" for assymetric samples, ie clear assymetry in groups, systematic differences, one group more 0s than the other
x <- aldex.clr(otu_table_is_df, sample_type_is_vec, mc.samples=128, denom="iqlr", verbose=TRUE)

#performs the Welch’s t and Wilcoxon rank test for the instance when there are only two conditions
x.tt <- aldex.ttest(x, paired.test=FALSE)
#estimate effect size and the within and between condition values in the case of two conditions
x.effect <- aldex.effect(x, include.sample.summary=FALSE, verbose=TRUE)
#t-test and effect data are merged into one object.
x.all <- data.frame(x.tt,x.effect)

#plot
#which type of plot is to be produced. MA is a Bland-Altman style plot; MW is a difference between to a variance within plot as described in the paper
#default Benjamini-Hochberg fdr cutoff q <= 0.1, if < 0.1, points are red
par(mfrow=c(1,2))
aldex.plot(x.all, type="MA", test="welch", cutoff=0.05)
aldex.plot(x.all, type="MW", test="welch", cutoff=0.05)

par(mfrow=c(1,2))
aldex.plot(x.all, type="MA", test="wilcox")
aldex.plot(x.all, type="MW", test="wilcox")

#summary(x.all)

#get OTUs with q<= 0.1
# identify which values are significant in both the t-tests Welch's and Wilcox
found.by.all <- which(x.all$we.eBH < 0.05 & x.all$wi.eBH < 0.05)
# identify which values are significant in either test (then plot underneath)
found.by.one <- which(x.all$we.eBH < 0.05 | x.all$wi.eBH < 0.05)
# plot the within and between variation of the data
plot(x.all$diff.win, x.all$diff.btw, pch=19, cex=0.3, col=rgb(0,0,0,0.3), xlab="Difference within", ylab="Difference between")
points(x.all$diff.win[found.by.one], x.all$diff.btw[found.by.one], pch=19, cex=0.5, col=rgb(0,0,1,0.5))
points(x.all$diff.win[found.by.all], x.all$diff.btw[found.by.all], pch=19, cex=0.5, col=rgb(1,0,0,1))
abline(0,1,lty=2)
abline(0,-1,lty=2)

#just get welch's significant results
sig_we <- x.all[(x.all$we.eBH < 0.05), ]
sig_we <- sig_we[order(sig_we$diff.btw, na.last=NA), ]
sig_we <- cbind(as(sig_we,"data.frame"), as(tax_table(qiimez)[rownames(sig_we), ], "matrix"))
sig_we_plot <- plot_results(sig_we, "Rank2", "Rank2", "18SV4 aldex2, protists, sand vs interstitial")
#execute to view plot
sig_we_plot
#save plot
ggsave("18SV4_protists_aldex2_sand_interstitial.pdf")


#output significant OTUs
header = as.matrix(t(c("OTU_ID", colnames(sig_we))))
write.table(header, "18SV4_protists_aldex2_sand_interstitial_sig.txt", row.names=F, append=F, quote= FALSE, sep="\t", col.names=F)
write.table(sig_we,"18SV4_protists_aldex2_sand_interstitial_sig.txt", row.names=T, append=T, quote= FALSE, sep="\t", col.names=F)


