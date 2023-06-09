#taxonomy bar plots
library(qiime2R)
library(phyloseq)
library(plyr)
library(ggplot2)
library(RColorBrewer)


#set working directory
setwd("/Users/vera/Projects/")

#import qiime2 files as a phyloseq object named qiimez
qiimez <- qza_to_phyloseq(
  features="asvtable.qza",
  tree="unrooted-tree.qza",
  taxonomy="taxonomy.qza",
  metadata = "metadata.txt"
)

#list rank_names used in taxonomy
rank_names(qiimez)
#list unique taxa names for a given rank, e.g. for Rank2
get_taxa_unique(qiimez, "Rank2")
#can make bar plot of total asvs for a taxonomic rank
plot_bar(qiimez, fill="Rank2")

# may want to re-label some of the taxa names
#e.g. for Rank 2 replace "Eukaryota_X" with "Eukaryota"
tax_table(qiimez)[,2][tax_table(qiimez)[,2]=="Eukaryota_X"] <- "Eukaryota"

#glom Rank2 ASVs together (so in plot_bar will not show lines for individuals ASVs)
#i.e. sums all ASVs that have the same classification at Rank2 
rank2glommed = tax_glom(qiimez, "Rank2")

#get relative abundance
rank2glommed.relA <- transform_sample_counts(rank2glommed, function(x) x / sum(x))

#melt data (i.e. covert data table into long format) and turn into dataframe to control levels
rank2glommed.relA.df <- psmelt(rank2glommed.relA)

#make custom colour palette, using colour brewer  
myPalette = c(brewer.pal(length(levels(rank2glommed.relA.df$Rank2)), "Paired"))

#plot
p <- ggplot(rank2glommed.relA.df, aes(x=SampleID, y=Abundance, fill=Rank2, order = as.factor(Rank2)))

p + 
  geom_bar(stat="identity", colour="black") +
  scale_fill_manual(values=myPalette) + 
  ylab("relative abundance") +
  theme_set(theme_bw()) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5), legend.title=element_blank(), axis.title.x=element_blank()) +
  ggtitle("relative abundance") +
  #facet_wrap(~Rank2, scales="free_y")

#save plot as pdf
ggsave("Rank2_barplot.pdf", width=10, height=5, units="in")



#before plotting, can get mean relative abundance, standard deviation for groups of samples
#the code below is customized for an analysis of beach samples, so you'll have to change the variable names as needed for your data

#get average relA for each beach_site
rank2_beach_site_mean <- lapply( split(rank2glommed.relA.df$Abundance, f=list(rank2glommed.relA.df$Rank2, rank2glommed.relA.df$beach_site)), mean )
rank2_beach_site_relA <- as.data.frame(unlist(rank2_beach_site_mean))
colnames(rank2_beach_site_relA) <- "meanRelA"

#get standard deviation of relA for each beach_site
rank2_beach_site_sd <- lapply( split(rank2glommed.relA.df$Abundance, f=list(rank2glommed.relA.df$Rank2, rank2glommed.relA.df$beach_site)), sd )
rank2_beach_site_relA$sdRelA <- unlist(rank2_beach_site_sd)

#need to recover beach_site and Rank2 data, they are merged together in rownames
split_names <- strsplit(rownames(rank2_beach_site_relA), "\\.")
rank2 <- sapply(split_names, function(x) x[1])
beach_site <- sapply(split_names, function(x) x[2])

#add to datatable, order levels - to get them plotted in a customized order
rank2_beach_site_relA$rank2 <- factor(rank2, levels=c("Alveolata", "Amoebozoa", "Apusozoa", "Archaeplastida", "Excavata", "Hacrobia", "Opisthokonta", "Rhizaria", "Stramenopiles", "Eukaryota", "Unassigned"))
rank2_beach_site_relA$beach_site <- factor(beach_site, 
                                           levels=c("WB_high_so", "2nd_high_so", "WB_mid_so", "2nd_mid_so", "WB_mid_gw_so", "WB_low_so", "2nd_low_so", "WB_swash_so", "2nd_swash_so",
                                                    "NB_high", "3rd_high", "7th_high", 
                                                    "NB_mid", "3rd_mid", "7th_mid", "7th_mid_bar", "WB_mid_gw",
                                                    "3rd_mid_sw", "7th_mid_sw",
                                                    "NB_low", "WB_low", "3rd_low", "7th_low",
                                                    "NB_swash", "WB_swash", "3rd_swash", "7th_swash",
                                                    "NB_sw", "WB_sw", "2nd_sw", "3rd_sw", "7th_sw")
                                            )   

# so, now use rank2_beach_site_relA data frame to plot mean relative abundance
# use facet wrap to plot as individual bars, and not stacked bar plot
