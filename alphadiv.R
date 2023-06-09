#Example code for analyzing alpha and beta diversity

#load relevant libraries
library(qiime2R)
library(phyloseq)
library(ggplot2)
library(plyr)
library(RColorBrewer)

#set working directory
setwd("/Users/vera/Projects/")

#qiime2 files imported as a phyloseq object named qiimez
qiimez <- qza_to_phyloseq(
  features="asvtable.qza",
  tree="unrooted-tree.qza",
  taxonomy="taxonomy.qza",
  metadata = "metadata.txt"
)


#to list rank_names used in taxonomy
rank_names(qiimez)
#list unique taxa names in Rank2
get_taxa_unique(qiimez, "Phylum")

#just get data for a particular group
qiimez_cyanos <- subset_taxa(qiimez, Phylum=="Cyanobacteria")
get_taxa_unique(qiimez_cyanos, "Class")
plot_bar(qiimez_cyanos, fill="Class")


#Alpha diversity

#without rarefying
#calculate diversity indices
qiimez.rich <- estimate_richness(qiimez, measures=c("Chao1", "Shannon"))
#add sample_type info (or other metadata info) to richness table
qiimez.rich$sample_type <- sample_data(qiimez)$sample_type

## plot results:
#all values plotted
plot_richness(qiimez.rich, measures=c("Chao1", "Shannon"), x="sample_type", color="sample_type") +
  geom_point(size=1, alpha=0.7) + 
  scale_colour_manual(values=c("beach interstitial"="seagreen3", "beach sand only"="sienna1", "sea water"="cornflowerblue")) +     # control colour by sample_type
  ggtitle("qiimez.r_alphadiv")
ggsave("qiimez.r_alphadiv.pdf", width=7, height=7)

## can also plot like this:
# can avoid using plot_richness function, and plot directly with ggplot
ggplot(phy.rich, aes(x=sample_type, y=Shannon, color=sample_type)) +
  geom_boxplot() + 
  #scale_colour_manual(values=c("beach interstitial"="seagreen3", "beach sand only"="sienna1", "sea water"="cornflowerblue")) +     # control colour by sample_type
  ggtitle("qiimez.r_alphadiv")


###typically will rarify data for even sampling depth, sampling number is to sample with fewest reads
## want to check number of reads per sample to check if rarefication is needed/appropriate
qiimez.rare <- rarefy_even_depth(qiimez)

#calculate diversity indices
qiimez.rarerich <- estimate_richness(qiimez.rare, measures=c("Chao1", "Shannon"))

#add sample_type info to richness table
qiimez.rarerich$sample_type <- sample_data(qiimez.rare)$sample_type

#plot results for one specific metric
ggplot(qiimez.rarerich, aes(x=sample_type_order, y=Shannon, fill=sample_type)) +
  geom_boxplot() +
  ylab("Shannon diversity index")+
  scale_x_discrete(labels=c("sand","interstitial","seawater")) +   # to control order to categories on the x-axis
  scale_y_continuous(limits=c(0.5, 7.5), breaks=c(1, 2, 3, 4, 5, 6, 7, 8)) +    # to control range of values on y-axis
  scale_fill_manual(values=c("#009E73", "#E69F00", "#56B4E9")) +    # to control colour
  theme_set(theme_bw()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x.bottom = element_text(hjust = 1)) +
  theme(axis.text.x = element_text(size=15, angle=45, colour="black"), axis.text.y = element_text(size=12, colour="black"), axis.title.y = element_text(size=15), axis.title.x = element_blank()) +
  #theme(legend.position="none") +     # to remove legend from plot
  ggtitle("qiimez.r.alpha_diversity")
#export plot as pdf
ggsave("qiimez.r.alpha_diversity.pdf", width=5, height=7, useDingbats=FALSE) 

#e.g. T-tests, check if alpha diversity is significantly different

int <- subset(qiimez.rarerich, sample_type == "beach interstitial")
sand <- subset(qiimez.rarerich, sample_type == "beach sand only")
sw <- subset(qiimez.rarerich, sample_type == "sea water")

t.test(int$Shannon, sand$Shannon)

#or can do anova, and post-hoc tests



