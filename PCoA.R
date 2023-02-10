# code to do ordination via PCoA, through phyloseq functions
# test for significance of groups using PERMANOVA
# plot PCoA
# and also add species score arrows to these plots 
# based on this forum post: https://github.com/joey711/phyloseq/issues/515)


library(qiime2R)
library(phyloseq)
library(ggplot2)


setwd("/Users/veratai/Projects/SunfishLakePlanktothrix/cyano16S_analysis/")

sunfish <- qza_to_phyloseq(
  features="sunfish_cyano16S_filteredtable_ASVmin7.qza",
  tree="unrooted-tree-sunfish-cyano16S.qza",
  taxonomy="sunfish_cyano16S_taxonomy.qza",
  metadata = "sunfish_metadata_VT.txt"
)

#useful to get best taxonomic hit for an ASV
#can try using microbiomeutilities get_best_hit function (not sure what it is called)


#example of subsetting data based on a taxonomic group
#just get oxyphoto data
get_taxa_unique(sunfish, "Class")
sunfish_oxyphoto <- subset_taxa(sunfish, Class=="Oxyphotobacteria")
get_taxa_unique(sunfish_oxyphoto, "Order")
#there are some NAs, will eventually need to find a proper classification for these
#make barplot of relative abundances
plot_bar(sunfish_oxyphoto, fill="Family")


#but will just work with whole dataset here:

#rarefy to normalize sampling depth for all samples
sunfish.r <- rarefy_even_depth(sunfish)


#can specify custom shapes and colours manually, for example:
shape.beach <- c(8, 1, 17, 18, 16)
colScale.sample_type <- c("#00BA38", "#F8766D", "#619CFF")

#Unweighted UniFrac, color is based on a category from metadata
ord.uf = ordinate(sunfish, "PCoA", "unifrac", weighted=FALSE)
plot_ordination(sunfish, ord.uf, color="IntervalMidpoint") + 
  geom_point(size=5, alpha=0.75) + 
  #scale_color_manual(values = colScale.sample_type) + 
  #scale_shape_manual(values=shape.beach) + 
  ggtitle ("sunfish16S.uwUF.pcoa")
#ggsave("sunfish16S_uwuf_pcoa.pdf")

#with bray-curtis data
ord.bc = ordinate(sunfish, "PCoA", "bray")
plot_ordination(sunfish, ord.bc, color="IntervalMidpoint") + 
  geom_point(size=5, alpha=0.75) + 
  #scale_color_manual(values = colScale.sample_type) + 
  #scale_shape_manual(values=shape.beach) + 
  ggtitle ("sunfish16S.BC.pcoa")
#ggsave("sunfish16S_bc_pcoa.pdf")


#run adonis (aka PERMANOVA), testing distance against sample_type, 999 permutations by default
library(vegan)
df = as(sample_data(sunfish.r), "data.frame")
d = distance(sunfish.r, "unifrac")
sample_type_uf_adonis = adonis(d ~ sample_type, df)
#show results
sample_type_uf_adonis

#adonis sensitive to differences in dispersion of data
#i.e differences in dispersion will results in significant adonis result that is not due to differences in mean/variance

#ideally beta dispersion is NOT significant between groups

#check dispersion of data
betadis <- betadisper(d, df$sample_type)
#show results
betadis

anova(betadis)     
plot(betadis)
boxplot(betadis)
plot(TukeyHSD(betadis))

# extract scores
# use plot_ordination but set justDF = TRUE, to get the dataframe that was used to build the plot.
# also use type = "biplot" to get both sample and species scores combined into one plot
# or can generate separate plots for samples and species using facet panels (type = "split")
# get coordinates of Axis1 and Axis2 (these are the scores), and associated metadata
# to get scores for samples, the column id.type will be Samples
# to get scores for species, or other taxonomic category, column id.type will be Taxa

#unweighted unifrac
sunfish.uf.df <- plot_ordination(sunfish, ord.uf, type = "biplot", justDF = TRUE) 
#bray-curtis
sunfish.bc.df <- plot_ordination(sunfish, ord.bc, type = "biplot", justDF = TRUE) 

#subset unifrac table to just get taxa scores
taxa_scores.uf.df <- sunfish.uf.df[sunfish.uf.df$id.type=="Taxa",]

#get scores that are furthest from 0,0
# calculate distances in the two dimensions (x=Axis.1 and y=Axis.2)
# from an origin of x=0, y=0
# this is Pythagorean theorem (c^2 = a^2 + b^2), so solve for c
taxa_scores.uf.df$distance <- sqrt(taxa_scores.uf.df$Axis.1^2 + taxa_scores.uf.df$Axis.2^2)

#sort for the top 20 distances (i.e. species correlated with explaining greatest variation in samples)
taxa_scores.uf.distsorted <- taxa_scores.uf.df[order(taxa_scores.uf.df$distance, decreasing = TRUE), ]
taxa_scores.uf.top30 <- head(taxa_scores.uf.distsorted, n = 30)

#make plot to add species score arrows
sp_score_plot <- plot_ordination(sunfish, ord.uf, color="IntervalMidpoint") + 
  geom_point(size=5, alpha=0.75) 
  #scale_color_manual(values = colScale.sample_type, name="transect site") + 
  #scale_shape_manual(values=shape.beach)

# Then add species scores as arrows
# use taxa scores to get coordinates for arrows

# Define the arrow aesthetic mapping
# which values will be the coordinates for the end of the arrow 
arrow_map <- aes(xend = Axis.1, yend = Axis.2, x = 0, y = 0, shape = NULL, color = NULL)
# where to draw the label, and what label to use
label_map <- aes(x = 1.2 * Axis.1, y = 1.2 * Axis.2, shape = NULL, color = NULL, label = Phylum)
# size of arrowhead
arrowhead = arrow(length = unit(0.02, "npc"))

# add arrows to plot, using top 30 arrows
sp_score_plot + 
  geom_segment( mapping = arrow_map, linewidth = .5, data = taxa_scores.uf.top30, color = "gray", arrow = arrowhead ) +  
  geom_text( mapping = label_map,  size = 4, data = taxa_scores.uf.top30, show.legend = FALSE) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("sunfish16S.uwUF.pcoa + species scores, phylum labeled")


#ggsave("sunfish16S.uwUF.pcoa_spscores.pdf", useDingbats=FALSE)
#***symbols not interpreted correctly in Adobe Illustrator - a font issue, so make sure to specify useDingbats=FALSE


