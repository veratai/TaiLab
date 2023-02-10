#with new classification using silva132 and sklearn, and no chloroplasts
#plot ordination for bacteria, need to re-do Noriko's capscale analysis
#without log-transformed metadata

library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(vegan)

#function to parse taxonomy, not greengenes format, split by ;
parse_taxonomy_vera <- function (char.vec){
  parse_taxonomy_default(strsplit(char.vec, ";", TRUE)[[1]])
}

setwd("/Users/vera/Projects/BeachMicrobiome/HakaiBeachJune2014/ClassificationQiime2")

#for 16S v4 EMP data, bacteria
#needed to change "taxonomy" column name in otu table to "Consensus Lineage"
#parse taxonomy as default, otherwise will look for greengenes format
#replaced taxonomy with taxonomy using SILVA v132 and scikit-learn

otufile = "./16SV4_all_otutable_4phyloseq_silva132_16SV4_NOCHLOROPLAST_taxonomy.txt"
mapfile = "/Users/vera/Projects/BeachMicrobiome/HakaiBeachJune2014/GoogleDrive/16SV4/201406_16SV4_metadata_map_rmlow.txt"
trefile = "/Users/vera/Projects/BeachMicrobiome/HakaiBeachJune2014/GoogleDrive/16SV4/16SV4_all/16SV4_all_tree.tre"
qiimex = import_qiime(otufile, mapfile, trefile, parseFunction=parse_taxonomy_vera, showProgress = FALSE)
#not sure why, but tree doesn't import correctly

#import separately
tree <- read_tree("/Users/vera/Projects/BeachMicrobiome/HakaiBeachJune2014/GoogleDrive/16SV4/16SV4_all/16SV4_all_tree.tre")

#then create new phyloseq object with tree added, along with other stuff
qiimez_with_zeros <- phyloseq(otu_table(qiimex), sample_data(qiimex), tax_table(qiimex), tree)

#****realized that some OTUs left in with 0 total reads 
#probably OTUs left in after removing samples with less than 8000 total reads
qiimez = prune_taxa(taxa_sums(qiimez_with_zeros) > 0, qiimez_with_zeros)

#rarefy
qiimez.r = rarefy_even_depth(qiimez)

#set shapes and colours manually to match Noriko's plots for prok and protists
shape.beach <- c(8, 1, 17, 18, 16)
colScale.sample_type <- c("#00BA38", "#F8766D", "#619CFF")
#from a ggplot, find info using ggbuild(p)$data

# ***Figure in manuscript, with same colour palette as alpha div plots
ordu.r = ordinate(qiimez.r, "PCoA", "unifrac", weighted=FALSE)
plot_ordination(qiimez.r, ordu.r, color="sample_type", shape="beach") + 
  geom_point(size=5, alpha=0.75) + 
  scale_color_manual(values = c("#009E73", "#E69F00", "#56B4E9")) + 
  scale_shape_manual(values=shape.beach) + 
  theme_set(theme_bw()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle ("16SV4.all.rarefied_uwUF silva132 sklearn, no chloroplasts")
#ggsave("16SV4_all_r_uwuf_pcoa_silva132sklearn_nochloroplasts.pdf")
#***symbols not interpreted correctly in Adobe Illustrator - a font issue, specify useDingbats=FALSE
ggsave("16SV4_all_r_uwuf_pcoa_silva132_nochloroplasts_wlegend.pdf", useDingbats=FALSE)

#*** for manuscript, same size panels, no legend
plot_ordination(qiimez.r, ordu.r, color="sample_type", shape="beach") + 
  geom_point(size=5, alpha=0.75) + 
  scale_color_manual(values = c("#009E73", "#E69F00", "#56B4E9")) + 
  scale_shape_manual(values=shape.beach) + 
  theme_set(theme_bw()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") +
  ggtitle ("16SV4.all.rarefied_uwUF silva132 sklearn, no chloroplasts")
#***symbols not interpreted correctly in Adobe Illustrator - a font issue, specify useDingbats=FALSE
ggsave("16SV4_all_r_uwuf_pcoa_silva132_nochloroplasts_nolegend.pdf", height=5, width=5,useDingbats=FALSE)

#run adonis, testing distance against sample_type, 999 permutations by default
library(vegan)
df = as(sample_data(qiimez.r), "data.frame")
d = distance(qiimez.r, "unifrac")
sample_type_uf_adonis = adonis(d ~ sample_type, df)
sample_type_uf_adonis
#results
#             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#sample_type  2    4.8314 2.41572  9.9937 0.18339  0.001 ***
#Residuals   89   21.5134 0.24172         0.81661           
#Total       91   26.3449                 1.00000             


#adonis sensitive to differences in dispersion of data
#i.e differences in dispersion will results in significant adonis result that is not due to differences in mean/variance

#check dispersion of data
betadis <- betadisper(d, df$sample_type)
betadis
#Average distance to median:
#beach interstitial    beach sand only          sea water 
#           0.5333             0.3944             0.4281 

anova(betadis)     #highly significant
plot(betadis)
boxplot(betadis)
plot(TukeyHSD(betadis))
#dispersion not significant different between beach sand and seawater, but sig different with beach interstitial



#CAP scale analysis, sand only
#on rarefied samples, qiimez.r

#change factors for transect_site
library(plyr)
sample_data(qiimez.r)$transect_site_desc <- revalue(sample_data(qiimez.r)$transect_site_desc, c("high"="high tide", "high_so"="high tide", "low"="low tide", "low_flat"="low tide", "low_so"="low tide", "mid"="mid tide", "mid_bar"="mid tide bar", "mid_sw"="mid tide groundwater", "mid_so"="mid tide", "mid_sw_so"="mid tide groundwater", "gw_seepage"="mid tide seawater", "swash"="swash", "swash_so"="swash"))
#change order of levels for transect_site
sample_data(qiimez.r)$transect_site_desc <- factor(sample_data(qiimez.r)$transect_site_desc, levels = c("high tide", "mid tide", "mid tide bar", "mid tide groundwater", "mid tide seawater", "low tide", "swash"))

##subset samples to do separate analyses for "beach_sand_only"
qiimez.r.SO <- subset_samples(qiimez.r, sample_type == "beach sand only")

# Remove data points with missing metadata
not_na <- subset_samples(qiimez.r.SO, moisture_content != "NA")
not_na <- subset_samples(not_na, SG_Median != "NA")
not_na <- subset_samples(not_na, nh4_n != "NA")
not_na <- subset_samples(not_na, no3_n != "NA")
not_na <- subset_samples(not_na, organic_matter_percent != "NA")
not_na <- subset_samples(not_na, po4_p != "NA")
#no chlorophyll data for West Beach sand only
#not_na <- subset_samples(not_na, chlorophyll != "NA")

uwuf_not_na <- distance(physeq = not_na, method = "UniFrac")

# CAP ordinate
cap_ord <- ordinate(physeq = not_na, method = "CAP", distance = uwuf_not_na, formula = ~ moisture_content + SG_Median + nh4_n  + no3_n + organic_matter_percent + po4_p)
#cap_ord to see results

#statistical test, 999 permutations by default
anova(cap_ord)
#these are correct:
#Model: capscale(formula = distance ~ moisture_content + SG_Median + nh4_n + no3_n + organic_matter_percent + po4_p, data = data)
#Df SumOfSqs      F Pr(>F)    
#Model     6   1.2143 1.3597  0.001 ***
#  Residual 18   2.6792                                 

# CAP plot
#set shapes and colours manually to match other plots
shape.beach <- c(8,16)
#colScale.sample_type <- c("#00BA38", "#F8766D", "#619CFF")

cap_plot <- plot_ordination(physeq = not_na, ordination = cap_ord, color = "transect_site_desc", axes = c(1,2)) +
  aes(shape = beach) + 
  geom_point(size=4, alpha=0.75) + 
  scale_colour_manual(values=c("#009E73","#CC79A7","#D55E00","#56B4E9","#0072B2"), name="transect site") + 
  #scale_color_manual(values = colScale.sample_type) + 
  #scale_color_manual(values = colScale.sample_type, breaks=c("high", "gw_seepage", "mid_bar", "mid", "low", "swash", "seawater"), labels=c("high", "gw_seepage", "waterbar", "mid", "low", "swash", "seawater")) + 
  scale_shape_manual(values=shape.beach)

# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")
# Add labels, make a data.frame
arrowdf <- data.frame(labels = c("moisture", "grain size", "NH4", "NO3", "orgC", "PO4"), arrowmat)
# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, yend = CAP2, x = 0, y = 0, shape = NULL, color = NULL)
label_map <- aes(x = 1.2 * CAP1, y = 1.2 * CAP2, shape = NULL, color = NULL, label = labels)
arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
cap_plot + 
  geom_segment( mapping = arrow_map, size = .5, data = arrowdf, color = "gray", arrow = arrowhead ) +  
  geom_text( mapping = label_map,  size = 4, data = arrowdf, show.legend = FALSE) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("16SV4.r.CAP_sandonly.uwUniFrac silva132 sklearn, no chloroplasts")
#***symbols not interpreted correctly in Adobe Illustrator - a font issue, specify useDingbats=FALSE
ggsave("16SV4_r_CAP_sandonly.uwUniFrac_silva132sklearn_nochloroplasts_wlegend.pdf", useDingbats=FALSE)

# for manuscript, same size panels, nolegend
cap_plot + 
  geom_segment( mapping = arrow_map, size = .5, data = arrowdf, color = "gray", arrow = arrowhead ) +  
  geom_text( mapping = label_map,  size = 4, data = arrowdf, show.legend = FALSE) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") +
  ggtitle("16SV4.r.CAP_sandonly.uwUniFrac silva132 sklearn, no chloroplasts")
#***symbols not interpreted correctly in Adobe Illustrator - a font issue, specify useDingbats=FALSE
ggsave("16SV4_r_CAP_sandonly.uwUniFrac_silva132sklearn_nochloroplasts_nolegend.pdf", height=5, width=5, useDingbats=FALSE)


#run adonis, testing distance against beach, 999 permutations by default
df = as(sample_data(qiimez.r.SO), "data.frame")
d = distance(qiimez.r.SO, "unifrac")
sandonly_uf_adonis = adonis(d ~ beach, df)
sandonly_uf_adonis
#results
#          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#beach      1    0.3355 0.33553  2.1672 0.08282  0.001 ***
#Residuals 24    3.7158 0.15482         0.91718           
#Total     25    4.0513                 1.00000   
#adonis sensitive to differences in dispersion of data
#i.e differences in dispersion will results in significant adonis result that is not due to differences in mean/variance

#check dispersion of data
betadis <- betadisper(d, df$beach)
betadis
#Average distance to median:
#2nd   West 
#0.3750 0.3797 

anova(betadis)
plot(betadis)
boxplot(betadis)
plot(TukeyHSD(betadis))
#dispersion not significantly different between 2nd and West beaches

#run adonis, testing distance against transect_site_desc, 999 permutations by default
df = as(sample_data(qiimez.r.SO), "data.frame")
d = distance(qiimez.r.SO, "unifrac")
sandonly_uf_adonis = adonis(d ~ transect_site_desc, df)
sandonly_uf_adonis
#results
#                   Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#transect_site_desc  4    0.8458 0.21145  1.3853 0.20877  0.001 ***
#Residuals          21    3.2055 0.15264         0.79123           
#Total              25    4.0513                 1.00000
#adonis sensitive to differences in dispersion of data
#i.e differences in dispersion will results in significant adonis result that is not due to differences in mean/variance

#check dispersion of data
betadis <- betadisper(d, df$transect_site_desc)
betadis
#Average distance to median:
#high tide          mid tide mid tide seawater          low tide 
#0.3534            0.3555            0.2955            0.3575 
#swash 
#0.3630 

anova(betadis)
plot(betadis)
boxplot(betadis)
plot(TukeyHSD(betadis))
#dispersion significantly different with mid-tide seawater



##subset samples to do separate analyses for "interstitial"
qiimez.r.int <- subset_samples(qiimez.r, sample_type == "beach interstitial")

# Remove data points with missing metadata
not_na <- subset_samples(qiimez.r.int, moisture_content != "NA")
not_na <- subset_samples(not_na, SG_Median != "NA")
not_na <- subset_samples(not_na, nh4_n != "NA")
not_na <- subset_samples(not_na, no3_n != "NA")
not_na <- subset_samples(not_na, organic_matter_percent != "NA")
not_na <- subset_samples(not_na, po4_p != "NA")
#considering chlorophyll
not_na <- subset_samples(not_na, chlorophyll != "NA")

uwuf_not_na <- distance(physeq = not_na, method = "UniFrac")

# CAP ordinate, adding in chlorophyll for this analysis
cap_ord <- ordinate(physeq = not_na, method = "CAP", distance = uwuf_not_na, formula = ~ moisture_content + SG_Median + nh4_n  + no3_n + organic_matter_percent + po4_p + chlorophyll)
#cap_ord to see results

#statistical test, 999 permutations by default
anova(cap_ord)
#these are correct:
#Model: capscale(formula = distance ~ moisture_content + SG_Median + nh4_n + no3_n + organic_matter_percent + po4_p + chlorophyll, data = data)
#Df SumOfSqs      F Pr(>F)    
#Model     7   3.1465 1.7554  0.001 ***
#Residual 38   9.7308                   

# CAP plot
#set shapes and colours manually to match other plots
shape.beach <- c(8, 1, 17, 18, 16)
#colScale.sample_type <- c("#00BA38", "#F8766D", "#619CFF")

#updating colours, transect_site desc for manuscript
cap_plot <- plot_ordination(physeq = not_na, ordination = cap_ord, color = "transect_site_desc", axes = c(1,2)) +
  aes(shape = beach) + 
  geom_point(size=4, alpha=0.75) + 
  scale_colour_manual(values=c("#009E73","#CC79A7","#FDBF6F", "#D55E00","#984EA3", "#56B4E9","#0072B2"), name="transect site") + 
  #scale_color_manual(values = colScale.sample_type) + 
  #scale_color_manual(values = colScale.sample_type, breaks=c("high", "gw_seepage", "mid_bar", "mid", "low", "swash", "seawater"), labels=c("high", "gw_seepage", "waterbar", "mid", "low", "swash", "seawater")) + 
  scale_shape_manual(values=shape.beach)

# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")
# Add labels, make a data.frame
arrowdf <- data.frame(labels = c("moisture", "grain size", "NH4", "NO3", "orgC", "PO4", "chl"), arrowmat)
# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, yend = CAP2, x = 0, y = 0, shape = NULL, color = NULL)
label_map <- aes(x = 1.5 * CAP1, y = 1.0 * CAP2, shape = NULL, color = NULL, label = labels)
arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
cap_plot + 
  geom_segment( mapping = arrow_map, size = .5, data = arrowdf, color = "gray", arrow = arrowhead ) +  
  geom_text( mapping = label_map,  size = 4, data = arrowdf, show.legend = FALSE) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("16SV4.r.CAP_interstitial.uwUniFrac silva132 sklearn, no chloroplasts")
#***symbols not interpreted correctly in Adobe Illustrator - a font issue, specify useDingbats=FALSE
ggsave("16SV4_r_CAP_interstitial.uwUniFrac_silva132sklearn_nochloroplasts_wlegend.pdf", useDingbats=FALSE)

# for manuscript, same size panels, no legend
cap_plot + 
  geom_segment( mapping = arrow_map, size = .5, data = arrowdf, color = "gray", arrow = arrowhead ) +  
  geom_text( mapping = label_map,  size = 4, data = arrowdf, show.legend = FALSE) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") +
  ggtitle("16SV4.r.CAP_interstitial.uwUniFrac silva132 sklearn, no chloroplasts")
#***symbols not interpreted correctly in Adobe Illustrator - a font issue, specify useDingbats=FALSE
ggsave("16SV4_r_CAP_interstitial.uwUniFrac_silva132sklearn_nochloroplasts_nolegend.pdf", height=5, width=5, useDingbats=FALSE)

#run adonis, testing distance against beach, 999 permutations by default
df = as(sample_data(qiimez.r.int), "data.frame")
d = distance(qiimez.r.int, "unifrac")
interstitial_uf_adonis = adonis(d ~ beach, df)
interstitial_uf_adonis
#results
#          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#beach      4    2.4222 0.60554  2.2681 0.16473  0.001 ***
#Residuals 46   12.2812 0.26698         0.83527           
#Total     50   14.7034                 1.00000  

#adonis sensitive to differences in dispersion of data
#i.e differences in dispersion will results in significant adonis result that is not due to differences in mean/variance

#check dispersion of data
betadis <- betadisper(d, df$beach)
betadis
#Average distance to median:
#  2nd    3rd    7th  North   West 
#0.4980 0.4649 0.4907 0.4827 0.5114 

anova(betadis)
plot(betadis)
boxplot(betadis)
plot(TukeyHSD(betadis))
#dispersion not significantly different between beaches

#run adonis, testing distance against transect_site_desc, 999 permutations by default
df = as(sample_data(qiimez.r.int), "data.frame")
d = distance(qiimez.r.int, "unifrac")
interstitial_uf_adonis = adonis(d ~ transect_site_desc, df)
interstitial_uf_adonis
#results
#                   Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#transect_site_desc  6    2.4987 0.41645  1.5014 0.16994  0.001 ***
#Residuals          44   12.2047 0.27738         0.83006           
#Total              50   14.7034                 1.00000                 

#adonis sensitive to differences in dispersion of data
#i.e differences in dispersion will results in significant adonis result that is not due to differences in mean/variance

#check dispersion of data
betadis <- betadisper(d, df$transect_site_desc)
betadis
#Average distance to median:
#           high tide             mid tide         mid tide bar mid tide groundwater 
#              0.4995               0.4827               0.0000               0.4187 
#mid tide seawater             low tide                swash 
#           0.4131               0.4899               0.5178 

anova(betadis)
plot(betadis)
boxplot(betadis)
plot(TukeyHSD(betadis))
#mid_bar, one sample, so 0 dispersion
#otherwise not significantly different between transect_sites


###checking results

#sand only PCoA, with rarefied data
# =qiimez.r.SO
# ***Figure in manuscript, with same colour palette as alpha div plots
ordu.r.SO = ordinate(qiimez.r.SO, "PCoA", "unifrac", weighted=FALSE)
plot_ordination(qiimez.r.SO, ordu.r.SO, color="transect_site_desc", shape="beach") + 
  geom_point(size=5, alpha=0.75) + 
  scale_colour_manual(values=c("#009E73","#CC79A7","#D55E00","#56B4E9","#0072B2"), name="transect site") + 
  scale_shape_manual(values=shape.beach) + 
  theme_set(theme_bw()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle ("16SV4.all.rarefied_uwUF silva132 sklearn, no chloroplasts SAND ONLY")
#***symbols not interpreted correctly in Adobe Illustrator - a font issue, specify useDingbats=FALSE
#ggsave("16SV4_all_r_uwuf_pcoa_silva132_nochloroplasts_SANDONLY.pdf", useDingbats=FALSE)

#examine frequency of OTUs
#http://evomics.org/wp-content/uploads/2016/01/phyloseq-Lab-01-Answers.html

#sum counts of OTUs (total reads, not rarefied)
library("data.table")
qiimeztotal_dt = data.table(tax_table(qiimez),
                 TotalCounts = taxa_sums(qiimez),
                 OTU = taxa_names(qiimez))
ggplot(qiimeztotal_dt, aes(TotalCounts)) + 
  geom_histogram(binwidth=100) + 
  xlim(0,2500) +
  ylim(0,10000)
  ggtitle("Histogram of Total Counts")


#examine cumulative sum to determine filtering threshold
# assign Total counts to a new data table
taxcumsum = qiimeztotal_dt[, .N, by = TotalCounts]
#sort data table by cumulative sums, and set key
setkey(taxcumsum, TotalCounts)
#add up OTU counts (cumsum(N)), and make new variable CumSum
taxcumsum[, CumSum := cumsum(N)]
# Define the plot
pCumSum = ggplot(taxcumsum, aes(TotalCounts, CumSum)) + 
    geom_point() +
    xlab("Total Counts (num reads") +
    ylab("Cumulative Sum of taxa") +
    ggtitle("Cummulative sum of taxa vs. Total counts (reads")
pCumSum
# a small number of reads are very abundant (high on x-axis)
# high on y-axis are the OTUs where 1 OTU has high number of reads
# y axis makes big jumps low on axis
# because lots of OTUs have low number of reads (no progress on x-axis)

# Zoom-in
pCumSum + xlim(0, 100)


#filter data, keep only most abundant OTUs
#from rarified data, get relative abundance of OTUs by sample
qiimez.r.SO_relA <- transform_sample_counts(qiimez.r.SO, function(x) x / sum(x))
#examine frequency distribution


#remove OTUs that have an overall rel A < 0.01 (1%) across samples
qiimez_relA_1perc = filter_taxa(qiimez_relA, function(x) x > 0.01, TRUE)




##check results when also removing "Mitochondria" OTUs
#for 16S v4 EMP data, bacteria
#needed to change "taxonomy" column name in otu table to "Consensus Lineage"
#parse taxonomy as default, otherwise will look for greengenes format
#replaced taxonomy with taxonomy using SILVA v132 and scikit-learn

otufile_m = "./16SV4_all_otutable_4phyloseq_silva132_16SV4_NOCHLORO_NOMITO_taxonomy.txt"
mapfile = "/Users/vera/Projects/BeachMicrobiome/HakaiBeachJune2014/GoogleDrive/16SV4/201406_16SV4_metadata_map_rmlow.txt"
trefile = "/Users/vera/Projects/BeachMicrobiome/HakaiBeachJune2014/GoogleDrive/16SV4/16SV4_all/16SV4_all_tree.tre"
qiimex_m = import_qiime(otufile_m, mapfile, trefile, parseFunction=parse_taxonomy_vera, showProgress = FALSE)

#import separately
tree <- read_tree("/Users/vera/Projects/BeachMicrobiome/HakaiBeachJune2014/GoogleDrive/16SV4/16SV4_all/16SV4_all_tree.tre")

#then create new phyloseq object with tree added, along with other stuff
qiimez_m <- phyloseq(otu_table(qiimex_m), sample_data(qiimex_m), tax_table(qiimex_m), tree)

#then create new phyloseq object with tree added, along with other stuff
qiimez_m_with_zeros <- phyloseq(otu_table(qiimex_m), sample_data(qiimex_m), tax_table(qiimex_m), tree)

#****realized that some OTUs left in with 0 total reads 
#probably OTUs left in after removing samples with less than 8000 total reads
qiimez_m = prune_taxa(taxa_sums(qiimez_m_with_zeros) > 0, qiimez_m_with_zeros)

### raw OTU abundances per sample
total_reads <- sample_sums(qiimez_m)
sum(total_reads)
summary(total_reads)

### raw OTU abundances per sample
total_reads <- sample_sums(qiimez_m)
sum(total_reads)
summary(total_reads)

#rarefy
qiimez_m.r = rarefy_even_depth(qiimez_m)

# check plot
# ***Figure in manuscript, with same colour palette as alpha div plots
ordu_m.r = ordinate(qiimez_m.r, "PCoA", "unifrac", weighted=FALSE)
plot_ordination(qiimez_m.r, ordu_m.r, color="sample_type", shape="beach") + 
  geom_point(size=5, alpha=0.75) + 
  scale_color_manual(values = c("#009E73", "#E69F00", "#56B4E9")) + 
  scale_shape_manual(values=shape.beach) + 
  theme_set(theme_bw()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle ("16SV4.all.rarefied_uwUF silva132 sklearn, no chloroplasts no mitochondria")


#change factors for transect_site
library(plyr)
sample_data(qiimez_m.r)$transect_site_desc <- revalue(sample_data(qiimez_m.r)$transect_site_desc, c("high"="high tide", "high_so"="high tide", "low"="low tide", "low_flat"="low tide", "low_so"="low tide", "mid"="mid tide", "mid_bar"="mid tide bar", "mid_sw"="mid tide groundwater", "mid_so"="mid tide", "mid_sw_so"="mid tide groundwater", "gw_seepage"="mid tide seawater", "swash"="swash", "swash_so"="swash"))
#change order of levels for transect_site
sample_data(qiimez_m.r)$transect_site_desc <- factor(sample_data(qiimez_m.r)$transect_site_desc, levels = c("high tide", "mid tide", "mid tide bar", "mid tide groundwater", "mid tide seawater", "low tide", "swash"))

##subset samples to do separate analyses for "beach_sand_only"
qiimez_m.r.SO <- subset_samples(qiimez_m.r, sample_type == "beach sand only")

# Remove data points with missing metadata
not_na <- subset_samples(qiimez_m.r.SO, moisture_content != "NA")
not_na <- subset_samples(not_na, SG_Median != "NA")
not_na <- subset_samples(not_na, nh4_n != "NA")
not_na <- subset_samples(not_na, no3_n != "NA")
not_na <- subset_samples(not_na, organic_matter_percent != "NA")
not_na <- subset_samples(not_na, po4_p != "NA")
#no chlorophyll data for West Beach sand only
#not_na <- subset_samples(not_na, chlorophyll != "NA")

uwuf_not_na <- distance(physeq = not_na, method = "UniFrac")

# CAP ordinate
cap_ord <- ordinate(physeq = not_na, method = "CAP", distance = uwuf_not_na, formula = ~ moisture_content + SG_Median + nh4_n  + no3_n + organic_matter_percent + po4_p)
#cap_ord to see results

# CAP plot
#set shapes and colours manually to match other plots
shape.beach <- c(8,16)
#colScale.sample_type <- c("#00BA38", "#F8766D", "#619CFF")

cap_plot <- plot_ordination(physeq = not_na, ordination = cap_ord, color = "transect_site_desc", axes = c(1,2)) +
  aes(shape = beach) + 
  geom_point(size=4, alpha=0.75) + 
  scale_colour_manual(values=c("#009E73","#CC79A7","#D55E00","#56B4E9","#0072B2"), name="transect site") + 
  #scale_color_manual(values = colScale.sample_type) + 
  #scale_color_manual(values = colScale.sample_type, breaks=c("high", "gw_seepage", "mid_bar", "mid", "low", "swash", "seawater"), labels=c("high", "gw_seepage", "waterbar", "mid", "low", "swash", "seawater")) + 
  scale_shape_manual(values=shape.beach)

# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")
# Add labels, make a data.frame
arrowdf <- data.frame(labels = c("moisture", "grain size", "NH4", "NO3", "orgC", "PO4"), arrowmat)
# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, yend = CAP2, x = 0, y = 0, shape = NULL, color = NULL)
label_map <- aes(x = 1.2 * CAP1, y = 1.2 * CAP2, shape = NULL, color = NULL, label = labels)
arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
cap_plot + 
  geom_segment( mapping = arrow_map, size = .5, data = arrowdf, color = "gray", arrow = arrowhead ) +  
  geom_text( mapping = label_map,  size = 4, data = arrowdf, show.legend = FALSE) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("16SV4.r.CAP_sandonly.uwUniFrac silva132 sklearn, no chloroplasts no mitochondria")

##subset samples to do separate analyses for "interstitial"
qiimez_m.r.int <- subset_samples(qiimez_m.r, sample_type == "beach interstitial")

# Remove data points with missing metadata
not_na <- subset_samples(qiimez_m.r.int, moisture_content != "NA")
not_na <- subset_samples(not_na, SG_Median != "NA")
not_na <- subset_samples(not_na, nh4_n != "NA")
not_na <- subset_samples(not_na, no3_n != "NA")
not_na <- subset_samples(not_na, organic_matter_percent != "NA")
not_na <- subset_samples(not_na, po4_p != "NA")
#considering chlorophyll
not_na <- subset_samples(not_na, chlorophyll != "NA")

uwuf_not_na <- distance(physeq = not_na, method = "UniFrac")

# CAP ordinate, adding in chlorophyll for this analysis
cap_ord <- ordinate(physeq = not_na, method = "CAP", distance = uwuf_not_na, formula = ~ moisture_content + SG_Median + nh4_n  + no3_n + organic_matter_percent + po4_p + chlorophyll)
#cap_ord to see results

# CAP plot
#set shapes and colours manually to match other plots
shape.beach <- c(8, 1, 17, 18, 16)
#colScale.sample_type <- c("#00BA38", "#F8766D", "#619CFF")

#updating colours, transect_site desc for manuscript
cap_plot <- plot_ordination(physeq = not_na, ordination = cap_ord, color = "transect_site_desc", axes = c(1,2)) +
  aes(shape = beach) + 
  geom_point(size=4, alpha=0.75) + 
  scale_colour_manual(values=c("#009E73","#CC79A7","#FDBF6F", "#D55E00","#984EA3", "#56B4E9","#0072B2"), name="transect site") + 
  #scale_color_manual(values = colScale.sample_type) + 
  #scale_color_manual(values = colScale.sample_type, breaks=c("high", "gw_seepage", "mid_bar", "mid", "low", "swash", "seawater"), labels=c("high", "gw_seepage", "waterbar", "mid", "low", "swash", "seawater")) + 
  scale_shape_manual(values=shape.beach)

# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")
# Add labels, make a data.frame
arrowdf <- data.frame(labels = c("moisture", "grain size", "NH4", "NO3", "orgC", "PO4", "chl"), arrowmat)
# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, yend = CAP2, x = 0, y = 0, shape = NULL, color = NULL)
label_map <- aes(x = 1.5 * CAP1, y = 1.0 * CAP2, shape = NULL, color = NULL, label = labels)
arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
cap_plot + 
  geom_segment( mapping = arrow_map, size = .5, data = arrowdf, color = "gray", arrow = arrowhead ) +  
  geom_text( mapping = label_map,  size = 4, data = arrowdf, show.legend = FALSE) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("16SV4.r.CAP_interstitial.uwUniFrac silva132 sklearn, no chloroplasts no mitochondria")
