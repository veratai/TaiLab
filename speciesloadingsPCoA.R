# Add species loadings arrows to PCoA


library(phyloseq)
library(ggplot2)
library(vegan)
library(tidyverse)
library(tidyr)
library(dplyr)


# Rarefy data
phy_site1.exp.r <- rarefy_even_depth(phy_site1_exp_format, rngseed = TRUE)

# PCoA of bray-curtis distances with rarified data
ordu.site1.exp.r = ordinate(phy_site1.exp.r, "PCoA", "bray", weighted=FALSE)

#create table of ordinating results using phyloseq object used to generate PCoA (phy) and the resulting ordination (ordu) 
site.123.exp.df <- plot_ordination(phy_site123.exp.r, ordu.site123.exp.r, type = "biplot", justDF = TRUE) 

#subset table to just get taxa scores
taxa_scores.site.123.exp.df <- site.123.exp.df[site.123.exp.df$id.type=="Taxa",]

#get scores that are furthest from 0,0
# calculate distances in the two dimensions (x=Axis.1 and y=Axis.2)
# from an origin of x=0, y=0
# this is Pythagorean theorem (c^2 = a^2 + b^2), so solve for c
taxa_scores.site.123.exp.df$distance <- sqrt(taxa_scores.site.123.exp.df$Axis.1^2 + taxa_scores.site.123.exp.df$Axis.2^2)

#sort based on distances
taxa_scores.site.123.exp.distsorted <- taxa_scores.site.123.exp.df[order(taxa_scores.site.123.exp.df$distance, decreasing = TRUE), ]
#get top 20, i.e longest distances, species correlated with explaining greatest variation in samples 
taxa_scores.site.123.exp.distsorted_top20 <- taxa_scores.site.123.exp.distsorted[1:20, ]


#make plot to add species score arrows
#first plot ordination
sp_score_plot <- plot_ordination(phy_site123.exp.r, ordu.site123.exp.r, color="season") + 
  geom_point(size=10, alpha=0.75) 

# Then add species scores as arrows
# use taxa scores to get coordinates for arrows

# Define the arrow aesthetic mapping
# which values will be the coordinates for the end of the arrow 
arrow_map <- aes(xend = Axis.1, yend = Axis.2, x = 0, y = 0, shape = NULL, color = NULL)

# where to draw the label, and what label to use
label_map <- aes(x = 1.2 * Axis.1, y = 1.2 * Axis.2, shape = NULL, color = NULL, label = best_hit)

# size of arrowhead
arrowhead = arrow(length = unit(0.02, "npc"))

# add arrows to plot,looking at all arrows
Site123ExpSpeciesLoading <- sp_score_plot + 
  geom_segment( mapping = arrow_map, linewidth = .5, data = taxa_scores.site.123.exp.distsorted_top20, color = "gray", arrow = arrowhead ) +  
  #geom_text( mapping = label_map,  size = 2, data = taxa_scores.site.123.exp.distsorted_top20, show.legend = FALSE) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=60), axis.title.y = element_text(size=60), 
        axis.text = element_text(size=60), legend.text = element_text(size=60),
        legend.title = element_text(size=60)) +
  ggtitle("Site 123 PCoA with top 20 species loadings")+
  scale_colour_manual(values=c("#fee5d9", "#fcbba1", "#fc9272", "#fb6a4a", "#ef3b2c", "#cb181d", "#a50f15",
                                        "#67000d", "#6baed6", "#4292c6", "#2171b5", "#08519c", "#08306b"))
                                        
ggsave("Site123SpeciesLoading_top20.pdf", plot = Site123ExpSpeciesLoading, width = 30, height = 15,
                                               dpi = 300, path="~/ThesisFinalCodeR/R_Home/FinalThesisFigures")






#CUSTOMIZING arrows based on angle, direction

#get the angle for each arrow
taxa_scores.site.123.exp.distsorted.angle <-  angle(taxa_scores.site.123.exp.distsorted,
                                                 x_col = "Axis.1",
                                                 y_col = "Axis.2", origin = c(0,0))

#filter the arrows based on angle of the arrow
#customize this filtering based on PCoA biplot
#in this case, excluding arrows between 50 to 170 degrees, and 225 to 315 degrees
#so keeping arrows between 315 and 50 degrees, and 170 to 225.
filtered_angle <- taxa_scores.site.123.exp.distsorted.angle[!((taxa_scores.site.123.exp.distsorted.angle$.degrees >= 50 &
                                                                 taxa_scores.site.123.exp.distsorted.angle$.degrees <= 170) |
                                                              (taxa_scores.site.123.exp.distsorted.angle$.degrees >= 225 &
                                                                 taxa_scores.site.123.exp.distsorted.angle$.degrees <= 315)), ]

#label arrows whether they are "positive" or "negative" on Axis 1
filtered_angle <- filtered_angle %>% 
  mutate(Direction = ifelse(Axis.1 > 0, "Positive Axis1", "Negative Axis1"))

#Subset filtered angle dataframe into two dataframes with only positive and negative axis1 values
filtered_angle_positive <- filtered_angle[filtered_angle$Direction == "Positive Axis1", ] 
filtered_angle_negative <- filtered_angle[filtered_angle$Direction == "Negative Axis1", ] 

# Site 123 filtered angle dataframe
filtered_angle_positive <- filtered_angle_positive[ -c(12:22,25)]
filtered_angle_negative <- filtered_angle_negative[ -c(12:22,25)]

threshold1 <- quantile(filtered_angle_positive$distance, 0.90)
threshold2 <- quantile(filtered_angle_negative$distance, 0.90)
#filtered_angle_top10 <- filtered_angle[filtered_angle$distance >= threshold, ]
filtered_angle_positive_top10 <- filtered_angle_positive[filtered_angle_positive$distance >= threshold1, ]
filtered_angle_negative_top10 <- filtered_angle_negative[filtered_angle_negative$distance >= threshold2, ]



#add arrows to plot, looking at the top ten percent of from the subsetted specific arrows
Site123ExpSpeciesLoadingTop10 <- sp_score_plot + 
  geom_segment( mapping = arrow_map, linewidth = .5, data = filtered_angle_top10, color = "gray", arrow = arrowhead ) +  
  #geom_text( mapping = label_map,  size = 2, data = taxa_scores.site.1.top30, show.legend = FALSE) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=60), axis.title.y = element_text(size=60), 
        axis.text = element_text(size=60), legend.text = element_text(size=60),
        legend.title = element_text(size=60)) +
  ggtitle("Site 123 exp taxa distribution")+
  scale_colour_manual(values=c("#fee5d9", "#fcbba1", "#fc9272", "#fb6a4a", "#ef3b2c", "#cb181d", "#a50f15",
                                        "#67000d", "#6baed6", "#4292c6", "#2171b5", "#08519c", "#08306b"))
                                        
ggsave("Site123SpeciesLoadingExpTop10.pdf", plot = Site123ExpSpeciesLoadingTop10 + 
         scale_colour_manual(values=c("#fee5d9", "#fcbba1", "#fc9272", "#fb6a4a", "#ef3b2c", "#cb181d",
                                               "#a50f15", "#67000d", "#6baed6", "#4292c6", "#2171b5", "#08519c", "#08306b")), width = 30, height = 15,
                                               dpi = 300, path="~/ThesisFinalCodeR/R_Home/FinalThesisFigures")


# add arrows that are only at the top 10 positive and negative axis 1 only
Site123ExpSpeciesLoadingPositiveNegativeTop10 <- sp_score_plot + 
  geom_segment( mapping = arrow_map, linewidth = .5, data = filtered_angle_positive_negative_top10sum, color = "gray", arrow = arrowhead ) +  
  #geom_text( mapping = label_map,  size = 2, data = taxa_scores.site.1.top30, show.legend = FALSE) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=60), axis.title.y = element_text(size=60), 
        axis.text = element_text(size=60), legend.text = element_text(size=60),
        legend.title = element_text(size=60)) +
  ggtitle("Site 123 exp taxa distribution")+
  scale_colour_manual(values=c("#fee5d9", "#fcbba1", "#fc9272", "#fb6a4a", "#ef3b2c", "#cb181d", "#a50f15",
                                        "#67000d", "#6baed6", "#4292c6", "#2171b5", "#08519c", "#08306b"))
                                        
ggsave("Site123SpeciesLoadingPositiveandNegativeTop10.pdf", plot = Site123ExpSpeciesLoadingPositiveNegativeTop10 + 
         scale_colour_manual(values=c("#fee5d9", "#fcbba1", "#fc9272", "#fb6a4a", "#ef3b2c", "#cb181d",
                                               "#a50f15", "#67000d", "#6baed6", "#4292c6", "#2171b5", "#08519c", "#08306b")), width = 30, height = 15,
                                               dpi = 300, path="~/ThesisFinalCodeR/R_Home/FinalThesisFigures")
                                               
                                               
                                               
#ANOTHER OPTION (Gavin Simpson, https://fromthebottomoftheheap.net/slides/intro-vegan-webinar-2020/intro-to-vegan.html#1):
#PCoA adding species scores
pco <- wcmdscale(vegdist(varespec), eig = TRUE)
## plot
plot(pco)
## get PCoA scores
scrs <- scores(pco, choices = 1:2)
## take WA of PCoA scores,
## weight by abundance
spp_scrs <- wascores(scrs, varespec,
                     expand = FALSE)
## add
points(spp_scrs, col = "red", pch = 19)

