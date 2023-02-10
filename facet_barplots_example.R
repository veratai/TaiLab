#data is imported as a phyloseq object called qiimez

#list rank_names used in taxonomy
rank_names(qiimez)
#list unique taxa names in Rank2
get_taxa_unique(qiimez, "Rank2")

#rename taxa names, if needed
#for Rank 2 replace "Eukaryota_X" with "Eukaryota"
tax_table(qiimez)[,2][tax_table(qiimez)[,2]=="Eukaryota_X"] <- "Eukaryota"

#glom Rank2 OTUs together, so will not show lines for individuals OTUs
rank2glommed = tax_glom(qiimez, "Rank2")

#get relative abundance
rank2glommed.relA <- transform_sample_counts(rank2glommed, function(x) x / sum(x))

#melt data and turn into dataframe to control levels
rank2glommed.relA.df <- psmelt(rank2glommed.relA)

#get average relA for each beach_site
rank2_beach_site_mean <- lapply( split(rank2glommed.relA.df$Abundance, f=list(rank2glommed.relA.df$Rank2, rank2glommed.relA.df$beach_site)), mean )
#make data frame for these data
rank2_beach_site_relA <- as.data.frame(unlist(rank2_beach_site_mean))
colnames(rank2_beach_site_relA) <- "meanRelA"

#get standard deviation
rank2_beach_site_sd <- lapply( split(rank2glommed.relA.df$Abundance, f=list(rank2glommed.relA.df$Rank2, rank2glommed.relA.df$beach_site)), sd )
#add sd to dataframe
rank2_beach_site_relA$sdRelA <- unlist(rank2_beach_site_sd)


#need to recover info of beach_site and Rank2 names that have been merged together in rownames of the dataframe
split_names <- strsplit(rownames(rank2_beach_site_relA), "\\.")
rank2 <- sapply(split_names, function(x) x[1])
beach_site <- sapply(split_names, function(x) x[2])

#add to datatable, order levels
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

#get transect_site info that is in the text of beach_site info
transect_site <- lapply(beach_site, function(x) {sub(("[[:alnum:]]+_"), "", x) } )
transect_site <- unlist(lapply(transect_site, function(x) {sub(("seawater"), "sw", x) } )  )
#order and add the transect_site names to dataframe
rank2_beach_site_relA$transect_site <- factor(transect_site, levels=c("high_so", "mid_so", "mid_gw_so", "low_so", "swash_so", "high", "mid", "mid_bar", "mid_gw", "mid_sw", "low", "swash", "sw"))


#get range of a specific colour, get min/max of colours from ColorBrewer "Paired" palette, eg. brewer.pal(n = 4, name = "Paired")
#greenRamp <- colorRampPalette(c("#00441B","#C7E9C0"))
#greenPalette <- greenRamp(5)
greenPalette <- c("#006837", "#1A9850", "#66BD63", "#A6D96A", "#D9EF8B") 
#orangeRamp <- colorRampPalette(c("#7F0000", "#FDD49E"))
#orangePalette <- orangeRamp(7)
orangePalette <- c("#7F2704", "#A63603", "#D94801", "#F16913", "#FD8D3C", "#FDAE6B", "#FDD0A2")
blue <- "#56B4E9"
# get palette of 8 colours, 5 green shades for sand samples, 7 orange shades for interstial, blue for seawater
transectPalette <- c(greenPalette, orangePalette, blue)


#remove taxa with no data, i.e. samples with no reads, i.e. samples that didn't work
#otherwise get column of NA data, not sure why this is happening this time.
rank2_beach_site_relA_noNA <- rank2_beach_site_relA[rank2_beach_site_relA$meanRelA != 'NaN',]


#plot, coloured by transect_site
p <- ggplot(rank2_beach_site_relA_noNA, aes(x=beach_site, y=meanRelA, fill=transect_site, order = as.factor(rank2)))

p + 
  geom_bar(stat="identity", colour="black") +
  scale_fill_manual(values=transectPalette, labels=c("high sand", "mid sand", "mid groundwater sand", "low sand", "swash sand", "high interstitial", "mid interstitial", "mid bar interstitial", "mid groundwater interstitial", "mid seawater interstitial", "low interstitial", "swash interstitial", "seawater")) +   
  #coord_cartesian(ylim = c(0, 0.05)) +
  ylab("relative abundance") +
  theme_set(theme_bw()) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=6), legend.title=element_blank(), axis.title.x=element_blank(), legend.position="bottom") +
  ggtitle("18SV4 protists, relative abundance") +
  facet_wrap(~rank2, scales="free_y")
ggsave("./plots_final/18SV4_protists_relA_rank2_barplot_final.pdf", width=10, height=5, units="in")

#plot with points +/- sd, and coloured by transect_site
p <- ggplot(rank2_beach_site_relA_noNA, aes(x=beach_site, y=meanRelA, colour=transect_site, order = as.factor(rank2)))

p + 
  geom_point(alpha=0.85) +
  geom_errorbar(aes(ymin=meanRelA-sdRelA, ymax=meanRelA+sdRelA), width=.1) +
  scale_colour_manual(values=transectPalette, labels=c("high sand", "mid sand", "mid groundwater sand", "low sand", "swash sand", "high interstitial", "mid interstitial", "mid bar interstitial", "mid groundwater interstitial", "mid seawater interstitial", "low interstitial", "swash interstitial", "seawater")) +   
  #coord_cartesian(ylim = c(0, 0.05)) +
  #scale_shape_manual(values=c(15,17,3,19,18,0,2,10,3,7,1,8,5,4)) +  #by default only 6 shapes used, need to add more and specifically define shape type
  ylab("relative abundance") +
  theme_set(theme_bw()) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5), legend.title=element_blank(), axis.title.x=element_blank(), legend.position="bottom") +
  ggtitle("18SV4 protists, relative abundance") +
  facet_wrap(~rank2, scales="free_y")
ggsave("./plots_final/18SV4_protists_relA_rank2_pointssd_final.pdf", width=10, height=5, units="in")
