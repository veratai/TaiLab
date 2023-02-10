#useful phyloseq functions

setwd("/Users/veratai/Projects/SunfishLakePlanktothrix/cyano16S_analysis/")

sunfish <- qza_to_phyloseq(
  features="sunfish_cyano16S_filteredtable_ASVmin7.qza",
  tree="unrooted-tree-sunfish-cyano16S.qza",
  taxonomy="sunfish_cyano16S_taxonomy.qza",
  metadata = "sunfish_metadata_VT.txt"
)


#list rank_names used in taxonomy
rank_names(sunfish)
#list unique taxa names in Rank2
get_taxa_unique(sunfish, "Phylum")

#just get cyano data
sunfish_cyanos <- subset_taxa(sunfish, Phylum=="Cyanobacteria")
get_taxa_unique(sunfish_cyanos, "Class")
plot_bar(sunfish_cyanos, fill="Class")

#just get oxyphoto data
sunfish_oxyphoto <- subset_taxa(sunfish, Class=="Oxyphotobacteria")
get_taxa_unique(sunfish_oxyphoto, "Order")
plot_bar(sunfish_oxyphoto, fill="Family")

##change factors for a metadata category, e.g. transect_site
#***for transect_site:  West Beach 1.5 mid_sw will change to mid_gw
#***for transect_site:  3rd 3 - gw seepage and 7th 2.5 gw seepage change to mid_sw
#and low tide flat should = low
#revalue
sample_data(qiimez)$transect_site_desc <- revalue(sample_data(qiimez)$transect_site_desc, c("gw_seepage" = "mid seawater interstitial", "high"="high interstitial", "high_so"="high sand", "low"="low interstitial", "low_so" = "low sand", "low_flat" = "low interstitial", "mid" = "mid interstitial", "mid_so" = "mid sand", "mid_bar" = "mid bar interstitial", "mid_sw" = "mid groundwater interstitial", "mid_sw_so" = "mid groundwater sand", "seawater" = "seawater", "swash" = "swash interstitial", "swash_so" = "swash sand") )
sample_data(qiimez)$transect_site_desc <- factor(sample_data(qiimez)$transect_site_desc, levels = c("high interstitial", "high sand", "mid interstitial", "mid sand", "mid bar interstitial", "mid groundwater interstitial", "mid groundwater sand", "mid seawater interstitial", "low interstitial", "low sand", "swash interstitial", "swash sand", "seawater") )

