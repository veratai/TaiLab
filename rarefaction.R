#plot rarefaction curves

library(phyloseq)
library(ggplot2)
library(plyr)
library(vegan)
library(RColorBrewer)

#function to parse taxonomy, not greengenes format, split by ;
parse_taxonomy_vera <- function (char.vec){
  parse_taxonomy_default(strsplit(char.vec, ";", TRUE)[[1]])
}

#this function needs vegan package
## Rarefaction curve, ggplot style
ggrare <- function(physeq, step = 10, label = NULL, color = NULL, plot = TRUE, parallel = FALSE, se = TRUE) {
  ## Args:
  ## - physeq: phyloseq class object, from which abundance data are extracted
  ## - step: Step size for sample size in rarefaction curves
  ## - label: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  ## - color: (Optional). Default ‘NULL’. Character string. The name of the
  ##          variable to map to colors in the plot. This can be a sample
  ##          variable (among the set returned by
  ##          ‘sample_variables(physeq)’ ) or taxonomic rank (among the set
  ##          returned by ‘rank_names(physeq)’).
  ##
  ##          Finally, The color scheme is chosen automatically by
  ##          ‘link{ggplot}’, but it can be modified afterward with an
  ##          additional layer using ‘scale_color_manual’.
  ## - color: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  ## - plot:  Logical, should the graphic be plotted.
  ## - parallel: should rarefaction be parallelized (using parallel framework)
  ## - se:    Default TRUE. Logical. Should standard errors be computed. 
  ## require vegan
  x <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) { x <- t(x) }
  
  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  
  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)
  
  ## Get sample data 
  if (!is.null(sample_data(physeq, FALSE))) {
    sdf <- as(sample_data(physeq), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  
  ## Add, any custom-supplied plot-mapped variables
  if( length(color) > 1 ){
    data$color <- color
    names(data)[names(data)=="color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  if( length(label) > 1 ){
    labels$label <- label
    names(labels)[names(labels)=="label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  p <- ggplot(data = data, aes_string(x = "Size", y = ".S", group = "Sample", color = color))
  p <- p + labs(x = "Sample Size", y = "Species Richness")
  if (!is.null(label)) {
    p <- p + geom_text(data = labels, aes_string(x = "x", y = "y", label = label, color = color),
                       size = 4, hjust = 0)
  }
  p <- p + geom_line()
  if (se) { ## add standard error if available
    p <- p + geom_ribbon(aes_string(ymin = ".S - .se", ymax = ".S + .se", color = NULL, fill = color), alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}


setwd("/Users/vera/Projects/BeachMicrobiome/HakaiBeachJune2014/GoogleDrive")

#for 16S v4 EMP data, bacteria
#needed to change "taxonomy" column name in otu table to "Consensus Lineage"
#parse taxonomy as default, otherwise will look for greengenes format

otufile = "./16Sv4/16SV4_all/16SV4_all_otutable_4phyloseq.txt"
mapfile = "./16SV4/201406_16SV4_metadata_map_rmlow.txt"
trefile = "./16SV4/16SV4_all/16SV4_all_tree.tre"
qiimex = import_qiime(otufile, mapfile, trefile, parseFunction=parse_taxonomy_vera, showProgress = FALSE)
#not sure why, but tree doesn't import correctly

#import separately
tree <- read_tree("./16SV4/16SV4_all/16SV4_all_tree.tre")

#then create new phyloseq object with tree added, along with other stuff
qiimez <- phyloseq(otu_table(qiimex), sample_data(qiimex), tax_table(qiimex), tree)

#change factors for transect_site
#***for transect_site:  West Beach 1.5 mid_sw will change to mid_gw
#***for transect_site:  3rd 3 - gw seepage and 7th 2.5 gw seepage change to mid_sw
#and low tide flat should = low
#revalue
sample_data(qiimez)$transect_site_desc <- revalue(sample_data(qiimez)$transect_site_desc, c("gw_seepage" = "mid seawater interstitial", "high"="high interstitial", "high_so"="high sand", "low"="low interstitial", "low_so" = "low sand", "low_flat" = "low interstitial", "mid" = "mid interstitial", "mid_so" = "mid sand", "mid_bar" = "mid bar interstitial", "mid_sw" = "mid groundwater interstitial", "mid_sw_so" = "mid groundwater sand", "seawater" = "seawater", "swash" = "swash interstitial", "swash_so" = "swash sand") )
sample_data(qiimez)$transect_site_desc <- factor(sample_data(qiimez)$transect_site_desc, levels = c("high interstitial", "high sand", "mid interstitial", "mid sand", "mid bar interstitial", "mid groundwater interstitial", "mid groundwater sand", "mid seawater interstitial", "low interstitial", "low sand", "swash interstitial", "swash sand", "seawater") )


brewer.pal(n = 12, name = "Paired")
#need more pale colours for interstitial samples
brewer.pal(n = 8, name = "Pastel2")
#"#B3E2CD" "#FDCDAC" "#CBD5E8" "#F4CAE4" "#E6F5C9" "#FFF2AE" "#F1E2CC" "#CCCCCC"

#pick custom palette so interstitials are pale, sand is bolder
palette_brewer <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#F4CAE4", "#E31A1C", "#FDBF6F", "#FDCDAC", "#FF7F00", "#CAB2D6", "#6A3D9A", "#B15928")

#plot rarefaction curves
#p <- ggrare(qiimez, step = 1000, color = "beach", label = "anonymized_name", se = FALSE)
p <- ggrare(qiimez, step = 1000, color = "transect_site_desc", se = FALSE)
p <- p + scale_colour_manual(values=palette_brewer) + facet_wrap(~beach)
plot(p)
ggsave("16SV4_rarefaction.pdf")


