library(tidyverse)
library(qiime2R)
library(ggplot2)
library(ggnewscale)


setwd("/Users/veratai/Documents/TaiLab/SOPs/DEICODE")

#read in ordination generation by DEICODE, implemented in qiime2
deicodeSunfish <- read_qza(file="sunfish_ordination.qza")

#read in taxonomy table from qiime2
tax<-read_qza("~/Projects/SunfishLakePlanktothrix/cyano16S_analysis/sunfish_cyano16S_cyanos_phytoref_oxyphoto99_taxonomy.qza")$data %>%
  rename(FeatureID=Feature.ID)

#give new names to Taxon, names are too long for plotting
#here could use the microbiomeutilities package and apply this to the phyloseq object to get the best taxon hit
# or an split the string to just get a particular rank of taxonomy

#Here just making up random names, because the taxonomy in these data is all similar, and wanted to demo using Taxon to specify colour
tax$Taxon <- rep(c("aaa", "bbb", "ccc", "ddd", "eee", "fff", "ggg", "hhh", "iii", "jjj"),times=c(15,15,15,15,15,15,15,15,15,13))


#read in metadata table
meta<-read_tsv(file = "~/Projects/SunfishLakePlanktothrix/cyano16S_analysis/sunfish_metadata_VT.txt") 
#rename(SAMPLEID=`sampleid`) %>%
#filter(sampleid!="#q2:types")

#make sure sample ID column name is exactly the same between ordination and metadata
#check colnames of deicodeSunfish$data$Vectors
colnames(deicodeSunfish$data$Vectors)
#check colnames of colnames(meta)
colnames(meta)
#set first column name of meta to be the same as first column name of ordination
#i.e. the SampleID column
colnames(meta)[1]  <- colnames(deicodeSunfish$data$Vectors)[1] 




#create the base plot with only the arrows
#in this code, deicodeSunfish$data$Species is joined with tax
#so that can plot FeatureID alongside specific columns in tax, i.e. Taxon 

baseplot<-
  ggplot() +
  theme_bw() +
  xlab(paste(round(100*deicodeSunfish$data$ProportionExplained[1],2),"%")) +
  ylab(paste(round(100*deicodeSunfish$data$ProportionExplained[2],2),"%")) +
  ggtitle("DEICODE biplot of Sunfish data") +
  geom_segment(data=deicodeSunfish$data$Species %>% 
                 mutate(a=sqrt(PC1^2+PC2^2)) %>%      # calculate the distance from the origin
                 top_n(10, a) %>%     # keep 8 furthest away points
                 mutate(PC1=PC1*0.3, PC2=PC2*0.3) %>%      # scale arrows linearly... is this ok? 
                 left_join(tax),      # joining tax table to ordination$Species table
               aes(x=0, xend=PC1, y=0, yend=PC2, color=Taxon),
               arrow = arrow(length = unit(0.3,"cm"))
  )

#overlay samples
#PC data is joined with metadata
baseplot +
  ggnewscale::new_scale_colour() +
  geom_point(
    data=deicodeSunfish$data$Vectors %>% left_join(meta),
    aes(x=PC1, y=PC2, colour=SampleID), 
    shape=19, size=4.0
  )

#save ggplot as pdf
ggsave("sunfish_deicode_plot_top10.pdf", useDingbats=FALSE)


