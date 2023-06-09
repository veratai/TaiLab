#explore data, summarize data
#using phyloseq
#useful websites:
# http://evomics.org/wp-content/uploads/2016/01/phyloseq-Lab-01-Answers.html#load-packages-import-data
# https://microucph.github.io/amplicon_data_analysis/html/phyloseq_operations.html
# https://joey711.github.io/phyloseq/preprocess.html
# updated June 2023

library(phyloseq)
library(qiime2R)
library(ggplot2)

#set the directory where files are located, or where files are to be written for this project
#e.g.
setwd("/Users/veratai/Projects/SunfishLakePlanktothrix/cyano16S_analysis/")

#import qiime2 files containing the associated metabarcode data as a phyloseq object

sunfish <- qza_to_phyloseq(
  features="sunfish_cyano16S_otutable.qza",
  tree="sunfish_cyano16S_tree.qza",
  taxonomy="sunfish_cyano16S_taxonomy.qza",
  metadata = "sunfish_metadata.txt"
)

#executing the phyloseq object (i.e. running sunfish) will give a summary of the data
#Did the data import as you expected?  Check number of samples, etc...
sunfish
#to view/use otu table
otu_table(sunfish)
#to view/use metadata
sample_data(sunfish)
#to view/use taxonomy
tax_table(sunfish)

#try running these commands, and see what info gets displayed
sample_variables(sunfish)
rank_names(sunfish) 

#***As you run more code, ALWAYS check the resulting object
#so that you understand what the code did
#so that you understand the information it contains
#are the results what you expected?

#get total number of reads per sample, as a data frame
reads_sample <- as.data.frame(sample_sums(sunfish))
#give column 1 a name
colnames(reads_sample)[1] <- 'total_reads'

#add metadata for a specific metadata column, e.g. location
reads_sample$location <- sample_data(sunfish)$location

#plot histogram showing frequency of samples with number of reads
reads_histo <- ggplot(reads_sample, aes(total_reads)) + geom_histogram() + ggtitle("Sequencing Depth")
#split plot based on metadata category
reads_histo + facet_wrap(~location)

#get summary stats for number of reads based on metadata category
tapply(reads_sample$total_reads, reads_sample$location, summary)


#get count of the number of otus per sample, i.e. otu/asv richness
#count number of rows, i.e. ASVs, where read count is > 0
asvs_sample <- apply(otu_table(sunfish),MARGIN = 2,function(x) sum(x > 0))
#can plot histogram as above, or get summary stats

#get count of the number of reads per otu/asv
reads_otus <- taxa_sums(sunfish)
#plot histogram, adjust breaks if needed
hist(reads_otus)
#can use to check or decide on filtering thresholds
#can also use this data to plot cumulative sums, 
#i.e. how many OTUs with increasing number of reads
#see:
# http://evomics.org/wp-content/uploads/2016/01/phyloseq-Lab-01-Answers.html#load-packages-import-data

#check prevalence, i.e. how many samples contain a particular ASV
#can be used as a filtering threshold, e.g. may want to remove any unique ASVs
#see code here:
# http://evomics.org/wp-content/uploads/2016/01/phyloseq-Lab-01-Answers.html#load-packages-import-data


#USEFUL FUNCTIONS 
#to organize, rename, check data

#subset data to a new phyloseq object, based data belonging to a metadata category (i.e. sample_data)
sunfish_north <- subset_samples(sunfish, location=="N")

#subset data to a new phyloseq object, based on data belonging to a particular taxonomic group (i.e. tax_table)
sunfish_cyanos <- subset_taxa(sunfish, Phylum=="Cyanobacteria")
#see the Classes remaining within the subsetted data
get_taxa_unique(sunfish_cyanos, "Class")

#find taxonomic data for a particular ASV
#convert phyloseq info to a dataframe
temp <- data.frame(tax_table(sunfish))
temp$asv <- row.names(temp)
#can also do it this way, with magrittr, tidyverse
temp <- tax_table(sunfish) %>%
  data.frame(asv = row.names(.))

#e.g. find taxonomic data for a particular asv
subset(temp, asv == "4ca22060a499cedff8b81b95d5978800")
#e.g. find taxonomic data for a particular genus
subset(temp, Genus == "Synechococcus")

#find ASV count data for a particular ASV
temp <- data.frame(otu_table(sunfish))
temp$asv <- row.names(temp)
subset(temp, asv == "4ca22060a499cedff8b81b95d5978800")


#to rename or order metadata
#revalue, i.e. rename levels
#for e.g.
sample_data(sunfish)$location <- revalue(sample_data(sunfish)$location, c("N" = "North", "S" = "South", "SSW" = "South", "E" = "East", "W" = "West") )
#set custom order
sample_data(sunfish)$location <- factor(sample_data(sunfish)$location, levels = c("North", "East", "South", "West") )


#CLEAN-UP TAXON NAMES

#check names for taxonomic ranks, may need to edit where taxa names are uninformative (e.g. unculturable)
#list rank_names used in taxonomy
rank_names(sunfish)

#For 16S data, from Silva
#7 ranks: "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"
#for 16S data, need to filter out "chloroplasts", "mitochondria"
#also want to filter out taxa that are "Eukaryota" and "Unassigned" at the highest rank, i.e. were not identified as either Bacteria or Archaea

#For 18S data, from PR2
#9 taxonomic levels (domain to species) (domain, division (this 2nd rank is new - but PR2 referring to it as subdivision?), kingdom?, phylum, class, order, family, genus, species)
#includes selection of Bac and Arc (as of v4.12.0), so would want to filter these out.
#plastids labeled as:  Eukaryota:plas
#mitochondria labeled as: Eukaryota:mito

#list unique taxa names in Kingdom, check if there are Unassigned or Eukaryote ASVs
get_taxa_unique(sunfish, "Kingdom")
#check how many Unassigned ASVs
sunfish_unassigned <- subset_taxa(sunfish, Kingdom=="Unassigned")
#count number of unassigned ASVs
unassigned_asvs <- taxa_sums(sunfish_unassigned)
sum(unassigned_asvs)
#calculate percentage of reads that are unassigned
(sum(unassigned_asvs)/ sum(taxa_sums(sunfish)) ) * 100

#list unique taxa names in Order, check if there are Chloroplast ASVs
get_taxa_unique(sunfish, "Order")

#list unique taxa names in Family, check if there are Mitochondria ASVs
get_taxa_unique(sunfish, "Family")

#e.g. for 16S data, classified with Silva database, filtering out ASVs classified as "Eukaryota", Unassigned", "Chloroplast" and "Mitochondria"
sunfish_prokonly <- subset_taxa(sunfish, Kingdom!="d__Eukaryota" & Kingdom!="Unassigned" & Order!="Chloroplast" & Family!="Mitochondria")
#check, are Chloroplast ASVs gone? 
get_taxa_unique(sunfish_prokonly, "Order")

##do something similar if need to filter out plastid or mitochondrial ASVs from 18S data


# to replace/correct taxonomic names
# e.g. for Rank2, find entries with "uncultured_bacteria" and replace with "Bacteria_X"
tax_table(sunfish)[,2][tax_table(sunfish,)[,2]=="uncultured_bacteria"] <-"Bacteria_X"
# see code below for more comprehensive re-naming for uncultured, unidentified ASVs

#e.g. replacing names with a conditional, e.g. depending on the name in another rank
#for Rank 2 replace Proteobacteria with Rank 3 names for Alphaproteobacteria, Deltaproteobacteria, or Gammaproteobacteria as appropriate
#Rank2 is in position 2 of the tax_table, if Rank3 is Alphaproteobacteria, replace Rank2 with Alphaproteobacteria
tax_table(sunfish)[,2][tax_table(sunfish,)[,3]=="Alphaproteobacteria"] <-"Alphaproteobacteria"
tax_table(sunfish)[,2][tax_table(sunfish,)[,3]=="Deltaproteobacteria"] <-"Deltaproteobacteria"
tax_table(sunfish)[,2][tax_table(sunfish,)[,3]=="Gammaproteobacteria"] <-"Gammaproteobacteria"

#may want to remove prefixes from 16S classification
#e.g. haven't tested this code yet...
tax_table(sunfish)[,1] <- gsub("d__", "", tax_table(sunfish)[,1])
#these ones may not be applicable anymore, prefixes already removed
tax_table(sunfish)[,2] <- gsub("p__", "", tax_table(sunfish)[,2])
tax_table(sunfish)[,3] <- gsub("c__", "", tax_table(sunfish)[,3])
tax_table(sunfish)[,4] <- gsub("o__", "", tax_table(sunfish)[,4])
tax_table(sunfish)[,5] <- gsub("f__", "", tax_table(sunfish)[,5])
tax_table(sunfish)[,6] <- gsub("g__", "", tax_table(sunfish)[,6])
tax_table(sunfish)[,7] <- gsub("s__", "", tax_table(sunfish)[,7])


#FIXING TAXON NAMES that are sometimes not actual taxon names, are not useful or informative
#e.g. uncultured_bacteria, unidentified_marine, gut_metagenome
# need to replace these, using name from taxonomic rank above it
#modified from code found here:
#https://github.com/joey711/phyloseq/issues/850

library("tidyverse")

#define function to replace "NA" taxa names AND based on string search (set by replace), with taxa name from the rank above
#if replace is not defined when calling function, will still replace NA taxa names by default with name from rank above
#copy and paste this function as is

name_bad_taxa <- function(ps_obj, include_rank = T, bad_label = "Unidentified_<tax> (<rank>)", replace = "NA"){
  # Check arguments
  if(!grepl("<tax>", bad_label)){
    stop("Error: include '<tax>' in the na_label")
  }
  if (include_rank){
    if(!grepl("<rank>", bad_label)){
      stop("Error: include_rank = TRUE; include '<rank>' in the na_label")
    }
  } else {
    if(grepl("<rank>", bad_label)){
      stop("Error: include_rank = FALSE; remove '<rank>' from the na_label")
    }
  }
  
  # Convert to long data, columns are row_name(i.e. ASVs), rank and taxa names
  # Ranks must be ordered according to taxonomic hierarchy (e.g. Rank1, Rank2, Rank3 works)
  taxa_long <- tax_table(ps_obj) %>%
    data.frame(asv = row.names(.)) %>%
    pivot_longer(!asv,
                 names_to = "rank",
                 values_to = "tax")
  
  # For items containing replace change to NA
  taxa_long$tax[grepl(replace, taxa_long$tax)] <- NA
  
  # Fill in NAs using the value above
  taxa_long <- taxa_long %>%
    mutate(na = is.na(tax)) %>%
    group_by(row_name) %>%
    fill(tax)
  
  # Create na_labels, replacing <tax> with tax (value above)
  taxa_long <- taxa_long %>%
    mutate(expr = ifelse(na,
                         bad_label,
                         tax), 
           bad_label = str_replace(expr, "<tax>", tax))
  
  # Add the last annotated rank
  if (include_rank){
    taxa_long <- taxa_long %>%
      mutate(last_rank = ifelse(na,
                                NA,
                                rank)) %>%
      fill(last_rank) %>%
      mutate(bad_label = str_replace(bad_label, "<rank>", last_rank))
  }
  # Convert back to tax_table
  taxa_mat <- taxa_long %>%
    select(asv, rank, bad_label) %>%
    pivot_wider(names_from = rank, values_from = bad_label) %>%
    as.matrix()
  row.names(taxa_mat) <- taxa_mat[,"row_name"]
  taxa_mat <- taxa_mat[,colnames(taxa_mat) != "row_name"]
  tax_table(ps_obj) <- taxa_mat
  
  return(ps_obj)
}


#set fix_taxa to find the string in the taxa names that need fixing
fix_taxa = "uncultured"

#copy phyloseq object to new variable name
sunfish_curated <- sunfish_prokonly
#run name_na_taxa function
sunfish_curated <- name_bad_taxa(sunfish_prokonly, bad_label = "unidentified_<tax>", replace=fix_taxa)

#run fix_taxa and name_bad_taxa function again to replace a different taxon name
#e.g. repeat with
fix_taxa = "unidentified"
#and repeat with
fix_taxa = "metagenome"

#check taxa names in each rank to see if there are other names that need replacing
get_taxa_unique(sunfish_prokonly, "Species")



