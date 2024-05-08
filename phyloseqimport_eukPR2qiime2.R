#importing qiime2 files as a phyloseq object for eukaryotic data classified with PR2
#or, for any other classification where there are more than the standard 7 ranks

library(phyloseq)
library(qiime2R)

setwd("/Users/veratai/Projects/HIV_Uganda")

#this will only import the first 7 taxonomic ranks, but PR2 has 10 ranks
#hiveuk <- qza_to_phyloseq(
#  features="hivall_eukF_otutable.qza",
#  taxonomy="hivall_eukF_taxonomy.qza",
#  metadata = "hivall_eukmetadata.txt",
#  tree = "tree.qza"
#)

#could import otu table and metadata (and tree) as a phyloseq object first
#tempphy <- qza_to_phyloseq(
#  features="hivall_eukF_otutable.qza",
#  metadata = "hivall_eukmetadata.txt",
#  tree = "tree.qza"
)
#then merge with parsed tax table (see below)
#merge_phyloseq(otu_table(tempphy), sample_data(tempphy), phy_tree(temphy), tax_table(taxT) )

#but going to do this by importing each qza file as separate objects first
features <- read_qza("hivall_eukF_otutable.qza")

metadata <- read.table("hivall_eukmetadata.txt", header=TRUE)
#for metadata, rownames need to be the sample names
rownames(metadata) <- metadata$sample.id

tax = read_qza("hivall_eukF_taxonomy.qza")

#including code for tree, but haven't tested this yet
tree = read_qza("tree.qza")


## to get all 10 taxonomic ranks imported
#tax is a list of lists - important info is in $data, 
#view first bit of data, i.e. ASV names
head(tax$data)

#view ASV IDs
head(tax$data$Feature.ID)
#view taxonomy strings
head(tax$data$Taxon)

#split the taxonomy string into vector of taxonomic names
tax_split <- strsplit(tax$data$Taxon, ";", TRUE)

#get ASV names as a vector
asvs <- tax$data$Feature.ID
#create vector for names of the 10 taxonomic ranks in PR2
tax_ranks <- c("Domain", "Supergroup", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Subspecies")
#create empty matrix
taxT = matrix(, nrow = length(tax_split), ncol = 10, dimnames=list(asvs, tax_ranks))
#loop through each taxonomy vector in the list, and add as a row to the matrix
#set number of ranks as 10, based on ranks in PR2
for (i in 1:length(tax_split)) {
  #set each vector to have 10 items, i.e. 10 taxonomic ranks
  #fill empty ranks with NA, adds NA so that vector is length 10
  length(tax_split[[i]]) <- 10
  #go through each taxon name
  for (j in 1:length(tax_split[[i]])) {
    #check if taxon name is NA
    if (is.na(tax_split[[i]][j])) {
      #replace with taxon name from previous rank
      tax_split[[i]][j] <- tax_split[[i]][j-1]
    }
  }
  taxT[i, ] <- tax_split[[i]]
}

#check that format is correct, can be read as a phyloseq tax_table
head(tax_table(taxT))

#make phyloseq object
hiveuk <- merge_phyloseq(otu_table(features$data, taxa_are_rows=TRUE), sample_data(metadata), phy_tree(tree$data), tax_table(taxT) )

#get some info
rank_names(hiveuk)
get_taxa_unique(hiveuk, "Domain")


