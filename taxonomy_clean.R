# another way to clean up taxonomy (i.e. add higher taxonomic ranks where lower ones could not be identified, clean up "uncultured", "unidentified" identifications)
#but try code in phyloseq_explore.R first

#based on here
#https://github.com/joey711/phyloseq/issues/850


#convert phyloseq object as a data frame
tax.clean <- data.frame(tax_table(sunfish))

#change everything to characters - is this necessary?
for (i in 1:7){ tax.clean[,i] <- as.character(tax[,i])}

#change NA identifications or __ to empty string
tax.clean[is.na(tax.clean)] <- ""
tax.clean[tax.clean=="__"] <- ""

#get better identification for ASVs ranking as "uncultured", "unidentified" taxon or "metagenome"
#add code for other names that need fixing are "archaeon", "gut_microbe"
#will need to check what needs fixing based on your specific data

#need to check this code...
#replaces taxon name with an empty string
tax.clean <- data.frame(lapply(tax.clean, function(x) {
  gsub("uncultured[_A-Za-z0-9]+", "", x)
}))
tax.clean <- data.frame(lapply(tax.clean, function(x) {
  gsub("unidentified[_A-Za-z0-9]+", "", x)
}))
tax.clean <- data.frame(lapply(tax.clean, function(x) {
  gsub("metagenome[_A-Za-z0-9]+", "", x)
}))

#for everything that is "", use the name for the previous taxon
#and will fill this downwards for the rest of the taxonomic levels
#I think...

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Unidentified_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Unidentified_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Unidentified_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Unidentified_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Unidentified_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Unidentified",tax.clean$Genus[i], sep = "_")
  }
}


#replace cleaned taxonomy back into phyloseq object
tax_table(sunfish) <- as.matrix(tax.clean)

