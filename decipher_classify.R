#!/usr/bin/env RScript

#classify sequences using a DECIPHER trained classifier
#https://rdrr.io/bioc/DECIPHER/f/inst/doc/ClassifySequences.pdf

#path to sequences to classify
fas <- "<<path to FASTA file>>"
# read the sequences into memory
test <- readDNAStringSet(fas)
#test sequences cannot contain gap (“-” or “.”) characters
#remove with the RemoveGaps function:
test <- RemoveGaps(test)

#load training set
trainingSet <- readRDS("<<path to trained_classifier>>")

#Classify
#The most important (optional) arguments are the type of output, the strand used in testing, 
#the confidence threshold of assignments, and the number of processors to use. 
#Here, we are going to request the "extended" (default) output type that allows for plotting the results
#but there is also a "collapsed" type that might be easier to export (see section 4.4 below). 
#Also, we know that all of the test sequences are in the same (“+” strand) orientation as the training sequences
#so we can specify to only look at the "top" strand rather than the default of "both" strands (i.e., both “+” and “-” strands). 
#This makes the classification process over twice as fast. 
#We could also set processors to NULL to use all available processors.

#The threshold of 60% is recommended at the default confidence threshold. 
#Confidence levels are informally defined as 70% (stringent), 60% (cautious), 50% (sensible), and 40% (lenient). 
#Using a threshold of 0% will report classifications down to all rank levels. 
#Note that the test sequences should generally be fully overlapped by the information in the training sequences. 
#In this way, the training sequences can be longer than the test sequences, but the reverse situation would result in lower confidences.

ids <- IdTaxa(test,
              trainingSet,
              type="extended",
              strand="top",
              threshold=60,
              processors=1)


#returned object ids is a list
ids
#may want to convert to a dataframe, or other format depending on downstream use



#see notes in #https://rdrr.io/bioc/DECIPHER/f/inst/doc/ClassifySequences.pdf
#to get examples of viewing the results

#plot a piechart of relative abundance of the taxonomic groups assigned to test sequences. 
#also displays the training taxonomic tree, with edges colored where they match the taxonomic groups shown in the pie chart. 
#Note, can omit the taxonomic tree by omitting the trainingSet in the function
plot(ids, trainingSet)

#also check decipher notes to see code to get taxonomy bar plot



#OPTIONAL - exporting data, via pr2 code

#have set taxonomic levels, as in UNITE database
taxo_levels <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

# Transform to a dataframe
n_seq <- length(ids)
df_rows <- list()

# Go through all the elements of the list

for(i in 1:n_seq){
  seq_name <- names(ids[i])
  taxonomy<- ids[[i]]$taxon
  confidence <- ids[[i]]$confidence
  df_rows[[i]] = data.frame(seq_name, taxonomy, confidence, taxo_level=c("Root", taxo_levels))
}

df <- reduce(df_rows, bind_rows) %>% 
  filter(taxo_level !="Root") %>% 
  pivot_wider(names_from = taxo_level, values_from = c(taxonomy, confidence))

# Save to file
saveRDS(df, "name_of_file_assigned")

#or export as table

