# normalizing data for microbiome analysis
# and also approaches for compositional data analysis

##rarefy counts
##using vegan
rrarefy(otutable, sample_size, se=FALSE, MARGIN=1)
#or
rarefy(otutable, sample_size)
##or with phyloseq object
qiimez.rare = rarefy_even_depth(qiimez)

#use for comparing alpha diversity
#use for examining beta diversity using principle coordinates analysis (PCoA), with Bray-Curtis, unweighted UniFrac, Jaccard distances, or NMDS  (can also do these with center-log transformed data)
# can be used in DESeq, edgeR analysis too, but sensitive to 0s



## converting to compositional data from Gloor lab, based on https://github.com/ggloor/CoDaSeq/blob/master/simple_biplot.R)
## i.e. using log-transformation so that abundances are not constrained to sum to 1 (or 100%)


#first need bioconductor package manager

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

#need ALDEx2 as dependency for CoDaSeq
BiocManager::install("ALDEx2")

#then need devtools to install CoDaSeq
install.packages('devtools')
devtools::install_github('ggloor/CoDaSeq/CoDaSeq')
#doesn't work with older versions of dependencies

#install from source, first need to download CoDaSeq_0.99.1.tar.gz from github (https://github.com/ggloor/CoDaSeq)
install.packages("~/source/CoDaSeq-master/CoDaSeq_0.99.1.tar.gz", type="source", repos=NULL)


library(CoDaSeq)

# replace 0 values with the count zero multiplicative method and output counts
# if your otutable has samples in columns and OTUs in rows, this function expects the 
#samples to be in rows and OTUs to be in columns
# so the dataset is turned sideways on input (using t = transpose), and then back again on output
# you need to know which orientation your data needs to be in for each tool

otutable.counts <- t(cmultRepl(t(otutable), method="CZM", output="counts"))

# convert to proportions by sample (columns) using the apply function
# using this to fiter data by relative abundance
otutable.prop <- apply(otutable.counts, 2, function(x){x/sum(x)})

# Optionally, make a dataset where only taxa that are more abundant than 0.1% in all samples are included
# remove all taxa that are less than 0.1\% abundant in any sample
otutable.counts.point1perc <- otutable.counts[apply(otutable.prop, 1, min) > 0.001,]

# use counts data to make compositional dataset (center-log ratio transformation = log (x/ geometric mean), geometric mean = (exp (mean log (x)))  (FYI, aldex.clr is another function that generates center log-ratio transformed data)
# by samples (columns), but then transpose so samples in rows, compositional counts of OTUS in columns
# I think transpose data so that it is easier to add metadata info for grouping samples when they are the rows, or apply works row by row

otutable.coda.point1perc <- t(apply(otutable.counts.point1perc, 2, function(x){log(x) - mean(log(x))}))

#using log-transformed data in principal components analysis (PCA), uses euclidean distance
otutable.pca.point1perc <- prcomp(otutable.coda.point1perc)

#then plot biplot, Gloor lab CoDaSeq package has function for colouredBiplot, but can definitely customize similarly without this function.

pdf("simple_biplot.pdf", height=5, width=5)
coloredBiplot(pcx.abund, col="black", cex=c(0.6, 1), xlabs.col=conds$cond,
    arrow.len=0.05,
    xlab=paste("PC1 ", round (sum(pcx.abund$sdev[1]^2)/mvar(d.clr.abund),3), sep=""),
    ylab=paste("PC2 ", round (sum(pcx.abund$sdev[2]^2)/mvar(d.clr.abund),3), sep=""),
    expand=0.8,var.axes=T, scale=1, main="Biplot")
dev.off()





## can also normalize counts, corrects differences in sampling depth (so similar to rarefaction)

#install metagenomeSeq through bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("metagenomeSeq")


# this is based on the cumulative sum scaling of counts for an OTU/ASV
library(metagenomeSeq)
library(biomformat)

#otutable in .biom format
b <- read_biom("otutable.biom") 
#convert to format for metagenomeSeq functions
b <- biom2MRexperiment(b)
#Cumulative Sum Scaling Percentile Selection - calculates the percentile for which to sum counts up to and scale by
p = cumNormStatFast(b) # p = 0.5
#Cumulative Sum Scaling Normalization - calculates each column's quantile and calculates the sum up to and including that quantile
b = cumNorm(b, p = p)

#normalize raw count data
otutable.CSS.counts = MRcounts(b, norm = TRUE, log = TRUE)
exportMat(mat, file = file.path("otutable.CSS.tsv"))
exportStats(b, file = file.path("otutable.CSSstats.tsv"))

#metagenomeseq also has function for looking at time series, might be worth looking at
fitTimeSeries()
#see vignette('fitTimeSeries')


library(plyr)
# another function to calculate geometric mean, here 0 values are ignored
# used in DESeq2, as a work around to estimateSizeFactors
# probably shouldn't be used otherwise

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

#for use in DESeq2, calculate geometric mean over the normalized counts for each ASV (which are the rows of the data, want to determine which samples have ASV counts that are significantly different)
geoMeans = apply(counts(dataframe), 1, gm_mean)

#see DESeq2_notes.txt - another work around could be to add pseudo-count to zeros, then calculate gm_mean

