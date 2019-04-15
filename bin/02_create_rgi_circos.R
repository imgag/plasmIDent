#!/usr/bin/env Rscript

######################
# Create CIRCOS data #
#  tracks from rgi   #
######################

require(data.table)

args <- commandArgs(trailingOnly=TRUE)
file <- args[1]
rgi <- fread(file)

#rgi[, Contig := tstrsplit(Contig, "_")[[1]],]
#lapply(list(test), function(x) gsub('_[0-9]*$', '', x))
rgi[, Contig := lapply(Contig,  function(x) gsub('_[0-9]*$', '', x))]

# Export RGI
fwrite(rgi[,.(Contig, Start, Stop, Best_Hit_ARO)], sep="\t", file ="rgi.txt", col.names = F)
fwrite(rgi[,.(Contig, Start, Stop)], sep="\t", file ="rgi_span.txt", col.names = F)
