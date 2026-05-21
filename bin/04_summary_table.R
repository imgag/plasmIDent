#!/usr/bin/env Rscript

######################
# Create CIRCOS data #
#  tracks from rgi   #
######################

require(data.table)
require(seqinr)

# Read data
args <- commandArgs(trailingOnly=TRUE)
fasta <- args[1]
rgi <- args[2]
cov <- args[3]
gc <- args[4]
padding <- as.integer(args[5])


seq <- read.fasta(fasta)
dt_rgi <- fread(rgi)

dt_cov <- fread(cov)
setnames(dt_cov, c("contig", "start", "stop", "cov"))
dt_cov[,contig := as.character(contig)]

dt_gc <- fread(gc)
setnames(dt_gc, c("contig", "start", "stop", "gc"))
dt_gc[,contig := as.character(contig)]

# Get contig ids from GC (already filtered)
id <- dt_gc[,unique(contig),]
n <- length(id)

# Calculate coverage normalized by all contigs != plasmids
baseCov <- dt_cov[!(contig %in% id), median(cov),]
val_cov <- dt_cov[(contig %in% id), .(cov = median(cov)/baseCov), by = contig]

# Extract resistance genes
getID <- function(s) {
  #s <- "VBDFDF_3_1_234"
  sub("_[^_]+$", "", s)
}

dt_ar <- dt_rgi[, .(
  contig = vapply(Contig, getID, character(1)),
  ar_genes = vapply(Best_Hit_ARO, toString, character(1)),
  aro = vapply(ARO, toString, character(1))
)]
dt_ar[, contig := as.character(contig)]

# Put it all together
dt <- dt_gc[,.(gc = mean(gc)), by = contig]
dt[, length := vapply(seq[contig], getLength, integer(1)) - 2 * padding]
dt <- dt[val_cov, on = "contig"]

dt[,cov := round(cov, 2),]
dt[,gc := round(gc, 2),]

dt <- merge(dt, dt_ar, by = "contig", all.x = T)
fwrite(dt, sep="\t", file ="contig_summary.txt")
