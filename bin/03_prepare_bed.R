#!/usr/bin/env Rscript

####################
#  Prepare reads   #
#  for CIRCOS data #
#  tracks from rgi #
####################

require(data.table)

args <- commandArgs(trailingOnly=TRUE)
bed <- args[1]
window <- as.numeric(args[2])
outFile <- file(args[3])
logTransform <- as.logical(args[4])
contigName <- args[5]
len <- as.numeric(args[6])

if (is.na(logTransform)) logTransform <- FALSE

s  <- 10

dt <- fread(bed)

# Write dummy file if bed is empty
if (nrow(dt) < 1) { 
  writeLines(c("0", "0", "0", "0"), sep="\t", outFile)
  close(outFile)
 
  } else {
  
  if (length(colnames(dt)) == 6 ) setnames(dt, c("chr", "start", "stop", "name", "qual", "strand"))
  if (length(colnames(dt)) == 4 ) setnames(dt, c("chr", "start", "stop", "strand"))
  
  dt[,chr := as.character(chr)]
  
  if (!is.na(contigName)) dt <- dt[chr == contigName,,]
  
  dt[strand == '+', strand := 1,]
  dt[strand == '-', strand := 0,]
  
  dt[,start := start - window]
  dt[,stop := stop - window]
  dt[start < 0, start := 0]
  dt[stop <0, stop := 0]
  dt <- dt[!(start == 0 & stop == 0),,]
  
  if (!is.na(len)) {
    dt[start > len, start := len]
    dt[stop > len, stop := len]
    dt <- dt[!(start == len & stop == len),,]
  }
  
  if (logTransform) {
    dt[,strand := log10(strand + 1),]
  }
  
  dt$start <- format(dt$start, drop0trailing = TRUE, scientific = FALSE)
  dt$stop <- format(dt$stop, drop0trailing = TRUE, scientific = FALSE)
  
  # Write to file
  write.table(dt[,.(chr, trimws(as.character(start)), trimws(as.character(stop)), strand)], quote=F, sep="\t", file=outFile, row.names = F, col.names=F)
}
