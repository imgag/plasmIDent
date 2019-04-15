#!/usr/bin/env Rscript
######################
# Create CIRCOS data #
# Calculate GC cont. #
######################

require(seqinr)
require(data.table)

# Calculates the GC content in a sliding window for a character vector
GCwindow <- function(seq, width, circular=T, bases=c("c", "g")) {
  n <- floor(length(seq)/width)
  gc <- double(n)
  gcskew <- double(n)
  w <- floor(width/2)
  j <- 1
  sumgc <- 0
  diffgc <- 0
  for (i in 1:length(seq)){
    
    if (seq[i] %in% c("c", "g")) {
      sumgc <- sumgc + 1
      diffgc <- diffgc + 1
    }
    if (seq[i] == "c") {
      diffgc <- diffgc - 2
    }
    if ((i %% width) == 0 ) {
      gc[j] <- sumgc/width
      gcskew[j] <- diffgc/sumgc
      sumgc <- 0
      diffgc <- 0
      j <- j + 1
    }  
  }
  return(list(gc, gcskew))
}

# Write a vector with numerical values as Circos data output
# Width parameter aggregates the data (mean value) to reduce the number of datapoints
# Also removes sequence padding
aggregateCircos <- function(v, width, chr = "C", len) {
  
  n <- ceiling(len/width)
  out <- character(n+2)
  
  gapVal <- mean(v[1], v[length(v)])
  out[1] <- paste(chr, 0, 1, gapVal, sep="\t")
  j <- 2
  
  for (i in 1:length(v)){
    end <- min(i * width)
    start <- max(0, (end - width))
    val <- v[i]
    
    if (end > padding) {
      if (start < (len + padding)) {
        prstart <- format(max(0, start-padding), scientific = F)
        prend <- format(min(len, end-padding), scientific = F)
        val <- format(val, scientific = F)
        out[j] <- paste(chr, prstart, prend, val, sep="\t") 
        j <- j+1
      }
    }
    
  }
  
  out[n+2] <- paste(chr, len-1, len, gapVal, sep="\t")
  
  return(out)
}



# Get command line parameters
args = commandArgs(trailingOnly=TRUE)

# Load Plasmid sequences
fasta <- read.fasta(args[1])
padding <- as.integer(args[2])

for (i in 1:length(fasta)) {
  print(paste("Reading contig", getName(fasta[[i]])))
  
  # Calculate GC content for Window sizes 50 and 1000
  GC50 <- GCwindow(getSequence(fasta)[[i]], 50)
  GC1000 <- GCwindow(getSequence(fasta)[[i]], 500)
  
  GCskewsum1000 <- lapply(GC1000[2], cumsum)
  GCskewsum50 <- lapply(GC50[2], cumsum)
  
  len <- length(getSequence(fasta[[i]])) - 2*padding
  
  print(paste("GC vectors calculated, len: ", len)) 
  # Write values to file
  gc1000 <- aggregateCircos(GC1000[[1]], 500, getName(fasta[[i]]), len)
  gcskew1000 <- aggregateCircos(GC1000[[2]], 500, getName(fasta[[i]]), len)
  gcskewsum1000 <- aggregateCircos(GCskewsum1000[[1]], 500, getName(fasta[[i]]), len)
  
  gc50 <- aggregateCircos(GC50[[1]], 50, getName(fasta[[i]]), len)
  gcskew50 <- aggregateCircos(GC50[[2]], 50, getName(fasta[[i]]), len)
  gcskewsum50 <- aggregateCircos(GCskewsum50[[1]], 50, getName(fasta[[i]]), len)
  
  # Shift values for padded fasta
  #print(vals1000)
  #print(vals50)
  print("Values aggregating, writing to file now ")

  lapply(gc50, write, "gc50.txt", append=TRUE)
  lapply(gc1000, write, "gc1000.txt", append=TRUE)
  lapply(gcskew50, write, "gcskew50.txt", append=TRUE)
  lapply(gcskew1000, write, "gcskew1000.txt", append=TRUE)
  lapply(gcskewsum50, write, "gcskewsum50.txt", append=TRUE)
  lapply(gcskewsum1000, write, "gcskewsum1000.txt", append=TRUE)

}

testplot <- function (val50, val1000, w1 = 50, w2 = 500, pad=padding) {
  len()
  plot((1:length(val50))*w1, val50, type="l", col="black")
  plot((1:length(val1000))*w2, val1000, type="l", add=TRUE, col="red")

}
    



