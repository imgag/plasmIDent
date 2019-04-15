#!/usr/bin/env Rscript
######################
# Create CIRCOS data #
# Calculate GC cont. #
######################

require(seqinr)

# Calculates the GC content in a sliding window for a character vector
GCwindow <- function(seq, width, circular=T, bases=c("c", "g")) {
  n <- length(seq)
  gc <- double(n)
  w <- floor(width/2)
  
  for (i in 1:length(seq)){

    # Extract subsequence in the current window
    # Start at the start of sequence if circular = T
    if (i < w) {
      window <- c(seq[(n-(w-i)):n], seq[1:i+w]) 
    } else if ((i+w) > n) {
      window <- c(seq[(-w + i):n], seq[1:(w-(n-i))])
    } else {
      window <- seq[(i-w):(i+w)]
    }
    
    # Calculate base occurences
    t <- table(window)
    b1 <- t[names(t)==bases[1]]
    b2 <- t[names(t)==bases[2]]
    
    gc[i] <- sum(b1, b2)/width
  }
  return(gc)
}

# Write a vector with numerical values as Circos data output
# Width parameter aggregates the data (mean value) to reduce the number of datapoints
aggregateCircos <- function(v, width, chr = "C") {
  
  n <- floor(length(v) / width)
  out <- character(n)
  
  for (i in 1:n){
   end <- min(i * width, length(v))
   start <- max(0, end - width)
   val <- mean(v[start:end])
   out[i] <- paste(chr, format(start, scientific=F), format(end, scientific=F), val, sep="\t") 
  }
  return(out)
}

# Get command line parameters
args = commandArgs(trailingOnly=TRUE)

# Load Plasmid sequences
fasta <- read.fasta(args[1])
#fasta <- read.fasta("data/testAssembly.fasta")

for (i in 1:length(fasta)) {
  print(paste("Reading contig", getName(fasta[[i]])))
  
  # Calculate GC content for Window sizes 50 and 1000
  GC50 <- GCwindow(getSequence(fasta)[[1]], 50)
  GC1000 <- GCwindow(getSequence(fasta)[[1]], 1000)
  
  # Write values to file
  vals1000 <- aggregateCircos(GC1000, 1000, getName(fasta[[i]]))
  vals50 <- aggregateCircos(GC50, 50, getName(fasta[[i]]))
  
  lapply(vals50, write, "gc50.txt", append=TRUE)
  lapply(vals1000, write, "gc1000.txt", append=TRUE)
  
}




