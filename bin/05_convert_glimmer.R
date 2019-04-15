#!/usr/bin/env Rscript
#######################
# Conver glimmer out  #
# to circos .bed file #
#######################

args <- commandArgs(trailingOnly=TRUE)
file <- args[1]

lines <- readLines(file)
contigID <- ""
splitLine <- ""
strand <- ""

require(stringr)


outFile <- file("genes.txt", "w")

for (line in lines) {
  if (startsWith(line, ">")) {
    contigID <- unlist(str_split(line, ">"))[2]
    #contigID <- str_replace(contigID, '_', '-')
    contigID <- str_split(contigID, ' ')[[1]][1]
  } else {
    line <- str_squish(line)
    splitLine <- unlist(str_split(line, " "))
    
    start <- as.integer(splitLine[2])
    stop <- as.integer(splitLine[3])
    
    if (start<stop) {
      # Forward: 1
      strand = 1
    } else {
      # Backwards: 0
      tmp <- stop
      stop <- start
      start <- tmp
      strand = 0
    }
    
    # Filter out circle closing genes since they cannot be visualized in contig
    if ((stop - start) < 5000 ) {
      print(paste(contigID, start, stop, strand))
      write(paste(contigID, start, stop, strand, sep="\t"), file=outFile, append = TRUE)
    }
  }
}

close(outFile)  
  
  
