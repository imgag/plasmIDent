#!/usr/bin/env Rscript

########################
# Create summyry table #
#											 #
########################

require(data.table)
require(reshape2)
require(stringr)

args <- commandArgs(trailingOnly=TRUE)

if (length(args) >= 1) {
	df.list <- lapply(args, fread)
	names(df.list) <- basename(args)
	data.table::rbindlist(df.list, idcol = "sample")
}