#!/usr/bin/env Rscript

########################
# Create summary table #

########################

require(data.table)
require(stringr)
require(knitr)

args <- commandArgs(trailingOnly=TRUE)

if (length(args) >= 1) {
	df.list <- lapply(args, fread)
	names(df.list) <- str_remove(basename(args), "_summary.csv")
	
	# tsv table
	df <- data.table::rbindlist(df.list, idcol = "sample_name")
	fwrite(df, sep = "\t", file = "rgi_summary_samples.tsv")
	
	# html table, use aro_index.tsv from the card db to match aro accession to ID
	#<a href="https://www.thesitewizard.com/" target="_blank">thesitewizard.com</a>
	# https://card.mcmaster.ca/ontology/43314
	aro <- fread("/aro_index.csv")
	aro_dt <- data.table(aro = str_remove(aro$`ARO Accession`, "^ARO:"), 
											id = aro$`CVTERM ID`)
	
	df$aro <- as.character(df$aro)
	dfm <- merge(df, aro_dt, by = "aro")
	dfm$id <- paste0("<a href='https://card.mcmaster.ca/ontology/", dfm$id, "'", " target='_blank'>", dfm$id, "</a>")
	
	cat(knitr::kable(dfm[, !c("aro")],
									 format = "html", 
									 escape = FALSE, 
									 align = "r", 
									 caption = "Summary table of all replicons with resistance genes"), 
			file = "rgi_summary_samples.html")
	
} else {
	stop("No arguments provided")
}