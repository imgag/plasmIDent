#!/usr/bin/env Rscript

########################
# Create summary table #
########################

require(data.table)

escape_html <- function(x) {
	x <- as.character(x)
	x[is.na(x)] <- ""
	x <- gsub("&", "&amp;", x, fixed = TRUE)
	x <- gsub("<", "&lt;", x, fixed = TRUE)
	x <- gsub(">", "&gt;", x, fixed = TRUE)
	gsub('"', "&quot;", x, fixed = TRUE)
}

write_html_table <- function(dt, html_path) {
	cols <- names(dt)
	header <- paste(sprintf("<th>%s</th>", escape_html(cols)), collapse = "")

	rows <- apply(dt, 1, function(row) {
		cells <- vapply(seq_along(row), function(idx) {
			value <- row[[idx]]
			if (cols[[idx]] == "id") {
				sprintf("<td>%s</td>", ifelse(is.na(value), "", value))
			} else {
				sprintf("<td>%s</td>", escape_html(value))
			}
		}, character(1))
		sprintf("<tr>%s</tr>", paste(cells, collapse = ""))
	})

	html <- c(
		"<!DOCTYPE html>",
		"<html lang='en'>",
		"<head>",
		"<meta charset='utf-8'>",
		"<title>Summary table of all replicons with resistance genes</title>",
		"<style>",
		"body { font-family: sans-serif; margin: 2rem; }",
		"table { border-collapse: collapse; width: 100%; }",
		"caption { font-size: 1.2rem; font-weight: 600; margin-bottom: 1rem; }",
		"th, td { border: 1px solid #d0d0d0; padding: 0.5rem; text-align: left; }",
		"th { background: #f3f3f3; }",
		"tr:nth-child(even) { background: #fafafa; }",
		"</style>",
		"</head>",
		"<body>",
		"<table>",
		"<caption>Summary table of all replicons with resistance genes</caption>",
		sprintf("<thead><tr>%s</tr></thead>", header),
		sprintf("<tbody>%s</tbody>", paste(rows, collapse = "")),
		"</table>",
		"</body>",
		"</html>"
	)

	writeLines(html, con = html_path)
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
	stop("Usage: 06_summary_samples.R <aro_index.tsv> <sample_summary.csv> [...]")
}

aro_index_path <- args[1]
summary_paths <- args[-1]

df_list <- lapply(summary_paths, fread)
names(df_list) <- sub("_summary\\.csv$", "", basename(summary_paths))

df <- rbindlist(df_list, idcol = "sample_name", fill = TRUE)
fwrite(df, sep = "\t", file = "rgi_summary_samples.tsv")

aro <- fread(aro_index_path)
aro_dt <- data.table(
	aro = sub("^ARO:", "", aro[["ARO Accession"]]),
	id = as.character(aro[["CVTERM ID"]])
)

df[, aro := as.character(aro)]
dfm <- merge(df[!is.na(aro) & nzchar(aro)], aro_dt, by = "aro", all.x = TRUE)
dfm[, id := ifelse(
	is.na(id) | !nzchar(id),
	"",
	sprintf("<a href='https://card.mcmaster.ca/ontology/%s' target='_blank'>%s</a>", id, id)
)]

write_html_table(dfm[, !"aro"], "rgi_summary_samples.html")