% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genes.R
\name{process_genes}
\alias{process_genes}
\title{Create a mapping table for gene symbols and IDs}
\usage{
process_genes()
}
\value{
data.table
}
\description{
Maps the gene symbols with HGNC and Entrez IDs
}
\examples{
\dontrun{
genes <- process_data("genes")
validate(genes)
res <- update_table(genes)
}
}
\references{
\url{ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/json/hgnc_complete_set.json}
}
