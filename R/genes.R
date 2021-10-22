#' Create a mapping table for gene symbols and IDs
#'
#' Maps the gene symbols with HGNC and Entrez IDs
#'
#' @references
#' \url{ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/json/hgnc_complete_set.json}
#'
#' @return data.table
#' @examples
#' \dontrun{
#' genes <- process_data("genes")
#' validate(genes)
#' res <- update_table(genes)
#' }
#' @export
process_genes <- function() {
  msg(UpdateAnno::hgncAlias2Symbol_version)

  genes <- UpdateAnno::hgncAlias2Symbol
  genes <- unique(genes[, .(SYMBOL, HGNC, ENTREZ)])
  data.table::setnames(genes, c("symbol", "hgnc", "entrez"))

  genes <- add_attributes(genes, "genes")

  genes
}
