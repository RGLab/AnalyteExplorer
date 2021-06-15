#' Create a metadata table for blood transcription modules
#'
#' Updates the gene symbols with Hugo in blood transcription modules table from
#' Li S, et al. This table is used to display metadata about the selected module
#' in the app app.
#'
#' @references
#' Li S, Rouphael N, Duraisingham S, et al.
#' Molecular signatures of antibody responses derived from a systems biology
#' study of five human vaccines.
#' Nat Immunol. 2014;15(2):195-204. doi:10.1038/ni.2789
#' \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3946932/bin/NIHMS540680-supplement-26.zip}
#'
#' @return data.table
#' @examples
#' \dontrun{
#' btm <- process_blood_transcription_modules()
#' check_table(btm)
#' res <- update_table("blood_transcript_modules", btm)
#' }
#' @export
process_blood_transcription_modules <- function() {
  file <- system.file("extdata/btm_annotation_table.csv", package = "UpdateAnno", mustWork = TRUE)
  metadata <- data.table::fread(file)
  btms <- UpdateAnno::emory_blood_transcript_modules

  metadata[, current_genes := sapply(ID, function(x) {
    i <- grep(sprintf(" \\(%s\\)$", x), names(btms))[1]
    paste(btms[[i]], collapse = ", ")
  })]
  columns_to_keep <- c(
    "ID",
    "Composite name",
    "current_genes",
    "Top matched Gene Ontology terms (number of matched genes)",
    "Module size",
    "Module category"
  )
  metadata <- metadata[, ..columns_to_keep]
  new_names <- c(
    "id",
    "name",
    "genes",
    "matched_gene_ontology_terms",
    "number_of_genes",
    "module_category"
  )
  data.table::setnames(metadata, colnames(metadata), new_names)

  metadata
}
