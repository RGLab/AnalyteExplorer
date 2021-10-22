#' Create a metadata table for gene signatures
#'
#' Updates the gene symbols with Hugo, maps pathogen to the disease types
#' studied, and standardizes timepoint unit to days in gene signatures table
#' from HIPC Dashboard. This table is used to filter the gene signatures and
#' display metadata about the selected gene signature (or gene) in the app.
#'
#' @references
#' https://github.com/floratos-lab/hipc-dashboard-pipeline/tree/master/reformatted_data
#'
#' @return data.table
#' @examples
#' \dontrun{
#' library(UpdateAnno)
#' gene_signatures <- process_data("gene_signatures")
#' validate(gene_signatures)
#' res <- update_table(gene_signatures)
#' }
#' @import UpdateAnno
#' @export
process_gene_signatures <- function() {
  gene_signatures <- fetch_gene_signatures()

  gene_signatures[, id := row_key]

  gene_signatures[, disease_studied := sapply(target_pathogen, update_target_pathogen)]

  gene_signatures[, genes := sapply(response_component_original, update_gene_symbols)]
  gene_signatures <- gene_signatures[!is.na(genes), ]

  gene_signatures[, timepoint := mapply(update_timepoint, time_point, time_point_units)]
  gene_signatures[, timepoint_units := sapply(time_point_units, update_timepoint_units)]
  gene_signatures[, timepoint_with_units := paste(timepoint, timepoint_units, sep = "-")]

  columns_to_keep <- c(
    "id",
    "genes",
    "disease_studied",
    "response_behavior",
    "timepoint_with_units",
    "publication_reference_id",
    "comparison",
    "cohort"
  )
  gene_signatures <- gene_signatures[, ..columns_to_keep]

  gene_signatures <- add_attributes(gene_signatures, "gene_signatures")

  gene_signatures
}


# Helper Functions -------------------------------------------------------------

# Read in the latest gene signatures from HIPC Dashboard
fetch_gene_signatures <- function() {
  file <- tempfile()
  utils::download.file("https://github.com/floratos-lab/hipc-dashboard-pipeline/raw/master/reformatted_data/gene_expression-recreated_template.RDS", file)
  gene_signatures <- readRDS(file)

  gene_signatures <- gene_signatures[7:nrow(gene_signatures), 2:ncol(gene_signatures)]
  columns_to_keep <- c(
    "row_key",
    "response_component_original",
    "target_pathogen",
    "response_behavior",
    "time_point",
    "time_point_units",
    "publication_reference_id",
    "comparison",
    "cohort"
  )
  gene_signatures <- gene_signatures[, columns_to_keep]
  # remove one row of technical replicate - pub27974398
  gene_signatures <- gene_signatures[!duplicated(gene_signatures), ]

  data.table::setDT(gene_signatures)

  gene_signatures
}

# Map pathogen to the disease types given
update_target_pathogen <- function(x) {
  if (grepl("influenza", x, ignore.case = TRUE)) {
    return("Influenza")
  } else if (grepl("meningitidis", x, ignore.case = TRUE)) {
    return("Meningitidis")
  } else if (grepl("yellow fever", x, ignore.case = TRUE)) {
    return("Yellow Fever")
  } else if (grepl("ebola", x, ignore.case = TRUE)) {
    return("Ebola")
  } else if (grepl("immunodeficiency virus", x, ignore.case = TRUE)) {
    return("HIV")
  } else if (grepl("falciparum", x, ignore.case = TRUE)) {
    return("Malaria")
  } else if (grepl("mycobacterium", x, ignore.case = TRUE)) {
    return("Tuberculosis")
  } else if (grepl("vaccinia|variola", x, ignore.case = TRUE)) {
    return("Smallpox")
  } else if (grepl("rubella", x, ignore.case = TRUE)) {
    return("Rubella")
  } else if (grepl("pneumoniae", x, ignore.case = TRUE)) {
    return("Pneumonia")
  } else if (grepl("measles", x, ignore.case = TRUE)) {
    return("Measles")
  } else if (grepl("hepatitis B", x, ignore.case = TRUE)) {
    return("Hepatitis B")
  } else if (grepl("leishmania", x, ignore.case = TRUE)) {
    return("Leishmaniasis")
  } else if (grepl("papillomavirus", x, ignore.case = TRUE)) {
    return("HPV")
  } else if (grepl("tularensis", x, ignore.case = TRUE)) {
    return("Tularemia")
  } else if (grepl("equine", x, ignore.case = TRUE)) {
    return("VEEV")
  } else if (grepl("alphaherpesvirus", x, ignore.case = TRUE)) {
    return("Herpes")
  } else {
    warning(x, " is not recognized...")
    return(NA)
  }
}

# Ensure that gene symbols are current HUGO
update_gene_symbols <- function(x) {
  current <- strsplit(x, "; ")[[1]]
  symbols <- UpdateAnno::mapAlias2Symbol(current)
  symbols <- symbols[!is.na(symbols)]
  if (length(symbols) == 0) {
    return(NA)
  } else {
    return(paste(symbols, collapse = ", "))
  }
}

# Standardize timepoints to days
update_timepoint <- function(value, unit) {
  if (value %in% c("N/A", "")) {
    return(NA)
  } else if (grepl("\\d{1,2} (to|or) \\d{1,2}", value)) {
    vals <- strsplit(value, " (to|or) ")[[1]]
    return(round(stats::median(as.numeric(vals))))
  } else if (unit == "Weeks") {
    return(7 * as.numeric(value))
  } else if (unit == "Months") {
    return(30 * as.numeric(value))
  } else if (unit == "Hours") {
    return(round(as.numeric(value) / 24))
  } else {
    return(value)
  }
}

# Standardize timepoint units to days
update_timepoint_units <- function(unit) {
  if (grepl("day|hour", unit, ignore.case = TRUE)) {
    return("Days")
  } else if (unit == "" | unit == "N/A") {
    return("NA")
  } else if (unit %in% c("Months", "Weeks")) {
    return("Days")
  } else {
    return(unit)
  }
}
