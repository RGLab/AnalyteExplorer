#' Process data
#'
#' @param data_name a character. data name to process: blood_transcript_modules,
#' gene_signatures, or gene_expression_summaries
#'
#' @return A data.table.
#' @examples
#' \dontrun{
#' data <- process_data("gene_expression_summaries")
#' }
#' @export
process_data <- function(data_name) {
  msg(sprintf("Processing '%s'", data_name))
  switch(
    data_name,
    "genes" = process_genes(),
    "blood_transcript_modules" = process_blood_transcript_modules(),
    "gene_signatures" = process_gene_signatures(),
    "gene_expression_summaries" = process_gene_expression_summaries(),
    "cohorts" = process_cohorts(),
    stop(data_name, " is not a valid data name")
  )
}


#' Validate data
#'
#' @param data a data.table.
#' @param schema_name a character.
#'
#' @return A logical.
#' @examples
#' \dontrun{
#' validate(data)
#' }
#' @export
validate <- function(data, schema_name = "analyte_explorer") {
  assertthat::assert_that(methods::is(data, "AnalyteExplorer"))
  data_name <- attr(data, "type")

  msg("Validating data")

  assertthat::assert_that(nrow(data) > 0)
  assertthat::assert_that(ncol(data) > 0)

  # check for NULLs
  is_na <- sapply(data, function(x) any(is.na(x)))
  if (any(is_na)) {
    columns <- paste0(names(is_na)[is_na], collapse = "', '")
    stop(sprintf("'%s' columns have null or empty values!", columns))
  }

  # check for valid column names
  colnames_data <- sort(colnames(data))
  empty_table <- suppressWarnings(Rlabkey::labkey.selectRows(
    baseUrl = get_url_base(),
    folderPath = get_url_path(),
    schemaName = "analyte_explorer",
    queryName = data_name,
    maxRows = 0,
    colNameOpt = "rname"
  ))
  colnames_server <- sort(colnames(empty_table))
  assertthat::assert_that(length(colnames_data) == length(colnames_server))
  assertthat::assert_that(identical(colnames_data, colnames_server))

  # check for valid column types
  query_details <- Rlabkey::labkey.getQueryDetails(
    baseUrl = get_url_base(),
    folderPath = get_url_path(),
    schemaName = "analyte_explorer",
    queryName = data_name
  )
  coltypes_server <- fix_types(query_details$type)
  names(coltypes_server) <- query_details$fieldName
  coltypes_server <- coltypes_server[colnames_data]
  coltypes_data <- sapply(data, mode)
  coltypes_data <- coltypes_data[colnames_data]
  assertthat::assert_that(identical(coltypes_data, coltypes_server))

  invisible(TRUE)
}


#' Update table
#'
#' @param data a data.table.
#' @param schema_name a character.
#'
#' @return A list.
#' @examples
#' \dontrun{
#' res <- update_table("gene_expression", data)
#' }
#' @export
update_table <- function(data, schema_name = "analyte_explorer") {
  assertthat::assert_that(methods::is(data, "AnalyteExplorer"))
  data_name <- attr(data, "type")

  msg("Backing up")
  backup <- suppressWarnings(Rlabkey::labkey.selectRows(
    baseUrl = get_url_base(),
    folderPath = get_url_path(),
    schemaName = schema_name,
    queryName = data_name,
    colNameOpt = "rname"
  ))

  if (nrow(backup) > 0) {
    msg("Truncating table")
    Rlabkey::labkey.truncateTable(
      baseUrl = get_url_base(),
      folderPath = get_url_path(),
      schemaName = schema_name,
      queryName = data_name
    )
  }

  res <- tryCatch(
    {
      msg(sprintf("Importing %s rows of data", nrow(data)))
      labkey.importData(
        baseUrl = get_url_base(),
        folderPath = get_url_path(),
        schemaName = schema_name,
        queryName = data_name,
        toImport = data
      )
    },
    error = function(e) {
      msg("Failed to import data!")
      if (nrow(backup) > 0) {
        msg("Restoring table")
        labkey.importData(
          baseUrl = get_url_base(),
          folderPath = get_url_path(),
          schemaName = schema_name,
          queryName = data_name,
          toImport = backup
        )
      }
      stop(e)
    },
    finally = {
      msg("Import completed")
    }
  )

  res
}
