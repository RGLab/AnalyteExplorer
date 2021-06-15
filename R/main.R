#' Process data
#'
#' @param data_name a character. data name to process: blood_transcript_modules,
#' gene_signatures, or gene_expression
#'
#' @return A data.table.
#' @examples
#' \dontrun{
#' data <- process_data("gene_expression")
#' }
#' @export
process_data <- function(data_name) {
  log_message(sprintf("Processing '%s'...", data_name))
  switch(
    data_name,
    "blood_transcript_modules" = process_blood_transcript_modules(),
    "gene_signatures" = process_gene_signatures(),
    "gene_expression" = process_gene_expression(),
    stop(data_name, " is not a valid data name...")
  )
}


#' Check data
#'
#' @param data_name a character. data name to process: blood_transcript_modules,
#' gene_signatures, or gene_expression
#' @param data a data.table.
#'
#' @return A logical.
#' @examples
#' \dontrun{
#' check_table("blood_transcription_modules", data)
#' check_table("gene_signatures", data)
#' check_table("gene_expression", data)
#' }
#' @export
check_table <- function(data_name, data) {
  log_message("Checking data...")

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
#' @param data_name a character. data name to process: blood_transcript_modules,
#' gene_expression, gene_signatures
#' @param data a data.table.
#'
#' @return A list.
#' @examples
#' \dontrun{
#' res <- update_table("gene_expression", data)
#' }
#' @export
update_table <- function(data_name, data) {
  schema_name <- "analyte_explorer"

  log_message("Backing up...")
  backup <- suppressWarnings(Rlabkey::labkey.selectRows(
    baseUrl = get_url_base(),
    folderPath = get_url_path(),
    schemaName = schema_name,
    queryName = data_name,
    colNameOpt = "rname"
  ))

  if (nrow(backup) > 0) {
    log_message("Truncating table...")
    Rlabkey::labkey.truncateTable(
      baseUrl = get_url_base(),
      folderPath = get_url_path(),
      schemaName = schema_name,
      queryName = data_name
    )
  }

  res <- tryCatch(
    {
      log_message(sprintf("Importing %s rows of data...", nrow(data)))
      labkey.importData(
        baseUrl = get_url_base(),
        folderPath = get_url_path(),
        schemaName = schema_name,
        queryName = data_name,
        toImport = data
      )
    },
    error = function(e) {
      log_message("Failed to import data!")
      if (nrow(backup) > 0) {
        log_message("Restoring table...")
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
      log_message("Import completed...")
    }
  )

  res
}
