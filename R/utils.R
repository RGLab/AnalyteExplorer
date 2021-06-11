get_url_base <- function() {
  mget("labkey.url.base", envir = .GlobalEnv, ifnotfound = "https://www.immunespace.org")[[1]]
}

get_url_path <- function() {
  mget("labkey.url.path", envir = .GlobalEnv, ifnotfound = "/AnalyteExplorer/")[[1]]
}

log_message <- function(msg) {
  message(sprintf("[%s] %s", Sys.time(), msg))
}

labkey.importData <- function(baseUrl, folderPath, schemaName, queryName, toImport) {
  temp_file <- tempfile(fileext = ".csv")
  data.table::fwrite(toImport, temp_file)

  url <- paste(baseUrl, "query", folderPath, "import.api", sep = "")
  config <- Rlabkey:::labkey.getRequestOptions(method = "POST")
  body <- list(
    schemaName = schemaName,
    queryName = queryName,
    file = httr::upload_file(temp_file)
  )
  header <- httr::add_headers(`Content-Type` = "multipart/form-data")

  response <- httr::POST(url = url, config = config, body = body, header)
  parsed <- Rlabkey:::processResponse(response, haltOnError = TRUE)

  parsed
}

fix_types <- function(types) {
  sapply(types, function(x) {
    switch(
      x,
      "Integer" = "numeric",
      "Number (Double)" = "numeric",
      "Text (String)" = "character",
      stop(sprintf("%s is an unknown type!", x))
    )
  }, USE.NAMES = FALSE)
}

save_debug <- function(object, name) {
  debug_dir <- getOption("debug_dir")
  if (!is.null(debug_dir)) {
    file <- sprintf("%s/%s.rds", debug_dir, name)
    log_message(sprintf("Saving to %s...", file))
    saveRDS(object, file)
  }
}

