#' Process gene expression data summarized by gene, blood transcription module,
#' and gene signature.
#'
#' @return data.table
#' @examples
#' \dontrun{
#' ge <- process_gene_expression()
#' check_table(ge)
#' res <- update_table("gene_expression", ge)
#' }
#' @export
process_gene_expression <- function() {
  log_message("Processing gene expression summary data...")

  metadata <- fetch_metadata()
  esets <- fetch_expression_sets(unique(metadata$study_accession))

  eset <- combine_expression_sets(esets)
  eset <- eset[, eset$participant_id %in% unique(metadata$biosample_participantid)]
  eset <- clean_expression_set(eset)
  eset <- add_metadata(eset)

  summary_by_gene <- summarize_expression_data(eset, "gene")
  summary_by_btm <- summarize_expression_data(eset, "blood transcription module")
  summary_by_gs <- summarize_expression_data(eset, "gene signature")

  # Combine summaries
  res <- rbind(summary_by_gene, summary_by_btm, summary_by_gs)
  res[is.infinite(mean_fold_change), mean_fold_change := NA]
  res$id <- seq_len(nrow(res))
  save_debug(res, "5_gene_expression")

  res
}


# HELPER FUNCTIONS -------------------------------------------------------------

# https://www.immunespace.org/query/Studies/executeQuery.view?schemaName=assay.ExpressionMatrix.matrix&query.queryName=InputSamples_computed
#' @importFrom data.table :=
fetch_metadata <- function() {
  log_message("Fetching gene expression metadata...")

  metadata <- Rlabkey::labkey.selectRows(
    baseUrl = get_url_base(),
    folderPath = "/Studies/",
    schemaName = "assay.ExpressionMatrix.matrix",
    queryName = "inputSamples_computed",
    colNameOpt = "rname",
    colSelect = c(
      "Biosample/ParticipantId",
      "Biosample/study_time_collected",
      "Biosample/study_time_collected_unit",
      "Run/DataOutputs/Name"
    ),
    colFilter = Rlabkey::makeFilter(c("Biosample/study_time_collected", "GREATER_THAN_OR_EQUAL", "0"))
  )
  data.table::setDT(metadata)

  metadata[, biosample_study_time_collected := mapply(
    convert_unit_to_day, biosample_study_time_collected, biosample_study_time_collected_unit
  )]
  metadata[, biosample_study_time_collected_unit := "Days"]

  metadata[, study_accession := gsub("SUB\\d{6}\\.", "SDY", biosample_participantid)]
  metadata[, expression_matrix := gsub("\\.tsv", "", run_dataoutputs_name)]

  metadata <- metadata[, has_data := 0 %in% unique(biosample_study_time_collected) &
                         length(unique(biosample_study_time_collected)) > 1,
                       by = .(biosample_participantid)
  ]
  metadata <- metadata[has_data == TRUE]

  metadata
}

fetch_expression_sets <- function(study_accessions) {
  log_message("Fetching and creating expression sets...")

  esets <- parallel::mclapply(study_accessions, function(study_accession) {
    log_message(study_accession)

    con <- ImmuneSpaceR::CreateConnection(study_accession)
    con$getGEMatrix(con$cache$GE_matrices$name)
  }, mc.cores = parallel::detectCores())
  save_debug(esets, "1_expressions_sets")

  esets
}

combine_expression_sets <- function(esets) {
  log_message("Combining expression sets...")

  eset <- ImmuneSpaceR:::.combineEMs(esets)
  save_debug(eset, "2a_expression_set")

  eset
}

clean_expression_set <- function(eset) {
  log_message("Cleaning combined expression set...")

  # Handle baseline dupes
  # remove values before dneg7
  eset <- eset[, eset$study_time_collected >= -7]

  # If pids have d0 and d7, remove dneg7.
  # If multiple baseline, select one
  baseline <- eset[, eset$study_time_collected <= 0]
  duplicates <- baseline$participant_id[duplicated(baseline$participant_id)]
  biosamples_to_remove <- sapply(duplicates, function(pid) {
    dup_entries <- baseline[, baseline$participant_id == pid]
    to_keep <- dup_entries$biosample_accession[dup_entries$study_time_collected == max(dup_entries$study_time_collected)][[1]]
    dup_entries$biosample_accession[dup_entries$biosample_accession != to_keep]
  })
  biosamples_to_remove <- unlist(unname(biosamples_to_remove))
  eset <- eset[, !eset$biosample_accession %in% biosamples_to_remove]

  save_debug(eset, "2b_expression_set_cleaned")

  eset
}

# Add metadata fields for hover text
add_metadata <- function(eset) {
  log_message("Adding metadata to combined expression set...")

  study_info <- Rlabkey::labkey.selectRows(
    baseUrl = get_url_base(),
    folderPath = "/Studies/",
    schemaName = "immport",
    queryName = "study",
    colNameOpt = "rname"
  )

  eset$study_accession <- gsub("SUB\\d{6}\\.", "SDY", eset$participant_id)
  split_cohort_type <- strsplit(eset$cohort_type, split = "_")
  eset$sample_type <- sapply(split_cohort_type, function(x) x[length(x)])
  eset$condition <- study_info$condition_studied[match(
    eset$study_accession,
    study_info$study_accession
  )]

  pd <- Biobase::pData(eset)
  pd <- map_condition(pd)
  if (all.equal(rownames(pd), colnames(Biobase::exprs(eset)))) {
    Biobase::pData(eset) <- pd
  } else {
    stop("Ensure ordering and matching of biosample IDs")
  }

  save_debug(eset, "2c_expression_set_with_metadata")

  eset
}

summarize_expression_data <- function(eset, by) {
  log_message(sprintf("Creating gene expression summary data by %s...", by))

  if (by == "blood transcription module") {
    blood_transcription_modules <- process_blood_transcription_modules()
    gene_sets <- blood_transcription_modules$genes
    names(gene_sets) <- blood_transcription_modules$id
  } else if (by == "gene signature") {
    gene_signatures <- process_gene_signatures()
    gene_sets <- gene_signatures$genes
    names(gene_sets) <- gene_signatures$id
  } else {
    gene_sets <- NULL
  }

  if (!is.null(gene_sets)) {
    em <- Biobase::exprs(eset)

    # summarize by gene set as average of genes included in gene sets
    em_list <- lapply(names(gene_sets), function(id) {
      gene_set <- gene_sets[id]
      gene_set <- strsplit(gene_set, split = ", ")[[1]]
      select_rows <- which(rownames(em) %in% gene_set)
      em_subset <- em[select_rows, ]
      if (!is.null(dim(em_subset))) {
        colMeans(em_subset, na.rm = TRUE)
      } else {
        log_message(sprintf("No genes selected for %s", id))
        em_subset
      }
    })

    em_by_gene_set <- data.frame(do.call(rbind, em_list), stringsAsFactors = FALSE)
    rownames(em_by_gene_set) <- names(gene_sets)

    eset <- Biobase::ExpressionSet(
      assayData = as.matrix(em_by_gene_set),
      phenoData = Biobase::AnnotatedDataFrame(Biobase::pData(eset))
    )
  }

  save_debug(eset, sprintf("3_expression_set_by_%s", by))

  create_summary_data(eset, by)
}

#' Convert ImmuneSpace condition-studied to a curated version
#'
#' @param pd phenotypic meta-data data.table with condition from ImmuneSpace
map_condition <- function(pd) {
  unmarked_studies <- list(
    influenza = c(
      "SDY80",
      "SDY144",
      "SDY224",
      "SDY296",
      "SDY301",
      "SDY364",
      "SDY368",
      "SDY387"
    ),
    hepatitis = c(
      "SDY89",
      "SDY299",
      "SDY690"
    ),
    smallpox = c(
      "SDY1370"
    ),
    herpes_zoster = c(
      "SDY984"
    )
  )

  pd$mapped_condition <- apply(pd, 1, function(x) {
    study <- x[["study_accession"]]
    condition <- x[["condition"]]
    cohort <- x[["cohort"]]

    if (study == "SDY180") {
      if (grepl("Saline", cohort)) {
        return("Healthy")
      } else if (grepl("Pneunomax23", cohort)) {
        return("Pneumonia")
      } else {
        return("Influenza")
      }
    } else if (study %in% unmarked_studies$influenza |
               grepl("influenza|H1N1", condition, ignore.case = TRUE)) {
      return("Influenza")
    } else if (study %in% unmarked_studies$hepatitis |
               grepl("Hepatitis", condition, ignore.case = TRUE)) {
      return("Hepatitis")
    } else if (study %in% unmarked_studies$smallpox |
               grepl("Smallpox|vaccinia", condition, ignore.case = TRUE)) {
      return("Smallpox")
    } else if (study %in% unmarked_studies$ppp) {
      return("Palmoplantar_Pustulosis")
    } else if (study %in% unmarked_studies$herpes_zoster) {
      return("Herpes_Zoster")
    } else if (grepl("healthy|normal|naive", condition, ignore.case = TRUE)) {
      return("Healthy")
    } else if (grepl("CMV", condition, ignore.case = TRUE)) {
      return("CMV")
    } else if (grepl("TB|tuberculosis", condition, ignore.case = TRUE)) {
      return("Tuberculosis")
    } else if (grepl("Yellow Fever", condition, ignore.case = TRUE)) {
      return("Yellow_Fever")
    } else if (grepl("Mening", condition, ignore.case = TRUE)) {
      return("Meningitis")
    } else if (grepl("Malaria", condition, ignore.case = TRUE)) {
      return("Malaria")
    } else if (grepl("HIV", condition, ignore.case = TRUE)) {
      return("HIV")
    } else if (grepl("Dengue", condition, ignore.case = TRUE)) {
      return("Dengue")
    } else if (grepl("ZEBOV", condition, ignore.case = TRUE)) {
      return("Ebola")
    } else if (grepl("JDM|Dermatomyositis", condition, ignore.case = TRUE)) {
      return("Dermatomyositis")
    } else if (grepl("West Nile", condition, ignore.case = TRUE)) {
      return("West_Nile")
    } else if (grepl("Zika", condition, ignore.case = TRUE)) {
      return("Zika")
    } else if (grepl("Varicella", condition, ignore.case = TRUE)) {
      return("Varicella_Zoster")
    } else {
      return("Unknown")
    }
  })

  pd
}

#' Create a summary data object for easy use with app
#'
#' @param eset eset object with cohort column and gene or btm expression
create_summary_data <- function(eset, analyte_type) {
  log_message("Creating summary data...")

  eset$study_cohort <- paste(eset$study_accession, eset$cohort, sep = "_")

  # Create data frame with summary statistics for each cohort*timepoint
  res <- lapply(unique(eset$study_cohort), function(study_cohort) {
    log_message(study_cohort)

    eset_subset <- eset[, eset$study_cohort == study_cohort]

    sample_type <- unique(eset_subset$sample_type)
    condition <- unique(eset_subset$mapped_condition)
    study_accession <- unique(eset_subset$study_accession)
    cohort <- unique(eset_subset$cohort)

    timepoints_table <- table(eset_subset$study_time_collected)
    timepoints <- as.numeric(names(timepoints_table)[timepoints_table > 2])

    baseline <- eset_subset[, eset_subset$study_time_collected == 0]

    subres <- lapply(timepoints, function(timepoint) {
      log_message(timepoint)

      if (timepoint == 0) {
        mean_fold_change <- 0
        sd_fold_change <- 0
        analyte_id <- rownames(baseline)
      } else {
        eset_post <- eset_subset[, eset_subset$study_time_collected == timepoint]

        sample_count <- table(eset_post$participant_id)
        duplicates <- names(sample_count)[sample_count > 1]

        if (length(duplicates) > 0) {
          log_message("Removing duplicates...")
          biosample_to_remove <- sapply(duplicates, function(pid) {
            dup_entries <- eset_post[, eset_post$participant_id == pid]
            max_day <- max(dup_entries$study_time_collected)
            to_keep <- dup_entries$biosample_accession[dup_entries$study_time_collected == max_day][[1]]
            curr_rm <- dup_entries$biosample_accession[dup_entries$biosample_accession != to_keep]
          })
          biosample_to_remove <- unlist(unname(biosample_to_remove))
          eset_post <- eset_post[, !eset_post$biosample_accession %in% biosample_to_remove]
        }

        shared <- intersect(baseline$participant_id, eset_post$participant_id)
        if (length(shared) < 3) {
          log_message("There are less than three shared participants between baseline and post...")
          return()
        }

        eset_post <- eset_post[, eset_post$participant_id %in% shared]
        eset_base <- baseline[, baseline$participant_id %in% shared]
        eset_post <- eset_post[, order(match(eset_post$participant_id, eset_base$participant_id))]

        em_post <- Biobase::exprs(eset_post)
        em_base <- Biobase::exprs(eset_base)
        em_post <- em_post[order(match(rownames(em_post), rownames(em_base))), ]

        if (!all.equal(dim(em_base), dim(em_base))) {
          stop("Dimensions do not match...")
        }

        fold_change <- em_post - em_base

        mean_fold_change <- rowMeans(fold_change, na.rm = TRUE)
        sd_fold_change <- apply(fold_change, 1, function(x) stats::sd(x, na.rm = TRUE))
        analyte_id <- rownames(fold_change)
      }

      data.frame(
        cohort = cohort,
        sample_type = sample_type,
        study_accession = study_accession,
        condition = condition,
        timepoint = timepoint,
        analyte_id = analyte_id,
        analyte_type = analyte_type,
        mean_fold_change = mean_fold_change,
        sd_fold_change = sd_fold_change,
        stringsAsFactors = FALSE
      )
    })

    do.call("rbind", subres)
  })

  summarized <- data.table::rbindlist(res)
  save_debug(summarized, sprintf("4_summarized_by_%s", analyte_type))

  summarized
}

convert_unit_to_day <- function(value, unit) {
  unit <- tolower(unit)
  if (unit %in% c("day", "days")) {
    return(value)
  } else if (unit %in% c("hour", "hours")) {
    return(value / 24)
  } else if (unit %in% c("week", "weeks")) {
    return(value * 7)
  } else if (unit %in% c("month", "months")) {
    return(value * 30)
  } else if (unit %in% c("year", "years")) {
    return(value * 365)
  } else {
    log_message("")
    return(NA)
  }
}
