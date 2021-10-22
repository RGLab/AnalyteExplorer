#' Process gene expression data summarized by gene, blood transcription module,
#' and gene signature.
#'
#' @return data.table
#' @examples
#' \dontrun{
#' summaries <- process_data("gene_expression_summaries")
#' validate(summaries)
#' res <- update_table(summaries)
#' }
#' @export
process_gene_expression_summaries <- function() {
  msg("Processing gene expression summary data")

  metadata <- fetch_metadata()
  esets <- fetch_expression_sets(unique(metadata$study_accession))

  eset <- combine_expression_sets(esets)
  eset <- eset[, eset$participant_id %in% unique(metadata$biosample_participantid)]
  eset <- clean_expression_set(eset)
  eset <- add_metadata(eset, metadata)

  summary_by_gene <- summarize_expression_data(eset, "gene")
  summary_by_btm <- summarize_expression_data(eset, "blood transcription module")
  summary_by_gs <- summarize_expression_data(eset, "gene signature")

  # Combine summaries
  res <- rbind(summary_by_gene, summary_by_btm, summary_by_gs)
  res[is.infinite(mean_fold_change), mean_fold_change := NA]
  res$id <- seq_len(nrow(res))

  res <- add_attributes(res, "gene_expression_summaries")

  save_debug(res, "5_gene_expression_summaries")

  res
}


# HELPER FUNCTIONS -------------------------------------------------------------

# https://www.immunespace.org/query/Studies/executeQuery.view?schemaName=assay.ExpressionMatrix.matrix&query.queryName=InputSamples_computed
#' @importFrom data.table :=
fetch_metadata <- function() {
  msg("Fetching gene expression metadata")

  metadata <- Rlabkey::labkey.selectRows(
    baseUrl = get_url_base(),
    folderPath = "/Studies/",
    schemaName = "assay.ExpressionMatrix.matrix",
    queryName = "inputSamples_computed",
    colNameOpt = "rname",
    colSelect = c(
      "Biosample/biosample_accession",
      "Biosample/ParticipantId",
      "Biosample/study_time_collected",
      "Biosample/study_time_collected_unit",
      "Run/DataOutputs/Name",
      "Biosample/arm_accession",
      "Biosample/type"
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
  msg("Fetching and creating expression sets")

  esets <- lapply(study_accessions, function(study_accession) {
    msg(study_accession)

    con <- ImmuneSpaceR::CreateConnection(study_accession)
    con$getGEMatrix(con$cache$GE_matrices$name)
  })

  # if (any(sapply(esets, function(x) is(x, "try-error")))) {
  #   stop("Something went wrong!")
  # }

  save_debug(esets, "1_expressions_sets")

  esets
}

combine_expression_sets <- function(esets) {
  msg("Combining expression sets")

  eset <- ImmuneSpaceR:::.combineEMs(esets)
  save_debug(eset, "2a_expression_set")

  eset
}

clean_expression_set <- function(eset) {
  msg("Cleaning combined expression set")

  # Handle baseline dupes
  # remove values before dneg7
  eset <- eset[, eset$study_time_collected >= -7]

  # If pids have d0 and d7, remove dneg7.
  # If multiple baseline, select one
  baseline <- eset[, eset$study_time_collected <= 0]
  duplicates <- baseline$participant_id[duplicated(baseline$participant_id)]
  msg(sprintf("Removing %s duplicates", length(duplicates)))
  biosamples_to_remove <- sapply(duplicates, function(pid) {
    msg(pid)
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
add_metadata <- function(eset, metadata) {
  msg("Adding metadata to combined expression set")

  metadata <- unique(metadata[, .(biosample_biosample_accession, biosample_arm_accession, biosample_type, study_accession)])
  data.table::setnames(metadata, c("biosample_accession", "arm_accession", "sample_type", "study_accession"))

  pd <- Biobase::pData(eset)
  data.table::setDT(pd)
  pd <- merge(pd, metadata, by = "biosample_accession")
  rownames(pd) <- pd$biosample_accession
  Biobase::pData(eset) <- pd
  Biobase::sampleNames(eset) <- eset$biosample_accession

  save_debug(eset, "2c_expression_set_with_metadata")

  eset
}

summarize_expression_data <- function(eset, by) {
  msg(sprintf("Creating gene expression summary data by %s", by))

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
      # intersect
      genes_in_em <- intersect(rownames(em), gene_set)
      n_genes_in_em <- length(genes_in_em)
      n_gene_set <- length(gene_set)

      msg(id)
      if (n_genes_in_em > 0) {
        if (n_genes_in_em != n_gene_set) {
          msg(sprintf("%s/%s genes in matrix", n_genes_in_em, n_gene_set))
        }
        if (n_genes_in_em == 1) {
          em[genes_in_em, ]
        } else {
          colMeans(em[genes_in_em, ], na.rm = TRUE)
        }
      } else {
        msg("No genes selected")
        NULL
      }
    })

    em_by_gene_set <- do.call(rbind, em_list)
    rownames(em_by_gene_set) <- names(gene_sets)[!sapply(em_list, is.null)]

    eset <- Biobase::ExpressionSet(
      assayData = em_by_gene_set[, eset$biosample_accession],
      phenoData = Biobase::AnnotatedDataFrame(Biobase::pData(eset))
    )
  }

  save_debug(eset, sprintf("3_expression_set_by_%s", by))

  create_summary_data(eset, by)
}

# Create a summary data object for easy use with app
create_summary_data <- function(eset, analyte_type) {
  msg("Creating summary data")

  eset$study_cohort_type <- paste(eset$study_accession, eset$arm_accession, eset$sample_type, sep = "_")

  # Create data frame with summary statistics for each cohort*timepoint
  res <- lapply(unique(eset$study_cohort_type), function(study_cohort_type) {
    msg(study_cohort_type)

    eset_subset <- eset[, eset$study_cohort_type == study_cohort_type]

    sample_type <- unique(eset_subset$sample_type)
    study_accession <- unique(eset_subset$study_accession)
    arm_accession <- unique(eset_subset$arm_accession)

    timepoints_table <- table(eset_subset$study_time_collected)
    timepoints <- as.numeric(names(timepoints_table)[timepoints_table > 2])

    baseline <- eset_subset[, eset_subset$study_time_collected == 0]

    subres <- lapply(timepoints, function(timepoint) {
      msg(sprintf("Day %s", timepoint))

      if (timepoint == 0) {
        mean_fold_change <- 0
        sd_fold_change <- 0
        analyte_id <- rownames(baseline)
      } else {
        eset_post <- eset_subset[, eset_subset$study_time_collected == timepoint]

        sample_count <- table(eset_post$participant_id)
        duplicates <- names(sample_count)[sample_count > 1]

        if (length(duplicates) > 0) {
          msg(sprintf("Removing %s duplicates", length(duplicates)))
          biosample_to_remove <- sapply(duplicates, function(pid) {
            msg(pid)
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
          msg("There are less than three shared participants between baseline and post")
          return()
        }

        eset_post <- eset_post[, eset_post$participant_id %in% shared]
        eset_base <- baseline[, baseline$participant_id %in% shared]
        eset_post <- eset_post[, order(match(eset_post$participant_id, eset_base$participant_id))]

        em_post <- Biobase::exprs(eset_post)
        em_base <- Biobase::exprs(eset_base)
        em_post <- em_post[order(match(rownames(em_post), rownames(em_base))), ]

        if (!all.equal(dim(em_base), dim(em_base))) {
          stop("Dimensions do not match")
        }

        fold_change <- em_post - em_base

        mean_fold_change <- rowMeans(fold_change, na.rm = TRUE)
        sd_fold_change <- apply(fold_change, 1, function(x) stats::sd(x, na.rm = TRUE))
        analyte_id <- rownames(fold_change)
      }

      data.frame(
        arm_accession = arm_accession,
        sample_type = sample_type,
        study_accession = study_accession,
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
    msg("")
    return(NA)
  }
}
