#' Create a metadata table for cohorts
#'
#' @return data.table
#' @examples
#' \dontrun{
#' cohorts <- process_data("cohorts")
#' validate(cohorts)
#' res <- update_table(cohorts)
#' }
#' @export
process_cohorts <- function() {
  cohorts <- Rlabkey::labkey.selectRows(
    baseUrl = get_url_base(),
    folderPath = "/Studies",
    schemaName = "immport",
    queryName = "arm_or_cohort",
    colNameOpt = "rname",
    colSelect = c("arm_accession", "description", "name", "study_accession")
  )
  data.table::setDT(cohorts)

  studies <- unique(Rlabkey::labkey.selectRows(
    baseUrl = get_url_base(),
    folderPath = "/Studies",
    schema = "assay.ExpressionMatrix.matrix",
    query = "SelectedRuns",
    colNameOpt = "rname",
    colSelect = "study"
  )$study)

  study_info <- Rlabkey::labkey.selectRows(
    baseUrl = get_url_base(),
    folderPath = "/Studies/",
    schemaName = "immport",
    queryName = "study",
    colNameOpt = "rname",
    colSelect = c("study_accession", "condition_studied")
  )
  data.table::setDT(study_info)

  research_focus <- Rlabkey::labkey.selectRows(
    baseUrl = get_url_base(),
    folderPath = "/Studies",
    schemaName = "immport",
    queryName = "study_categorization",
    colNameOpt = "rname",
    colSelect = c("study_accession", "research_focus")
  )
  data.table::setDT(research_focus)

  cohorts <- cohorts[study_accession %in% studies]
  cohorts <- merge(cohorts, research_focus, by = "study_accession", all.x = TRUE)
  cohorts[, research_focus := ifelse(is.na(research_focus), "Uknown", research_focus)]
  cohorts <- merge(cohorts, study_info, by = "study_accession", all.x = TRUE)
  cohorts <- map_condition(cohorts)

  cohorts <- add_attributes(cohorts, "cohorts")

  cohorts
}

# Convert ImmPort condition_studied column to a curated version
map_condition <- function(pd) {
  unmarked_studies <- list(
    influenza = c(
      "SDY80",
      "SDY144",
      "SDY212",
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

  pd$condition_studied <- apply(pd, 1, function(x) {
    study <- x[["study_accession"]]
    condition <- x[["condition_studied"]]
    cohort <- x[["name"]]

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
      return("Palmoplantar Pustulosis")
    } else if (study %in% unmarked_studies$herpes_zoster) {
      return("Herpes Zoster")
    } else if (grepl("healthy|normal|naive", condition, ignore.case = TRUE)) {
      return("Healthy")
    } else if (grepl("CMV", condition, ignore.case = TRUE)) {
      return("Cytomegalovirus")
    } else if (grepl("TB|tuberculosis", condition, ignore.case = TRUE)) {
      return("Tuberculosis")
    } else if (grepl("Yellow Fever", condition, ignore.case = TRUE)) {
      return("Yellow Fever")
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
      return("West Nile")
    } else if (grepl("Zika", condition, ignore.case = TRUE)) {
      return("Zika")
    } else if (grepl("Varicella", condition, ignore.case = TRUE)) {
      return("Varicella Zoster")
    } else if (grepl("Pertussis", condition, ignore.case = TRUE)) {
      return("Pertussis")
    } else {
      return("Unknown")
    }
  })

  pd
}
