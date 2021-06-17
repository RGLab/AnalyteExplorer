library(plotly)
library(DT)
library(RColorBrewer)

# Helpers
get_figure_list <- function(pd) {
  selected_conditions <- unique(pd$condition)
  selected_analyte <- unique(pd$analyte_id)

  fig_list <- list()

  # color-palette
  pal <- RColorBrewer::brewer.pal(n = 7, "Dark2")
  pal <- sample(pal, length(selected_conditions))

  # Axes
  day_range <- range(pd$timepoint)
  fc_range <- range(pd$mean_fold_change)

  for (i in 1:length(selected_conditions)) {
    condition <- as.character(selected_conditions[[i]])
    pd_subset <- pd[pd$condition == condition, ]

    # Trendline should only be for points at which there is sufficient data
    min_cohorts_for_trend <- 0.4
    timepoints <- table(pd_subset$timepoint)
    timepoints <- timepoints[timepoints > min_cohorts_for_trend * length(unique(pd_subset$cohort))]
    timepoints <- sort(names(timepoints))
    mean_values <- sapply(timepoints, function(pt) {
      pd_pt <- pd_subset[pd_subset$timepoint == pt]
      mean(pd_pt$mean_fold_change)
    })

    pd_trend <- data.frame(
      cohort = "Average",
      cell_type = "NA",
      study = "Trend",
      analyte = selected_analyte,
      condition = "Trend",
      timepoint = as.double(timepoints),
      mean_fold_change = mean_values,
      sd_fold_change = 0,
      stringsAsFactors = FALSE
    )

    pd_subset <- merge(pd_subset, pd_trend, all = TRUE)
    pd_subset <- pd_subset[order(pd_subset$study_accession), ] # must draw trend at end to visible
    cohorts <- unique(pd_subset$cohort)

    p <- plot_ly()
    color_map <- c(Trend = "#5b5c5b")
    color_map[[condition]] <- pal[[i]]

    for (selected_cohort in cohorts) {
      p <- add_trace(p,
        data = pd_subset[which(pd_subset$cohort == selected_cohort), ],
        x = ~timepoint,
        y = ~mean_fold_change,
        color = ~condition,
        colors = color_map,
        text = selected_cohort,
        customdata = ~study_accession,
        hovertemplate = paste(
          "<b>Cohort</b>: %{text}",
          "<br><b>Study</b>: %{customdata}",
          "<br><b>Timepoint</b>: %{x}",
          "<br><b>log2-FC</b>: %{y:.2f}",
          "<extra></extra>"
        ),
        type = "scatter",
        mode = "lines+markers"
      )
    }

    fig_list[[condition]] <- p %>% layout(
      showlegend = FALSE,
      yaxis = list(
        title = condition,
        range = fc_range
      ),
      xaxis = list(range = day_range)
    )
  }

  fig <- subplot(fig_list,
    shareY = TRUE,
    nrows = length(selected_conditions),
    titleX = FALSE
  )

  return(fig)
}

shinyServer(function(input, output, session) {

  #---------------------------------------------------------
  #                       FOR TESTING
  #---------------------------------------------------------

  # input <- list(
  #   analyte_type = "Gene Signature",
  #   gs_disease_studied = "Herpes",
  #   gs_timepoint_with_units = "ALL",
  #   gs_response_behavior = "ALL",
  #   analyte_selection = "Early patterns of gene expression correlate with the humoral immune response to influenza vaccination in humans",
  #   condition = "Influenza"
  # )


  #---------------------------------------------------------
  #                       MAIN
  #---------------------------------------------------------
  # STATE
  filters_stored <- list(
    disease_studied = "ALL",
    timepoint_with_units = "ALL",
    response_behavior = "ALL"
  )

  # INIT
  for (name in names(filters_stored)) {
    ui_element <- paste0("gs_", name)
    updateSelectizeInput(session,
      ui_element,
      choices = c("ALL", unique(gs[[name]])),
      server = TRUE
    )
  }


  # On Analyte Type Selection
  observeEvent(input$analyte_type, {
    if (input$analyte_type == "Gene") {
      options <- unique(ge_gene$analyte_id)
    } else if (input$analyte_type == "Blood Transcription Module") {
      options <- btm$id
    } else if (input$analyte_type == "Gene Signature") {
      options <- gs$id
    }

    updateSelectizeInput(session,
      "analyte_selection",
      choices = options,
      server = TRUE
    )
  })

  # GeneSignatures Conditional Filters
  observeEvent(input$apply_filters, {
    filters_current <- list(
      disease_studied = input$gs_disease_studied,
      timepoint_with_units = input$gs_timepoint_with_units,
      response_behavior = input$gs_response_behavior
    )

    filters_unchanged <- filters_current[unlist(filters_current) == unlist(filters_stored)]

    gs_current <- gs

    for (fn in names(filters_current)) {
      if (filters_current[fn] != "ALL") {
        gs_current <- gs_current[gs_current[[fn]] == filters_current[[fn]], ]
      }
    }

    # Update all unchanged filters
    for (fn in names(filters_unchanged)) {
      ui_element <- paste0("gs_", fn)
      if (nrow(gs_current) > 0) {
        choices <- c("ALL", unique(gs_current[[fn]]))
      } else {
        choices <- "No Choices - Reset Filters!"
      }
      updateSelectizeInput(session,
        ui_element,
        choices = choices,
        server = TRUE
      )
    }

    # Update the analyteSelection
    if (nrow(gs_current) > 0) {
      choices <- unique(gs_current$uid)
    } else {
      choices <- "No Choices - Reset Filters!"
    }
    updateSelectizeInput(session,
      "analyte_selection",
      choices = choices,
      server = TRUE
    )
  })

  observeEvent(input$reset_filters, {
    filters_stored <- list(
      disease_studied = "ALL",
      timepoint_with_units = "ALL",
      response_behavior = "ALL"
    )

    for (nm in names(filters_stored)) {
      ui_element <- paste0("gs_", nm)
      updateSelectizeInput(session,
        ui_element,
        choices = c("ALL", unique(gs[[nm]])),
        server = TRUE
      )
    }

    updateSelectizeInput(session,
      "analyte_selection",
      choices = unique(gs$id),
      server = TRUE
    )
  })

  # Generate plot
  plot_data <- reactiveValues(data = NULL)
  metadata <- reactiveValues(data = NULL)

  observeEvent(input$submit, {
    pubmed_url <- "https://pubmed.ncbi.nlm.nih.gov/"

    if (input$analyte_type == "Gene") {
      data <- ge_gene

      selected_gene <- paste0("(;|)", input$analyte_selection[1], "(;|)")
      print(input$analyte_selection)
      print(selected_gene)
      selected <- grepl(selected_gene, gs$genes)
      print(selected)

      if (sum(selected) > 0) {
        info <- gs[selected, ]
        print(info)
        info$link <- paste0(
          '<a href="',
          paste0(pubmed_url, info$publication_reference),
          '">',
          info$publication_reference,
          "</a>"
        )
      } else {
        info <- data.frame()
      }

      metadata_options <- list(pageLength = 5)
    } else if (input$analyte_type == "Blood Transcription Module") {
      data <- ge_btm

      info <- btm[btm$id == input$analyte_selection[1], ]
      metadata_options <- list(
        pageLength = 1,
        searching = FALSE,
        paging = FALSE,
        info = FALSE,
        ordering = FALSE
      )
    } else if (input$analyte_type == "Gene Signature") {
      data <- ge_gs
      info <- gs[gs$id == input$analyte_selection[1], ]
      info$link <- paste0(
        '<a href="',
        paste0(pubmed_url, info$publication_reference),
        '">',
        info$publication_reference,
        "</a>"
      )
      metadata_options <- list(
        pageLength = 1,
        searching = FALSE,
        paging = FALSE,
        info = FALSE,
        ordering = FALSE
      )
    }

    plot_data$data <- data[data$analyte_id == input$analyte_selection &
      data$condition %in% input$condition_studied, ]

    metadata$data <- datatable(info,
      options = metadata_options,
      rownames = FALSE,
      escape = FALSE
    )
  })


  #---------------------------------------------------------
  #                       OUTPUTS
  #---------------------------------------------------------

  output$line_plots <- plotly::renderPlotly({
    if (is.null(plot_data$data)) {
      return()
    }

    fig <- get_figure_list(plot_data$data)
  })

  output$metadata <- renderDataTable(metadata$data)
})
