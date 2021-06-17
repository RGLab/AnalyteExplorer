shinyUI(fluidPage(
  pageWithSidebar(

    # Application title
    headerPanel("Analyte Explorer"),

    # Sidebar with four inputs
    sidebarPanel(
      h3("Overview"),
      p("This analysis allows the user to visualize log fold change of a gene or gene set over time by cohort. In the case of genes and blood transcription modules, the fold change is the delta of log2-transformed expression values, while gene signatures uses fold change of the geometric mean. A trendline is drawn for points with sufficient data.
              For all cohorts other than 'Healthy', the exposure process was 'vaccination'."),
      selectizeInput(
        inputId = "condition_studied",
        label = "Select Conditions Studied for Plotting:",
        choices = list(
          "Influenza" = "Influenza",
          "Meningitis" = "Meningitis",
          "Herpes Zoster" = "Herpes_Zoster",
          "Yellow Fever" = "Yellow_Fever",
          "Malaria" = "Malaria",
          "Tuberculosis" = "Tuberculosis",
          "Ebola" = "Ebola",
          "HIV" = "HIV",
          "Palmoplantar Pustulosis" = "Palmoplantar_Pustulosis",
          "Dermatomyositis" = "Dermatomyositis",
          "Hepatitis" = "Hepatitis",
          "Varicella Zoster" = "Varicella_Zoster",
          "Smallpox" = "Smallpox",
          "Healthy" = "Healthy",
          "Pneumonia" = "Pneumonia"
        ),
        options = list(
          maxItems = 5,
          placeholder = "Select Conditions Studied"
        )
      ),
      radioButtons(
        inputId = "analyte_type",
        label = "Select type of feature:",
        choices = list(
          "Gene",
          "Blood Transcription Module",
          "Gene Signature"
        ),
        selected = "Gene"
      ),
      p(),
      conditionalPanel(
        h4("Gene Signature Filters:"),
        condition = "input.analyte_type == 'Gene Signature'",
        selectizeInput(
          inputId = "gs_disease_studied",
          label = "Disease Studied:",
          choices = NULL,
          options = list(
            maxItems = 1,
            placeholder = "Select Disease Studied"
          )
        ),
        selectizeInput(
          inputId = "gs_timepoint_with_units",
          label = "Timepoint",
          choices = NULL,
          options = list(
            maxItems = 1,
            placeholder = "Select Timepoints"
          )
        ),
        selectizeInput(
          inputId = "gs_response_behavior",
          label = "Response Behavior",
          choices = NULL,
          options = list(
            maxItems = 1,
            placeholder = "Select Response Behavior"
          )
        ),
        actionButton("apply_filters", "Apply Filters"),
        actionButton("reset_filters", "Reset Filters"),
        p("Note: Apply filters one at a time to update dropdowns")
      ),
      strong("Selected Gene or Gene Set"),
      selectizeInput(
        inputId = "analyte_selection",
        label = NULL,
        choices = NULL,
        options = list(
          maxItems = 1,
          placeholder = "Select Gene or Gene Set"
        )
      ),
      actionButton("submit", "Submit"),
      h3(),
      p("Shiny app development by Evan Henrich @ Fred Hutch.")
    ),

    # OUTPUT
    mainPanel(
      tabsetPanel(
        type = "tabs",
        tabPanel("Plots", plotly::plotlyOutput("line_plots", height = "700px")),
        tabPanel("Metadata", DT::dataTableOutput("metadata"))
      )
    )
  )
))
