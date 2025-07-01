# GenTile Shiny App
# Load required libraries with error handling
required_packages <- c("shiny", "shinydashboard", "DT", "plotly", "shinyWidgets", "rlang")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("Package", pkg, "is required but not installed. Please install it with: install.packages('", pkg, "')", sep = ""))
  }
}

library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(shinyWidgets)
library(rlang)  # For %||% operator

# Set genTile repository path (DEFINE BEFORE SOURCING)
GENTILE_PATH <- "/Users/nathanaelandrews/wrk/github/genTile"

# Source helper functions
source("functions/pipeline_runner.R")
source("functions/plotting_functions.R")

# UI
ui <- dashboardPage(
  dashboardHeader(title = "GenTile: CRISPR Guide Design"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Input & Parameters", tabName = "input", icon = icon("upload")),
      menuItem("Results", tabName = "results", icon = icon("table")),
      menuItem("Visualization", tabName = "plots", icon = icon("chart-line"))
    )
  ),
  
  dashboardBody(
    # Custom CSS for styling
    tags$head(
      tags$style(HTML("
        .progress-text { font-weight: bold; margin-top: 10px; }
        .parameter-box { background-color: #f8f9fa; padding: 15px; border-radius: 5px; margin: 10px 0; }
        .mode-selection { background-color: #e3f2fd; padding: 10px; border-radius: 5px; }
      "))
    ),
    
    tabItems(
      # Input Tab
      tabItem(tabName = "input",
        fluidRow(
          # Left column - Input
          box(
            title = "Input Data", status = "primary", solidHeader = TRUE, width = 6,
            
            # Input method selection
            radioButtons("input_method", "Input Type (auto-detected):",
                        choices = c("Paste gene list" = "paste_genes",
                                   "Upload gene file" = "upload_genes",
                                   "Paste coordinates" = "paste_coords",
                                   "Upload coordinates file" = "upload_coords"),
                        selected = "paste_genes"),
            
            # Conditional input panels
            conditionalPanel(
              condition = "input.input_method == 'paste_genes'",
              textAreaInput("gene_text", "Gene List (one per line):",
                           placeholder = "TP53\nBRCA1\nEGFR\n...",
                           rows = 8)
            ),
            
            conditionalPanel(
              condition = "input.input_method == 'upload_genes'",
              fileInput("gene_file", "Upload Gene File",
                       accept = c(".txt", ".csv"))
            ),
            
            conditionalPanel(
              condition = "input.input_method == 'paste_coords'",
              textAreaInput("coord_text", "Coordinates (name,chr:pos,strand):",
                           placeholder = "tp53_region,chr17:7687546,-\nbrca1_region,chr17:43170245,-\n...",
                           rows = 8)
            ),
            
            conditionalPanel(
              condition = "input.input_method == 'upload_coords'",
              fileInput("coord_file", "Upload Coordinates File",
                       accept = c(".txt", ".csv"))
            ),
            
            br(),
            div(class = "parameter-box",
                h4("Detected Input Type:"),
                verbatimTextOutput("input_detection")
            )
          ),
          
          # Right column - Parameters
          box(
            title = "Parameters", status = "info", solidHeader = TRUE, width = 6,
            
            # CAGE data options (only for gene input)
            conditionalPanel(
              condition = "input.input_method.includes('genes')",
              h4("TSS Source"),
              selectInput("cell_line", "Cell Line for CAGE data:",
                         choices = c("K562" = "K562", 
                                    "HeLa" = "HeLa",
                                    "HepG2" = "HepG2", 
                                    "GM12878" = "GM12878",
                                    "H1" = "H1",
                                    "IMR90" = "IMR90",
                                    "HUVEC" = "HUVEC",
                                    "keratinocyte" = "keratinocyte",
                                    "CD14_monocyte" = "CD14_monocyte",
                                    "muscle_satellite" = "muscle_satellite",
                                    "osteoblast" = "osteoblast",
                                    "dermal_fibroblast" = "dermal_fibroblast",
                                    "SK-N-SH" = "SK-N-SH"),
                         selected = "K562"),
              br()
            ),
            
            # Sequence parameters
            h4("Sequence Parameters"),
            sliderInput("upstream", "Upstream distance (bp):",
                       min = 500, max = 5000, value = 1500, step = 100),
            sliderInput("downstream", "Downstream distance (bp):", 
                       min = 200, max = 2000, value = 500, step = 50),
            br(),
            
            # Guide selection modes
            div(class = "mode-selection",
                h4("Guide Selection Modes"),
                checkboxGroupInput("selection_modes", NULL,
                                  choices = c("Tiling (spaced coverage)" = "tiling",
                                             "CRISPRi (knockdown)" = "crispri", 
                                             "CRISPRa (activation)" = "crispra"),
                                  selected = "tiling"),
                
                # Mode-specific parameters
                conditionalPanel(
                  condition = "input.selection_modes.includes('tiling')",
                  sliderInput("zone_size", "Tiling exclusion zone (bp):",
                             min = 25, max = 200, value = 50, step = 25)
                ),
                
                conditionalPanel(
                  condition = "input.selection_modes.includes('crispri') || input.selection_modes.includes('crispra')",
                  sliderInput("target_guides", "Target guides per gene (CRISPRi/a):",
                             min = 1, max = 20, value = 5, step = 1)
                )
            ),
            br(),
            
            # Filtering options
            h4("Filtering Options"),
            checkboxInput("use_k562_variants", "Filter K562 variants", FALSE),
            checkboxInput("filter_restriction", "Filter restriction sites", FALSE),
            conditionalPanel(
              condition = "input.filter_restriction",
              selectizeInput("restriction_enzymes", "Enzymes to avoid:",
                            choices = NULL,
                            multiple = TRUE,
                            options = list(placeholder = "Type to search (e.g., BsaI)..."))
            )
          )
        ),
        
        # Run button and progress
        fluidRow(
          box(
            title = "Run Analysis", status = "success", solidHeader = TRUE, width = 12,
            
            div(style = "text-align: center;",
                actionButton("run_analysis", "Run GenTile Pipeline", 
                           class = "btn-primary btn-lg",
                           style = "margin: 20px;"),
                br(),
                div(class = "progress-text",
                    textOutput("progress_text")),
                br(),
                progressBar("progress_bar", value = 0, display_pct = TRUE, status = "primary")
            )
          )
        )
      ),
      
      # Results Tab
      tabItem(tabName = "results",
        fluidRow(
          box(
            title = "Guide Selection Results", status = "primary", solidHeader = TRUE, width = 12,
            
            # Summary statistics
            fluidRow(
              column(2, valueBoxOutput("total_guides", width = NULL)),
              column(2, valueBoxOutput("tiling_guides", width = NULL)), 
              column(2, valueBoxOutput("crispri_guides", width = NULL)),
              column(2, valueBoxOutput("crispra_guides", width = NULL)),
              column(2, valueBoxOutput("avg_guides_per_target", width = NULL)),
              column(2, valueBoxOutput("min_guides_per_target", width = NULL))
            ),
            
            # Results table with selection controls
            div(
              style = "margin-bottom: 10px;",
              actionButton("select_all", "Select All", class = "btn-sm btn-outline-primary"),
              actionButton("select_none", "Select None", class = "btn-sm btn-outline-secondary"),
              span(style = "margin-left: 20px; font-style: italic;", 
                   textOutput("selection_info", inline = TRUE))
            ),
            DTOutput("results_table"),
            
            # Download buttons
            br(),
            downloadButton("download_guides", "Download Selected Guides (.txt)", class = "btn-info"),
            downloadButton("download_bed", "Download BED file (.bed)", class = "btn-info")
          )
        )
      ),
      
      # Visualization Tab
      tabItem(tabName = "plots",
        fluidRow(
          box(
            title = "Guide Visualization", status = "primary", solidHeader = TRUE, width = 12,
            
            # Gene selector for visualization
            conditionalPanel(
              condition = "output.results_available",
              selectInput("selected_gene", "Select gene for visualization:",
                         choices = NULL),
              br(),
              
              # Genome track plot
              plotlyOutput("genome_plot", height = "400px"),
              br(),
              
            ),
            
            conditionalPanel(
              condition = "!output.results_available",
              h3("No results available", style = "text-align: center; color: #888;"),
              p("Run the analysis first to see visualizations.", style = "text-align: center;")
            )
          )
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Reactive values to store data
  values <- reactiveValues(
    pipeline_results = NULL,
    guides_data = NULL,
    input_genes = NULL,
    current_status = "Ready"
  )
  
  # Input detection
  output$input_detection <- renderText({
    input_text <- ""
    
    if (input$input_method == "paste_genes") {
      input_text <- input$gene_text
    } else if (input$input_method == "paste_coords") {
      input_text <- input$coord_text
    }
    
    if (nchar(input_text) == 0) {
      return("No input provided")
    }
    
    # Simple detection logic
    lines <- trimws(strsplit(input_text, "\n")[[1]])
    lines <- lines[nchar(lines) > 0]
    
    if (length(lines) == 0) return("No valid input")
    
    # Check if it looks like coordinates (contains commas and colons)
    coord_pattern <- any(grepl("^[^,]+,chr[^:]+:[0-9]+,[+-]$", lines))
    
    if (coord_pattern) {
      return(paste("Coordinates format detected -", length(lines), "entries"))
    } else {
      return(paste("Gene list detected -", length(lines), "genes"))
    }
  })
  
  # Update restriction enzyme choices
  observe({
    if (input$filter_restriction) {
      # Load enzyme database
      enzyme_db_path <- file.path(GENTILE_PATH, "data/reference/restriction_enzymes/enzymes.tsv")
      if (file.exists(enzyme_db_path)) {
        enzymes <- read.delim(enzyme_db_path, stringsAsFactors = FALSE)
        enzyme_choices <- setNames(enzymes$enzyme_name, enzymes$enzyme_name)
        updateSelectizeInput(session, "restriction_enzymes", choices = enzyme_choices)
      }
    }
  })
  
  # Progress tracking
  output$progress_text <- renderText({
    values$current_status
  })
  
  # Main analysis runner
  observeEvent(input$run_analysis, {
    
    # Validate inputs
    if (input$input_method %in% c("paste_genes", "paste_coords")) {
      input_text <- if (input$input_method == "paste_genes") input$gene_text else input$coord_text
      if (nchar(trimws(input_text)) == 0) {
        showNotification("Please provide input data", type = "error", duration = 5)
        return()
      }
    }
    
    if (length(input$selection_modes) == 0) {
      showNotification("Please select at least one guide selection mode", type = "error", duration = 5)
      return()
    }
    
    # Update progress
    updateProgressBar(session, "progress_bar", value = 0)
    values$current_status <- "Starting pipeline..."
    
    # Run pipeline (this will be implemented in pipeline_runner.R)
    tryCatch({
      results <- run_gentile_pipeline(
        input_method = input$input_method,
        input_data = list(
          gene_text = input$gene_text,
          coord_text = input$coord_text,
          gene_file = input$gene_file,
          coord_file = input$coord_file
        ),
        parameters = list(
          cell_line = input$cell_line,
          upstream = input$upstream,
          downstream = input$downstream,
          selection_modes = input$selection_modes,
          zone_size = input$zone_size,
          target_guides = input$target_guides,
          use_k562_variants = input$use_k562_variants,
          filter_restriction = input$filter_restriction,
          restriction_enzymes = input$restriction_enzymes
        ),
        progress_callback = function(message, progress) {
          values$current_status <- message
          updateProgressBar(session, "progress_bar", value = progress)
        }
      )
      
      values$pipeline_results <- results
      values$guides_data <- results$guides
      values$current_status <- "Pipeline completed successfully!"
      updateProgressBar(session, "progress_bar", value = 100)
      
      showNotification("Analysis completed successfully!", type = "message", duration = 3)
      
    }, error = function(e) {
      values$current_status <- paste("Error:", e$message)
      showNotification(paste("Pipeline failed:", e$message), type = "error", duration = 10)
    })
  })
  
  # Results available flag
  output$results_available <- reactive({
    !is.null(values$guides_data)
  })
  outputOptions(output, "results_available", suspendWhenHidden = FALSE)
  
  # Summary statistics
  output$total_guides <- renderValueBox({
    count <- if (!is.null(values$guides_data)) nrow(values$guides_data) else 0
    valueBox(
      value = count,
      subtitle = "Total Guides",
      icon = icon("dna"),
      color = "blue"
    )
  })
  
  output$tiling_guides <- renderValueBox({
    count <- if (!is.null(values$guides_data)) {
      sum(values$guides_data$tiling_guide == "TRUE", na.rm = TRUE)
    } else 0
    valueBox(
      value = count,
      subtitle = "Tiling Guides", 
      icon = icon("grip-lines"),
      color = "green"
    )
  })
  
  output$crispri_guides <- renderValueBox({
    count <- if (!is.null(values$guides_data)) {
      sum(values$guides_data$crispri_guide == "TRUE", na.rm = TRUE)
    } else 0
    valueBox(
      value = count,
      subtitle = "CRISPRi Guides",
      icon = icon("arrow-down"),
      color = "orange"
    )
  })
  
  output$crispra_guides <- renderValueBox({
    count <- if (!is.null(values$guides_data)) {
      sum(values$guides_data$crispra_guide == "TRUE", na.rm = TRUE)
    } else 0
    valueBox(
      value = count,
      subtitle = "CRISPRa Guides",
      icon = icon("arrow-up"),
      color = "purple"
    )
  })
  
  # Additional summary stats
  output$avg_guides_per_target <- renderValueBox({
    avg <- if (!is.null(values$guides_data)) {
      targets <- unique(values$guides_data$gene)
      round(nrow(values$guides_data) / length(targets), 1)
    } else 0
    valueBox(
      value = avg,
      subtitle = "Avg Guides/Target",
      icon = icon("calculator"),
      color = "navy"
    )
  })
  
  output$min_guides_per_target <- renderValueBox({
    min_val <- if (!is.null(values$guides_data)) {
      guide_counts <- table(values$guides_data$gene)
      min(guide_counts)
    } else 0
    valueBox(
      value = min_val,
      subtitle = "Min Guides/Target",
      icon = icon("arrow-down"),
      color = "maroon"
    )
  })
  
  # Results table
  output$results_table <- renderDT({
    if (is.null(values$guides_data)) return(NULL)
    
    # Display key columns for the table
    display_data <- values$guides_data[, c("guide_id", "gene", "Hsu2013", 
                                          "tiling_guide", "crispri_guide", "crispra_guide")]
    
    datatable(display_data, 
              options = list(pageLength = 25, scrollX = TRUE),
              filter = 'top',
              selection = 'multiple')
  })
  
  # Selection controls for data table
  observeEvent(input$select_all, {
    if (!is.null(values$guides_data)) {
      dataTableProxy("results_table") %>% 
        selectRows(seq_len(nrow(values$guides_data)))
    }
  })
  
  observeEvent(input$select_none, {
    dataTableProxy("results_table") %>% selectRows(NULL)
  })
  
  # Selection info
  output$selection_info <- renderText({
    if (is.null(values$guides_data)) return("")
    
    selected_rows <- input$results_table_rows_selected
    selected_count <- if(is.null(selected_rows)) 0 else length(selected_rows)
    total_count <- nrow(values$guides_data)
    
    if (selected_count == 0) {
      "No rows selected (downloads will include all data)"
    } else {
      paste(selected_count, "of", total_count, "guides selected")
    }
  })
  output$download_guides <- downloadHandler(
    filename = function() {
      paste0("genTile_guides_", Sys.Date(), ".txt")
    },
    content = function(file) {
      if (is.null(values$guides_data)) {
        write("No data available", file)
        return()
      }
      
      # Get selected rows or all if none selected
      selected_rows <- input$results_table_rows_selected
      if (is.null(selected_rows) || length(selected_rows) == 0) {
        data_to_download <- values$guides_data
      } else {
        data_to_download <- values$guides_data[selected_rows, ]
      }
      
      write.table(data_to_download, file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )
  
  output$download_bed <- downloadHandler(
    filename = function() {
      paste0("genTile_guides_", Sys.Date(), ".bed")
    },
    content = function(file) {
      if (is.null(values$pipeline_results$bed_data)) {
        write("No BED data available", file)
        return()
      }
      
      # Get selected rows or all if none selected
      selected_rows <- input$results_table_rows_selected
      if (is.null(selected_rows) || length(selected_rows) == 0) {
        bed_to_download <- values$pipeline_results$bed_data
      } else {
        # Filter BED data based on selected guide IDs
        selected_guides <- values$guides_data$guide_id[selected_rows]
        bed_to_download <- values$pipeline_results$bed_data[
          values$pipeline_results$bed_data$guide_id %in% selected_guides, 
        ]
      }
      
      write.table(bed_to_download, file, sep = "\t", row.names = FALSE, 
                 col.names = FALSE, quote = FALSE)
    }
  )
  
  # Gene selector for visualization (fix the naming issue)
  observe({
    if (!is.null(values$guides_data)) {
      # Extract proper gene names from the data
      genes <- unique(values$guides_data$gene)
      # Remove any sequence-like entries and keep only gene names
      genes <- genes[!grepl("^[ATCG]+$", genes)]  # Remove pure DNA sequences
      genes <- genes[nchar(genes) < 50]  # Remove very long entries
      
      updateSelectInput(session, "selected_gene", choices = genes, selected = genes[1])
    }
  })
  
  # Visualization plots
  output$genome_plot <- renderPlotly({
    if (is.null(values$guides_data) || is.null(input$selected_gene)) {
      return(plot_ly() %>% add_text(text = "No data available for visualization", 
                                   x = 0.5, y = 0.5, showlegend = FALSE))
    }
    
    create_genome_plot(values$guides_data, input$selected_gene, values$pipeline_results$bed_data)
  })
  
}

# Run the app
shinyApp(ui = ui, server = server)
