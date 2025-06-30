# Plotting Functions for GenTile Shiny App
# functions/plotting_functions.R

library(plotly)
library(dplyr)
library(rlang)  # For %||% operator

# Create genome track visualization for a single gene
create_genome_plot <- function(guides_data, selected_gene) {
  
  # Filter guides for selected gene
  gene_guides <- guides_data[guides_data$gene == selected_gene, ]
  
  if (nrow(gene_guides) == 0) {
    return(plot_ly() %>% 
           add_text(text = paste("No guides found for", selected_gene), 
                   x = 0.5, y = 0.5, showlegend = FALSE))
  }
  
  # Extract genomic coordinates and TSS from contig
  # Format: GENE::chr:start-end(strand)::TSS_chr:position
  contig_info <- parse_contig_info(gene_guides$contig[1])
  
  if (is.null(contig_info)) {
    return(plot_ly() %>% 
           add_text(text = "Could not parse genomic coordinates", 
                   x = 0.5, y = 0.5, showlegend = FALSE))
  }
  
  # Calculate absolute coordinates for guides
  gene_guides$abs_start <- contig_info$region_start + gene_guides$start
  gene_guides$abs_end <- contig_info$region_start + gene_guides$stop
  
  # Create base plot
  p <- plot_ly()
  
  # Add gene region track
  p <- p %>% add_trace(
    x = c(contig_info$region_start, contig_info$region_end),
    y = c(4, 4),
    type = "scatter",
    mode = "lines",
    line = list(color = "#333333", width = 8),
    name = paste("Gene", selected_gene),
    showlegend = TRUE
  )
  
  # Add TSS marker
  p <- p %>% add_trace(
    x = c(contig_info$tss_pos, contig_info$tss_pos),
    y = c(3.7, 4.3),
    type = "scatter",
    mode = "lines",
    line = list(color = "#FF0000", width = 3),
    name = "TSS",
    showlegend = TRUE
  )
  
  # Add CRISPRa region highlight (-400 to -50 from TSS)
  if (any(gene_guides$crispra_guide == "TRUE")) {
    crispra_start <- contig_info$tss_pos - 400
    crispra_end <- contig_info$tss_pos - 50
    
    p <- p %>% add_trace(
      x = c(crispra_start, crispra_end, crispra_end, crispra_start),
      y = c(1.8, 1.8, 2.2, 2.2),
      type = "scatter",
      mode = "lines",
      fill = "toself",
      fillcolor = "rgba(153, 50, 204, 0.2)",
      line = list(color = "rgba(153, 50, 204, 0.4)", width = 1),
      name = "CRISPRa region (-400 to -50bp)",
      showlegend = TRUE
    )
  }
  
  # Add CRISPRi region highlight (-50 to +300 from TSS)
  if (any(gene_guides$crispri_guide == "TRUE")) {
    crispri_start <- contig_info$tss_pos - 50
    crispri_end <- contig_info$tss_pos + 300
    
    p <- p %>% add_trace(
      x = c(crispri_start, crispri_end, crispri_end, crispri_start),
      y = c(1.8, 1.8, 2.2, 2.2),
      type = "scatter",
      mode = "lines",
      fill = "toself",
      fillcolor = "rgba(220, 20, 60, 0.2)",
      line = list(color = "rgba(220, 20, 60, 0.4)", width = 1),
      name = "CRISPRi region (-50 to +300bp)",
      showlegend = TRUE
    )
  }
  
  # Add guide tracks by mode with 5-color system
  mode_data <- classify_guide_modes(gene_guides)
  
  for (mode in names(mode_data)) {
    mode_guides <- mode_data[[mode]]
    
    if (nrow(mode_guides) > 0) {
      p <- p %>% add_trace(
        x = mode_guides$abs_start,
        y = rep(get_mode_y_position(mode), nrow(mode_guides)),
        type = "scatter",
        mode = "markers",
        marker = list(
          color = get_mode_color(mode),
          size = 12,
          symbol = get_mode_symbol(mode),
          line = list(color = "#000000", width = 1)
        ),
        text = paste(
          "Guide:", mode_guides$target,
          "<br>Gene:", mode_guides$gene,
          "<br>Score:", round(mode_guides$Hsu2013, 1),
          "<br>Mode:", mode,
          "<br>Position:", mode_guides$abs_start
        ),
        hoverinfo = "text",
        name = mode,
        showlegend = TRUE
      )
    }
  }
  
  # Format plot
  p <- p %>% layout(
    title = paste("CRISPR Guide Design for", selected_gene),
    xaxis = list(
      title = "Genomic Position",
      tickformat = ",d"
    ),
    yaxis = list(
      title = "",
      range = c(0.5, 5),
      tickvals = c(1, 2, 3, 4),
      ticktext = c("Guides", "Target Regions", "TSS", "Gene"),
      showgrid = FALSE
    ),
    hovermode = "closest",
    showlegend = TRUE,
    legend = list(x = 1.02, y = 1)
  )
  
  return(p)
}

# Parse contig information to extract coordinates and TSS
parse_contig_info <- function(contig_string) {
  # Format: GENE::chr:start-end(strand)::TSS_chr:position
  parts <- strsplit(contig_string, "::")[[1]]
  
  if (length(parts) < 3) return(NULL)
  
  # Parse region: chr:start-end(strand)
  region_match <- regexpr("([^:]+):([0-9]+)-([0-9]+)\\(([+-])\\)", parts[2], perl = TRUE)
  if (region_match == -1) return(NULL)
  
  region_text <- regmatches(parts[2], region_match)
  region_parts <- regmatches(region_text, gregexpr("[^:()+-]+", region_text, perl = TRUE))[[1]]
  
  chr <- region_parts[1]
  region_start <- as.numeric(region_parts[2])
  region_end <- as.numeric(region_parts[3])
  strand <- ifelse(grepl("\\+", region_text), "+", "-")
  
  # Parse TSS: TSS_chr:position
  tss_match <- regexpr("TSS_[^:]+:([0-9]+)", parts[3], perl = TRUE)
  if (tss_match == -1) return(NULL)
  
  tss_text <- regmatches(parts[3], tss_match)
  tss_pos <- as.numeric(gsub(".*:([0-9]+)", "\\1", tss_text))
  
  return(list(
    chr = chr,
    region_start = region_start,
    region_end = region_end,
    strand = strand,
    tss_pos = tss_pos
  ))
}

# Classify guides by mode for 5-color system
classify_guide_modes <- function(gene_guides) {
  
  mode_data <- list()
  
  # Tiling only (Sea Green)
  mode_data[["Tiling only"]] <- gene_guides[
    gene_guides$tiling_guide == "TRUE" & 
    gene_guides$crispri_guide == "FALSE" & 
    gene_guides$crispra_guide == "FALSE", 
  ]
  
  # CRISPRi only (Crimson)
  mode_data[["CRISPRi only"]] <- gene_guides[
    gene_guides$tiling_guide == "FALSE" & 
    gene_guides$crispri_guide == "TRUE" & 
    gene_guides$crispra_guide == "FALSE", 
  ]
  
  # CRISPRa only (Dark Orchid)
  mode_data[["CRISPRa only"]] <- gene_guides[
    gene_guides$tiling_guide == "FALSE" & 
    gene_guides$crispri_guide == "FALSE" & 
    gene_guides$crispra_guide == "TRUE", 
  ]
  
  # Tiling + CRISPRi (Orange)
  mode_data[["Tiling + CRISPRi"]] <- gene_guides[
    gene_guides$tiling_guide == "TRUE" & 
    gene_guides$crispri_guide == "TRUE", 
  ]
  
  # Tiling + CRISPRa (Medium Slate Blue)
  mode_data[["Tiling + CRISPRa"]] <- gene_guides[
    gene_guides$tiling_guide == "TRUE" & 
    gene_guides$crispra_guide == "TRUE", 
  ]
  
  return(mode_data)
}

# Get color for each mode
get_mode_color <- function(mode) {
  colors <- list(
    "Tiling only" = "#2E8B57",      # Sea Green
    "CRISPRi only" = "#DC143C",     # Crimson
    "CRISPRa only" = "#9932CC",     # Dark Orchid
    "Tiling + CRISPRi" = "#FF6B35", # Orange
    "Tiling + CRISPRa" = "#7B68EE"  # Medium Slate Blue
  )
  return(colors[[mode]] %||% "#708090")
}

# Get symbol for each mode
get_mode_symbol <- function(mode) {
  symbols <- list(
    "Tiling only" = "circle",
    "CRISPRi only" = "triangle-down",
    "CRISPRa only" = "triangle-up", 
    "Tiling + CRISPRi" = "diamond",
    "Tiling + CRISPRa" = "square"
  )
  return(symbols[[mode]] %||% "circle")
}

# Get Y position for each mode (stack them)
get_mode_y_position <- function(mode) {
  positions <- list(
    "Tiling only" = 1.5,
    "CRISPRi only" = 1.3,
    "CRISPRa only" = 1.1,
    "Tiling + CRISPRi" = 0.9,
    "Tiling + CRISPRa" = 0.7
  )
  return(positions[[mode]] %||% 1.0)
}$crispra_guide == "TRUE")) {
    crispra_start <- tss_info$tss_pos - 400
    crispra_end <- tss_info$tss_pos - 50
    
    p <- p %>% add_trace(
      x = c(crispra_start, crispra_end, crispra_end, crispra_start),
      y = c(1.8, 1.8, 2.2, 2.2),
      type = "scatter",
      mode = "lines",
      fill = "toself",
      fillcolor = "rgba(153, 50, 204, 0.2)",
      line = list(color = "rgba(153, 50, 204, 0.4)", width = 1),
      name = "CRISPRa region",
      showlegend = TRUE
    )
  }
  
  # Add CRISPRi region highlight (-50 to +300 from TSS)
  if (!is.null(tss_info) && any(gene_guides$crispri_guide == "TRUE")) {
    crispri_start <- tss_info$tss_pos - 50
    crispri_end <- tss_info$tss_pos + 300
    
    p <- p %>% add_trace(
      x = c(crispri_start, crispri_end, crispri_end, crispri_start),
      y = c(1.8, 1.8, 2.2, 2.2),
      type = "scatter",
      mode = "lines",
      fill = "toself",
      fillcolor = "rgba(220, 20, 60, 0.2)",
      line = list(color = "rgba(220, 20, 60, 0.4)", width = 1),
      name = "CRISPRi region",
      showlegend = TRUE
    )
  }
  
  # Add guide tracks by mode
  mode_levels <- c("Tiling only", "CRISPRi only", "CRISPRa only", 
                   "Tiling + CRISPRi", "Tiling + CRISPRa")
  
  for (i in 1:length(mode_levels)) {
    mode <- mode_levels[i]
    mode_guides <- filter_guides_by_mode(gene_guides, mode)
    
    if (nrow(mode_guides) > 0) {
      y_pos <- 1.5 - (i * 0.2)  # Stack different modes
      
      p <- p %>% add_trace(
        x = mode_guides$abs_start,
        y = rep(y_pos, nrow(mode_guides)),
        type = "scatter",
        mode = "markers",
        marker = list(
          color = mode_guides$guide_color,
          size = 10,
          symbol = get_guide_symbol(mode),
          line = list(color = "#000000", width = 1)
        ),
        text = paste(
          "Guide:", mode_guides$guide_id,
          "<br>Score:", round(mode_guides$Hsu2013, 1),
          "<br>Mode:", mode,
          "<br>Position:", mode_guides$abs_start
        ),
        hoverinfo = "text",
        name = mode,
        showlegend = TRUE
      )
    }
  }
  
  # Format plot
  p <- p %>% layout(
    title = paste("CRISPR Guide Design for", selected_gene),
    xaxis = list(
      title = "Genomic Position",
      tickformat = ",d"
    ),
    yaxis = list(
      title = "",
      range = c(0.5, 5),
      tickvals = c(1, 2, 3, 4),
      ticktext = c("Guides", "Target Regions", "CAGE/TSS", "Gene"),
      showgrid = FALSE
    ),
    hovermode = "closest",
    showlegend = TRUE,
    legend = list(x = 1.02, y = 1)
  )
  
  return(p)
}

# Helper function to extract TSS information from guides
extract_tss_from_guides <- function(gene_guides) {
  
  # Try to extract TSS from contig header
  if ("contig" %in% names(gene_guides) && nrow(gene_guides) > 0) {
    contig_parts <- strsplit(gene_guides$contig[1], "::")[[1]]
    
    if (length(contig_parts) >= 3) {
      tss_part <- contig_parts[3]
      if (grepl("TSS_", tss_part)) {
        tss_match <- regexpr("TSS_[^:]+:([0-9]+)", tss_part, perl = TRUE)
        if (tss_match > 0) {
          tss_pos <- as.numeric(regmatches(tss_part, tss_match))
          tss_pos <- as.numeric(gsub(".*:([0-9]+)", "\\1", tss_pos))
          
          return(list(tss_pos = tss_pos))
        }
      }
    }
  }
  
  # Fallback: use middle of guides region
  if ("abs_start" %in% names(gene_guides)) {
    mid_pos <- mean(c(min(gene_guides$abs_start), max(gene_guides$abs_end)), na.rm = TRUE)
    return(list(tss_pos = mid_pos))
  }
  
  return(NULL)
}

# Helper function to filter guides by mode
filter_guides_by_mode <- function(guides_data, mode) {
  
  switch(mode,
    "Tiling only" = guides_data[
      guides_data$tiling_guide == "TRUE" & 
      guides_data$crispri_guide == "FALSE" & 
      guides_data$crispra_guide == "FALSE", 
    ],
    "CRISPRi only" = guides_data[
      guides_data$tiling_guide == "FALSE" & 
      guides_data$crispri_guide == "TRUE" & 
      guides_data$crispra_guide == "FALSE", 
    ],
    "CRISPRa only" = guides_data[
      guides_data$tiling_guide == "FALSE" & 
      guides_data$crispri_guide == "FALSE" & 
      guides_data$crispra_guide == "TRUE", 
    ],
    "Tiling + CRISPRi" = guides_data[
      guides_data$tiling_guide == "TRUE" & 
      guides_data$crispri_guide == "TRUE", 
    ],
    "Tiling + CRISPRa" = guides_data[
      guides_data$tiling_guide == "TRUE" & 
      guides_data$crispra_guide == "TRUE", 
    ],
    data.frame()  # Default empty
  )
}

# Helper function to get marker symbols for different modes
get_guide_symbol <- function(mode) {
  switch(mode,
    "Tiling only" = "circle",
    "CRISPRi only" = "triangle-down",
    "CRISPRa only" = "triangle-up", 
    "Tiling + CRISPRi" = "diamond",
    "Tiling + CRISPRa" = "square",
    "circle"  # Default
  )
}

# Create score distribution plot
create_score_distribution_plot <- function(guides_data) {
  
  if (!"Hsu2013" %in% names(guides_data) || nrow(guides_data) == 0) {
    return(plot_ly() %>% 
           add_text(text = "No score data available", 
                   x = 0.5, y = 0.5, showlegend = FALSE))
  }
  
  # Prepare data with mode labels
  plot_data <- guides_data %>%
    mutate(
      mode = case_when(
        tiling_guide == "TRUE" & crispri_guide == "TRUE" ~ "Tiling + CRISPRi",
        tiling_guide == "TRUE" & crispra_guide == "TRUE" ~ "Tiling + CRISPRa", 
        tiling_guide == "TRUE" ~ "Tiling only",
        crispri_guide == "TRUE" ~ "CRISPRi only",
        crispra_guide == "TRUE" ~ "CRISPRa only",
        TRUE ~ "Other"
      )
    )
  
  # Create histogram
  p <- plot_ly(
    data = plot_data,
    x = ~Hsu2013,
    color = ~mode,
    colors = c(
      "Tiling only" = "#2E8B57",
      "CRISPRi only" = "#DC143C", 
      "CRISPRa only" = "#9932CC",
      "Tiling + CRISPRi" = "#FF6B35",
      "Tiling + CRISPRa" = "#7B68EE"
    ),
    type = "histogram",
    alpha = 0.7,
    nbinsx = 20
  ) %>%
  layout(
    title = "Guide Score Distribution by Mode",
    xaxis = list(title = "Hsu2013 Score"),
    yaxis = list(title = "Number of Guides"),
    barmode = "overlay",
    showlegend = TRUE
  )
  
  return(p)
}

# Create guides per gene summary plot
create_guides_per_gene_plot <- function(guides_data) {
  
  if (nrow(guides_data) == 0) {
    return(plot_ly() %>% 
           add_text(text = "No guide data available", 
                   x = 0.5, y = 0.5, showlegend = FALSE))
  }
  
  # Count guides per gene by mode
  summary_data <- guides_data %>%
    group_by(gene) %>%
    summarise(
      Tiling = sum(tiling_guide == "TRUE", na.rm = TRUE),
      CRISPRi = sum(crispri_guide == "TRUE", na.rm = TRUE),
      CRISPRa = sum(crispra_guide == "TRUE", na.rm = TRUE),
      Total = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(Total))
  
  # Create stacked bar plot
  p <- plot_ly(summary_data, x = ~gene, y = ~Tiling, type = "bar", name = "Tiling", 
               marker = list(color = "#2E8B57")) %>%
    add_trace(y = ~CRISPRi, name = "CRISPRi", marker = list(color = "#DC143C")) %>%
    add_trace(y = ~CRISPRa, name = "CRISPRa", marker = list(color = "#9932CC")) %>%
    layout(
      title = "Guide Count by Gene and Mode",
      xaxis = list(title = "Gene", tickangle = -45),
      yaxis = list(title = "Number of Guides"),
      barmode = "stack",
      showlegend = TRUE
    )
  
  return(p)
}
