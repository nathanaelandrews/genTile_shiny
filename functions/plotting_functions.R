# Plotting Functions for GenTile Shiny App
# functions/plotting_functions.R

library(plotly)
library(dplyr)
library(rlang)  # For %||% operator

# =============================================================================
# MAIN PLOTTING FUNCTIONS
# =============================================================================

# Create genome track visualization for a single gene
create_genome_plot <- function(guides_data, selected_gene, bed_data = NULL) {
  
  # Filter guides for selected gene
  gene_guides <- guides_data[guides_data$gene == selected_gene, ]
  
  if (nrow(gene_guides) == 0) {
    return(create_empty_plot(paste("No guides found for", selected_gene)))
  }
  
  # Parse genomic coordinates and TSS from contig
  contig_info <- parse_contig_info(gene_guides$contig[1])
  
  if (is.null(contig_info)) {
    return(create_empty_plot("Could not parse genomic coordinates"))
  }
  
  # Calculate absolute coordinates for guides
  gene_guides <- add_absolute_coordinates(gene_guides, contig_info)
  
  # Set plot range: TSS +/- 2000bp for better visualization
  plot_start <- contig_info$tss_pos - 2000
  plot_end <- contig_info$tss_pos + 2000
  
  # Create the plot step by step
  p <- plot_ly()
  
  # Add gene region track
  p <- add_gene_track(p, contig_info, selected_gene, plot_start, plot_end)
  
  # Add TSS marker
  p <- add_tss_marker(p, contig_info)
  
  # Add functional region highlights
  p <- add_functional_regions(p, gene_guides, contig_info)
  
  # Add guide markers as rectangles
  p <- add_guide_rectangles(p, gene_guides)
  
  # Format the final plot
  p <- format_genome_plot(p, selected_gene, plot_start, plot_end)
  
  return(p)
}

# This function has been removed as requested

# Create guides per gene summary plot
create_guides_per_gene_plot <- function(guides_data) {
  
  if (nrow(guides_data) == 0) {
    return(create_empty_plot("No guide data available"))
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

# =============================================================================
# COORDINATE PARSING AND DATA PREPARATION
# =============================================================================

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

# Add absolute genomic coordinates to guides data
add_absolute_coordinates <- function(gene_guides, contig_info) {
  gene_guides$abs_start <- contig_info$region_start + gene_guides$start
  gene_guides$abs_end <- contig_info$region_start + gene_guides$stop
  return(gene_guides)
}

# =============================================================================
# PLOT BUILDING FUNCTIONS
# =============================================================================

# Create empty plot with message
create_empty_plot <- function(message) {
  plot_ly() %>% 
    add_text(text = message, x = 0.5, y = 0.5, showlegend = FALSE) %>%
    layout(
      xaxis = list(showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE),
      yaxis = list(showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
    )
}

# Add gene region track to plot
add_gene_track <- function(p, contig_info, gene_name, plot_start, plot_end) {
  # Show gene extent within the plot range
  gene_start <- max(contig_info$region_start, plot_start)
  gene_end <- min(contig_info$region_end, plot_end)
  
  p %>% add_trace(
    x = c(gene_start, gene_end),
    y = c(3.5, 3.5),
    type = "scatter",
    mode = "lines",
    line = list(color = "#333333", width = 8),
    name = paste("Gene", gene_name),
    showlegend = TRUE
  )
}

# Add TSS marker to plot
# Add TSS marker to plot
add_tss_marker <- function(p, contig_info) {
  p %>% add_trace(
    x = c(contig_info$tss_pos, contig_info$tss_pos),
    y = c(3.2, 3.8),
    type = "scatter",
    mode = "lines",
    line = list(color = "#FF0000", width = 3),
    name = "TSS",
    showlegend = TRUE
  )
}

# Add variants track at the top of the plot
add_variants_track <- function(p, contig_info, plot_start, plot_end) {
  
  # Path to K562 variants file
  variants_file <- "/Users/nathanaelandrews/wrk/github/genTile/data/reference/vcf/K562_variants_hg38_sorted.bed"
  
  # Check if variants file exists
  if (!file.exists(variants_file)) {
    cat("Variants file not found: ", variants_file, "\n")
    return(p)
  }
  
  tryCatch({
    # Read variants in the region of interest
    variants <- read.delim(variants_file, header = FALSE, stringsAsFactors = FALSE)
    colnames(variants)[1:3] <- c("chr", "start", "end")
    
    # Debug: print info to console
    cat("Plot chromosome: ", contig_info$chr, "\n")
    cat("Plot region: ", plot_start, "-", plot_end, "\n")
    cat("Total variants in file: ", nrow(variants), "\n")
    
    # Filter variants to the current chromosome and plot region
    chr_variants <- variants[
      variants$chr == contig_info$chr & 
      variants$start >= plot_start & 
      variants$end <= plot_end,
    ]
    
    cat("Variants in plot region: ", nrow(chr_variants), "\n")
    
    # If we have variants in this region, add them to the plot
    if (nrow(chr_variants) > 0) {
      # Limit to first 100 variants to avoid plot overload
      if (nrow(chr_variants) > 100) {
        chr_variants <- chr_variants[1:100, ]
        cat("Showing first 100 variants\n")
      }
      
      # Add each variant as a separate vertical line (more reliable)
      for (i in 1:nrow(chr_variants)) {
        var_pos <- chr_variants$start[i]
        
        p <- p %>% add_trace(
          x = c(var_pos, var_pos),
          y = c(5.0, 5.2),
          type = "scatter",
          mode = "lines",
          line = list(color = "#FFA500", width = 3),
          name = if (i == 1) paste("K562 Variants (", nrow(chr_variants), ")") else "",
          showlegend = if (i == 1) TRUE else FALSE,
          hoverinfo = "text",
          text = paste("Variant at", format(var_pos, big.mark = ","))
        )
      }
      
      cat("Added ", nrow(chr_variants), " variant markers to plot\n")
    } else {
      cat("No variants found in this region\n")
    }
    
  }, error = function(e) {
    cat("Error loading variants track: ", e$message, "\n")
  })
  
  return(p)
}

# Add functional region highlights (CRISPRa and CRISPRi zones)
add_functional_regions <- function(p, gene_guides, contig_info) {
  
  # Add CRISPRa region highlight (-400 to -50 from TSS) - moved above gene
  if (any(gene_guides$crispra_guide == "TRUE")) {
    crispra_start <- contig_info$tss_pos - 400
    crispra_end <- contig_info$tss_pos - 50
    
    p <- p %>% add_trace(
      x = c(crispra_start, crispra_end, crispra_end, crispra_start),
      y = c(4.3, 4.3, 4.7, 4.7),
      type = "scatter",
      mode = "lines",
      fill = "toself",
      fillcolor = "rgba(153, 50, 204, 0.2)",
      line = list(color = "rgba(153, 50, 204, 0.4)", width = 1, connectgaps = TRUE),
      name = "CRISPRa region (-400 to -50bp)",
      showlegend = TRUE
    )
  }
  
  # Add CRISPRi region highlight (-50 to +300 from TSS) - moved above gene
  if (any(gene_guides$crispri_guide == "TRUE")) {
    crispri_start <- contig_info$tss_pos - 50
    crispri_end <- contig_info$tss_pos + 300
    
    p <- p %>% add_trace(
      x = c(crispri_start, crispri_end, crispri_end, crispri_start),
      y = c(4.3, 4.3, 4.7, 4.7),
      type = "scatter",
      mode = "lines",
      fill = "toself",
      fillcolor = "rgba(220, 20, 60, 0.2)",
      line = list(color = "rgba(220, 20, 60, 0.4)", width = 1),
      name = "CRISPRi region (-50 to +300bp)",
      showlegend = TRUE
    )
  }
  
  return(p)
}

# Add guide rectangles to plot (20bp each)
add_guide_rectangles <- function(p, gene_guides) {
  
  # Classify guides by mode
  mode_data <- classify_guides_by_mode(gene_guides)
  
  # Rectangle width for guides (20bp)
  guide_width <- 20
  
  # Add each mode as separate trace with rectangles
  for (mode in names(mode_data)) {
    mode_guides <- mode_data[[mode]]
    
    if (nrow(mode_guides) > 0) {
      y_pos <- get_mode_y_position(mode)
      
      # Create all rectangles for this mode in a single trace
      all_x <- c()
      all_y <- c()
      all_text <- c()
      
      for (i in 1:nrow(mode_guides)) {
        guide_start <- mode_guides$abs_start[i]
        guide_end <- guide_start + guide_width
        
        # Add rectangle coordinates with NA to separate shapes
        rect_x <- c(guide_start, guide_end, guide_end, guide_start, guide_start, NA)
        rect_y <- c(y_pos - 0.05, y_pos - 0.05, y_pos + 0.05, y_pos + 0.05, y_pos - 0.05, NA)
        
        all_x <- c(all_x, rect_x)
        all_y <- c(all_y, rect_y)
        
        # Add hover text for each point in rectangle (will show for any part of rectangle)
        hover_text <- create_guide_hover_text(mode_guides[i,], mode)
        all_text <- c(all_text, rep(hover_text, 6))
      }
      
      p <- p %>% add_trace(
        x = all_x,
        y = all_y,
        type = "scatter",
        mode = "lines",
        fill = "toself",
        fillcolor = get_mode_color(mode),
        line = list(color = "#000000", width = 1),
        text = all_text,
        hoverinfo = "text",
        name = mode,
        showlegend = TRUE
      )
    }
  }
  
  return(p)
}

# Format the final genome plot
format_genome_plot <- function(p, gene_name, plot_start, plot_end) {
  p %>% layout(
    title = paste("CRISPR Guide Design for", gene_name),
    xaxis = list(
      title = "Genomic Position",
      tickformat = ",d",
      range = c(plot_start, plot_end)
    ),
    yaxis = list(
      title = "",
      range = c(0.8, 5.4),
      tickvals = c(2.2, 3.5, 4.2, 5.1),
      ticktext = c("Guides", "Gene", "Target Regions", "Variants"),
      showgrid = FALSE
    ),
    hovermode = "closest",
    showlegend = TRUE,
    legend = list(x = 1.02, y = 1)
  )
}

# =============================================================================
# MODE CLASSIFICATION AND STYLING
# =============================================================================

# Classify guides by mode for 5-color system
classify_guides_by_mode <- function(gene_guides) {
  
  mode_data <- list()
  
  # Tiling (Sea Green)
  mode_data[["Tiling"]] <- gene_guides[
    gene_guides$tiling_guide == "TRUE" & 
    gene_guides$crispri_guide == "FALSE" & 
    gene_guides$crispra_guide == "FALSE", 
  ]
  
  # CRISPRi (Crimson)
  mode_data[["CRISPRi"]] <- gene_guides[
    gene_guides$tiling_guide == "FALSE" & 
    gene_guides$crispri_guide == "TRUE" & 
    gene_guides$crispra_guide == "FALSE", 
  ]
  
  # CRISPRa (Dark Orchid)
  mode_data[["CRISPRa"]] <- gene_guides[
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
  
  # Remove empty modes
  mode_data <- mode_data[sapply(mode_data, nrow) > 0]
  
  return(mode_data)
}

# Create mode label for single guide
classify_guide_mode_label <- function(tiling, crispri, crispra) {
  case_when(
    tiling == "TRUE" & crispri == "TRUE" ~ "Tiling + CRISPRi",
    tiling == "TRUE" & crispra == "TRUE" ~ "Tiling + CRISPRa", 
    tiling == "TRUE" ~ "Tiling",
    crispri == "TRUE" ~ "CRISPRi",
    crispra == "TRUE" ~ "CRISPRa",
    TRUE ~ "Other"
  )
}

# Get color palette for modes
get_mode_color_palette <- function() {
  c(
    "Tiling" = "#2E8B57",            # Sea Green
    "CRISPRi" = "#DC143C",           # Crimson
    "CRISPRa" = "#9932CC",           # Dark Orchid
    "Tiling + CRISPRi" = "#FF6B35",  # Orange
    "Tiling + CRISPRa" = "#7B68EE",  # Medium Slate Blue
    "Other" = "#708090"              # Slate Gray
  )
}

# Get color for specific mode
get_mode_color <- function(mode) {
  colors <- get_mode_color_palette()
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

# Get Y position for each mode (better distributed spacing)
get_mode_y_position <- function(mode) {
  positions <- list(
    "Tiling" = 2.8,
    "CRISPRi" = 2.4,
    "CRISPRa" = 2.0,
    "Tiling + CRISPRi" = 1.6,
    "Tiling + CRISPRa" = 1.2
  )
  return(positions[[mode]] %||% 1.0)
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Create hover text for guides
create_guide_hover_text <- function(guides, mode) {
  paste(
    "Guide:", guides$target,
    "<br>Gene:", guides$gene,
    "<br>Hsu2013 Score:", round(guides$Hsu2013, 1),
    "<br>Mode:", mode,
    "<br>Position:", format(guides$abs_start, big.mark = ",")
  )
}
