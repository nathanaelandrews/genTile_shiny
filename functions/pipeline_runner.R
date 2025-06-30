# Pipeline Runner Functions for GenTile Shiny App
# functions/pipeline_runner.R

# Set genTile repository path
GENTILE_PATH <- "/Users/nathanaelandrews/wrk/github/genTile"

# Main pipeline function
run_gentile_pipeline <- function(input_method, input_data, parameters, progress_callback = NULL) {
  
  # Create session-specific temporary directory
  session_id <- paste0("gentile_", Sys.getpid(), "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  temp_dir <- file.path(tempdir(), session_id)
  dir.create(temp_dir, recursive = TRUE)
  
  # Ensure cleanup on exit
  on.exit({
    unlink(temp_dir, recursive = TRUE)
  })
  
  if (!is.null(progress_callback)) {
    progress_callback("Preparing input data...", 5)
  }
  
  # Step 1: Prepare input file
  input_file <- prepare_input_data(input_method, input_data, temp_dir)
  
  if (!is.null(progress_callback)) {
    progress_callback("Extracting genomic sequences...", 15)
  }
  
  # Step 2: Extract sequences
  sequences_file <- extract_sequences(input_file, input_method, parameters, temp_dir)
  
  if (!is.null(progress_callback)) {
    progress_callback("Designing guides with FlashFry...", 35)
  }
  
  # Step 3: Design guides
  guides_file <- design_guides(sequences_file, temp_dir)
  
  if (!is.null(progress_callback)) {
    progress_callback("Selecting optimal guides...", 70)
  }
  
  # Step 4: Select guides
  selected_guides <- select_guides(guides_file, parameters, temp_dir)
  
  if (!is.null(progress_callback)) {
    progress_callback("Processing results...", 90)
  }
  
  # Step 5: Process and return results
  results <- process_results(selected_guides, temp_dir)
  
  if (!is.null(progress_callback)) {
    progress_callback("Complete!", 100)
  }
  
  return(results)
}

# Helper function to prepare input data
prepare_input_data <- function(input_method, input_data, temp_dir) {
  
  input_file <- file.path(temp_dir, "input.txt")
  
  if (input_method == "paste_genes") {
    # Parse pasted gene text
    genes <- trimws(strsplit(input_data$gene_text, "\n")[[1]])
    genes <- genes[nchar(genes) > 0]
    
    if (length(genes) == 0) {
      stop("No valid genes provided")
    }
    
    if (length(genes) > 50) {
      stop("Too many genes. Maximum is 50 genes per analysis.")
    }
    
    writeLines(genes, input_file)
    
  } else if (input_method == "upload_genes") {
    # Handle uploaded gene file
    if (is.null(input_data$gene_file)) {
      stop("No gene file uploaded")
    }
    
    genes <- readLines(input_data$gene_file$datapath)
    genes <- trimws(genes)
    genes <- genes[nchar(genes) > 0 & !grepl("^#", genes)]
    
    if (length(genes) == 0) {
      stop("No valid genes found in uploaded file")
    }
    
    if (length(genes) > 50) {
      stop("Too many genes. Maximum is 50 genes per analysis.")
    }
    
    writeLines(genes, input_file)
    
  } else if (input_method == "paste_coords") {
    # Parse pasted coordinates
    coords <- trimws(strsplit(input_data$coord_text, "\n")[[1]])
    coords <- coords[nchar(coords) > 0]
    
    if (length(coords) == 0) {
      stop("No valid coordinates provided")
    }
    
    if (length(coords) > 50) {
      stop("Too many coordinates. Maximum is 50 entries per analysis.")
    }
    
    # Validate coordinate format
    valid_pattern <- "^[^,]+,chr[^:]+:[0-9]+,[+-]$"
    invalid_coords <- coords[!grepl(valid_pattern, coords)]
    
    if (length(invalid_coords) > 0) {
      stop(paste("Invalid coordinate format. Expected: name,chr:pos,strand\nInvalid entries:", 
                paste(invalid_coords[1:min(3, length(invalid_coords))], collapse = ", ")))
    }
    
    writeLines(coords, input_file)
    
  } else if (input_method == "upload_coords") {
    # Handle uploaded coordinates file
    if (is.null(input_data$coord_file)) {
      stop("No coordinates file uploaded")
    }
    
    coords <- readLines(input_data$coord_file$datapath)
    coords <- trimws(coords)
    coords <- coords[nchar(coords) > 0 & !grepl("^#", coords)]
    
    if (length(coords) == 0) {
      stop("No valid coordinates found in uploaded file")
    }
    
    if (length(coords) > 50) {
      stop("Too many coordinates. Maximum is 50 entries per analysis.")
    }
    
    # Validate coordinate format
    valid_pattern <- "^[^,]+,chr[^:]+:[0-9]+,[+-]$"
    invalid_coords <- coords[!grepl(valid_pattern, coords)]
    
    if (length(invalid_coords) > 0) {
      stop(paste("Invalid coordinate format. Expected: name,chr:pos,strand\nInvalid entries:", 
                paste(invalid_coords[1:min(3, length(invalid_coords))], collapse = ", ")))
    }
    
    writeLines(coords, input_file)
    
  } else {
    stop("Unknown input method")
  }
  
  return(input_file)
}

# Extract genomic sequences
extract_sequences <- function(input_file, input_method, parameters, temp_dir) {
  
  # Determine if using genes or coordinates
  is_gene_input <- grepl("genes", input_method)
  
  # Build get_sequence.sh command
  cmd <- file.path(GENTILE_PATH, "bin/get_sequence.sh")
  
  sequences_file <- file.path(temp_dir, "sequences.fa")
  
  if (is_gene_input) {
    # Gene mode
    args <- c(
      "-g", shQuote(input_file),
      "-r", shQuote(file.path(GENTILE_PATH, "data/reference/genome/hg38/hg38.fa")),
      "-c", parameters$cell_line,
      "-u", parameters$upstream,
      "-d", parameters$downstream,
      "-o", shQuote(sequences_file)
    )
  } else {
    # Position mode
    args <- c(
      "-p", shQuote(input_file),
      "-r", shQuote(file.path(GENTILE_PATH, "data/reference/genome/hg38/hg38.fa")),
      "-u", parameters$upstream,
      "-d", parameters$downstream,
      "-o", shQuote(sequences_file)
    )
  }
  
  # Execute command
  result <- system2(cmd, args, stdout = TRUE, stderr = TRUE)
  
  if (!is.null(attr(result, "status")) && attr(result, "status") != 0) {
    stop(paste("Sequence extraction failed:", paste(result, collapse = "\n")))
  }
  
  if (!file.exists(sequences_file)) {
    stop("Sequence extraction failed: output file not created")
  }
  
  return(sequences_file)
}

# Design guides using FlashFry
design_guides <- function(sequences_file, temp_dir) {
  
  cmd <- file.path(GENTILE_PATH, "bin/design_guides.sh")
  guides_file <- file.path(temp_dir, "guides.scored.txt")
  
  args <- c(
    "-i", shQuote(sequences_file),
    "-o", shQuote(guides_file)
  )
  
  # Execute command
  result <- system2(cmd, args, stdout = TRUE, stderr = TRUE)
  
  if (!is.null(attr(result, "status")) && attr(result, "status") != 0) {
    stop(paste("Guide design failed:", paste(result, collapse = "\n")))
  }
  
  if (!file.exists(guides_file)) {
    stop("Guide design failed: output file not created")
  }
  
  return(guides_file)
}

# Select optimal guides
select_guides <- function(guides_file, parameters, temp_dir) {
  
  cmd <- file.path(GENTILE_PATH, "bin/select_guides.sh")
  selected_file <- file.path(temp_dir, "selected_guides.txt")
  
  # Build arguments based on selected modes
  args <- c(
    "-f", shQuote(guides_file),
    "-o", shQuote(selected_file)
  )
  
  # Add mode flags
  if ("tiling" %in% parameters$selection_modes) {
    args <- c(args, "-t")
    args <- c(args, "-z", parameters$zone_size)
  }
  
  if ("crispri" %in% parameters$selection_modes) {
    args <- c(args, "-i")
    args <- c(args, "-n", parameters$target_guides)
  }
  
  if ("crispra" %in% parameters$selection_modes) {
    args <- c(args, "-a")
    args <- c(args, "-n", parameters$target_guides)
  }
  
  # Add filtering options
  if (parameters$use_k562_variants) {
    k562_variants <- file.path(GENTILE_PATH, "data/reference/vcf/K562_variants_hg38_sorted.bed")
    if (file.exists(k562_variants)) {
      args <- c(args, "-V", shQuote(k562_variants))
    }
  }
  
  if (parameters$filter_restriction && length(parameters$restriction_enzymes) > 0) {
    enzyme_list <- paste(parameters$restriction_enzymes, collapse = ",")
    args <- c(args, "-R", enzyme_list)
  }
  
  # Execute command
  result <- system2(cmd, args, stdout = TRUE, stderr = TRUE)
  
  if (!is.null(attr(result, "status")) && attr(result, "status") != 0) {
    stop(paste("Guide selection failed:", paste(result, collapse = "\n")))
  }
  
  if (!file.exists(selected_file)) {
    stop("Guide selection failed: output file not created")
  }
  
  return(selected_file)
}

# Process results and return structured data
process_results <- function(selected_guides_file, temp_dir) {
  
  # Read the selected guides
  guides_data <- read.delim(selected_guides_file, stringsAsFactors = FALSE)
  
  # Read the BED file if it exists
  bed_file <- file.path(temp_dir, "selected_guides.bed")
  bed_data <- NULL
  if (file.exists(bed_file)) {
    bed_data <- read.delim(bed_file, header = FALSE, stringsAsFactors = FALSE)
    names(bed_data) <- c("chr", "start", "end", "guide_id", "score", "strand")
  }
  
  # Create color assignments for guides based on modes
  guides_data$guide_color <- assign_guide_colors(guides_data)
  
  # Extract target/gene information for grouping
  guides_data$gene <- extract_gene_names(guides_data)
  
  results <- list(
    guides = guides_data,
    bed_data = bed_data,
    summary = create_summary_stats(guides_data)
  )
  
  return(results)
}

# Extract proper gene names from the data
extract_gene_names <- function(guides_data) {
  
  # The FlashFry output has gene names in the 'contig' column
  # Format: GENENAME::chr:pos-pos(strand)::TSS_info
  # We want to extract just the gene name (first part before first ::)
  if ("contig" %in% names(guides_data)) {
    # Extract gene name from contig (everything before first ::)
    gene_names <- sub("::.*", "", guides_data$contig)
  } else if ("guide_id" %in% names(guides_data)) {
    # Fallback: Extract gene name from guide_id (everything before first underscore)
    gene_names <- sub("_.*", "", guides_data$guide_id)
  } else if ("target" %in% names(guides_data)) {
    # Last resort - but this would be wrong as target is the gRNA sequence
    gene_names <- guides_data$target
  } else {
    # Ultimate fallback
    gene_names <- paste0("Target_", seq_len(nrow(guides_data)))
  }
  
  return(gene_names)
}

# Assign colors based on guide modes (5-color system)
assign_guide_colors <- function(guides_data) {
  
  colors <- character(nrow(guides_data))
  
  for (i in 1:nrow(guides_data)) {
    tiling <- guides_data$tiling_guide[i] == "TRUE"
    crispri <- guides_data$crispri_guide[i] == "TRUE"  
    crispra <- guides_data$crispra_guide[i] == "TRUE"
    
    if (tiling && crispri) {
      colors[i] <- "#FF6B35"  # Orange - Tiling + CRISPRi
    } else if (tiling && crispra) {
      colors[i] <- "#7B68EE"  # Medium Slate Blue - Tiling + CRISPRa
    } else if (tiling) {
      colors[i] <- "#2E8B57"  # Sea Green - Tiling only
    } else if (crispri) {
      colors[i] <- "#DC143C"  # Crimson - CRISPRi only
    } else if (crispra) {
      colors[i] <- "#9932CC"  # Dark Orchid - CRISPRa only
    } else {
      colors[i] <- "#708090"  # Slate Gray - No mode (shouldn't happen)
    }
  }
  
  return(colors)
}

# Create summary statistics
create_summary_stats <- function(guides_data) {
  
  total_guides <- nrow(guides_data)
  tiling_count <- sum(guides_data$tiling_guide == "TRUE", na.rm = TRUE)
  crispri_count <- sum(guides_data$crispri_guide == "TRUE", na.rm = TRUE)
  crispra_count <- sum(guides_data$crispra_guide == "TRUE", na.rm = TRUE)
  
  # Count overlapping modes
  tiling_crispri <- sum(guides_data$tiling_guide == "TRUE" & guides_data$crispri_guide == "TRUE", na.rm = TRUE)
  tiling_crispra <- sum(guides_data$tiling_guide == "TRUE" & guides_data$crispra_guide == "TRUE", na.rm = TRUE)
  
  # Count unique genes/targets
  unique_targets <- length(unique(guides_data$gene))
  
  # Score statistics (assuming Hsu2013 column exists)
  score_stats <- NULL
  if ("Hsu2013" %in% names(guides_data)) {
    score_stats <- list(
      mean_score = mean(guides_data$Hsu2013, na.rm = TRUE),
      median_score = median(guides_data$Hsu2013, na.rm = TRUE),
      min_score = min(guides_data$Hsu2013, na.rm = TRUE),
      max_score = max(guides_data$Hsu2013, na.rm = TRUE)
    )
  }
  
  summary <- list(
    total_guides = total_guides,
    unique_targets = unique_targets,
    mode_counts = list(
      tiling = tiling_count,
      crispri = crispri_count,
      crispra = crispra_count,
      tiling_crispri = tiling_crispri,
      tiling_crispra = tiling_crispra
    ),
    score_stats = score_stats
  )
  
  return(summary)
}
