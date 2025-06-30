#!/usr/bin/env Rscript

# Run GenTile Shiny App
# Usage: Rscript run_app.R [port]

# Get port from command line (default 3838)
args <- commandArgs(trailingOnly = TRUE)
port <- if(length(args) > 0) as.numeric(args[1]) else 3838

# Run the app
shiny::runApp("app.R", host = "127.0.0.1", port = port, launch.browser = TRUE)
