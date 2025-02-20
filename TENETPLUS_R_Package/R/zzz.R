.onAttach <- function(libname, pkgname) {
  # Load the crayon package for colors
  if (!requireNamespace("crayon", quietly = TRUE)) {
    return()  # Avoid errors if crayon is missing
  }

  # Define colored ASCII logo using crayon functions
  logo <- paste(
    crayon::cyan(""),
    crayon::black("  ████████╗ ███████╗ ███╗   ██╗ ███████╗ ████████╗     ██╗       "),
    crayon::black("  ╚══██╔══╝ ██╔════╝ ████╗  ██║ ██╔════╝ ╚══██╔══╝     ██║         "),
    crayon::black("     ██║    █████╗   ██╔██╗ ██║ █████╗      ██║    ██████████╗"),
    crayon::black("     ██║    ██╔══╝   ██║╚██╗██║ ██╔══╝      ██║    ╚═══██╔═══╝"),
    crayon::black("     ██║    ███████╗ ██║ ╚████║ ███████╗    ██║        ██║"),
    crayon::black("     ╚═╝    ╚══════╝ ╚═╝  ╚═══╝ ╚══════╝    ╚═╝        ╚═╝"),
    sep = "\n"
  )

  # Print the logo on package load
  packageStartupMessage(logo)
}



.onLoad <- function(libname, pkgname) {
  # Set the bash script path inside the package
  bashFile <- system.file("bash", "TENET_Plus_for_py", package = "TENETPLUS")
  if (file.exists(bashFile)) {
        # Set executable permission
        system(paste("chmod +x", shQuote(bashFile)))
  } else {
        warning("TENET_Plus_for_py script not found in package.")
  }
  if (!file.exists(bashFile) || bashFile == "") {
    stop("Error: TENETPLUS bash script not found! Ensure 'TENET_Plus_for_py' is inside 'inst/bash/'.")
  }

  options(TENETPLUS.bashFile = bashFile)

  # Check and install missing dependencies
  required_packages <- c("ComplexHeatmap", "circlize")

  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    message("Missing required packages: ", paste(missing_packages, collapse = ", "))
    
    if ("ComplexHeatmap" %in% missing_packages) {
      message("Installing 'ComplexHeatmap' from Bioconductor...")
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install("ComplexHeatmap")
    }
    
    if ("circlize" %in% missing_packages) {
      install.packages("circlize", dependencies = TRUE)
    }
  }
}

