.onAttach <- function(libname, pkgname) {
  # Load the crayon package for colors
  if (!requireNamespace("crayon", quietly = TRUE)) {
    return()  # Avoid errors if crayon is missing
  }
  
  # Define colored ASCII logo using crayon functions
  logo <- paste(
    crayon::black("     ")
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
  options(TENETPLUS.bash_script = system.file("bash", "TENET_Plus_for_py", package = "TENETPLUS"))
}

