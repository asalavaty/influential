#=============================================================================
#
#    Run shiny apps (examples)
#
#=============================================================================

#' @title Run shiny app
#' @description Run shiny apps included in the influential R package.
#' Also, a web-based \href{https://influential.erc.monash.edu/}{Influential Software Package} with a convenient
#' user-interface (UI) has been developed for the comfort of all users including those without a coding background.
#' @param shinyApp The name of the shiny app you want to run. You can get the exact name of the available
#' shiny apps via the following command.
#' \emph{list.files(system.file("ShinyApps", package = "influential"))}. Please also note this function is
#' case-sensitive.
#' @return A shiny app.
#' @keywords runShinyApp
#' @export runShinyApp
#' @examples
#' \dontrun{
#' runShinyApp(shinyApp = "IVI")
#' }
runShinyApp <- function(shinyApp) {
  
  cli::cli_alert_info("Checking the requirements for running the shiny app.\n")
  
  # Check if the BiocManager is installed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    utils::install.packages("BiocManager", dependencies = TRUE)
  }
  
  # Loop through each package
  for (pkg in c("shiny", "shinythemes", "shinyWidgets", "shinyjs",
                "shinycssloaders", "colourpicker", "DT", "magrittr", "janitor",
                "ranger", "influential", "ggplot2", "igraph")) {
    # Check if the package namespace is available
    if (!requireNamespace(pkg, quietly = TRUE)) {
      # Install the package if it's not available
      BiocManager::install(pkg, update = TRUE, ask = FALSE, quiet = TRUE)
      cli::cli_alert_info(paste0("The package '", pkg, "' required for running the shiny app is installed.\n"))
    }
  }
  
  # locate all the shiny app examples that exist
  validExamples <- list.files(system.file("ShinyApps", package = "influential"))
  
  validExamplesMsg <-
    paste0(
      "Valid shiny apps are: '",
      paste(validExamples, collapse = "', '"),
      "'")
  
  # if an invalid shiny app is given, throw an error
  if (missing(shinyApp) || !nzchar(shinyApp) || !shinyApp %in% validExamples) {
    cli::cli_abort(
      c(
        "Please run {.fn influential::runShinyApp} with a valid shiny app name as an argument.",
        "i" = validExamplesMsg
      ),
      call = FALSE
    )
  }
  
  cli::cli_alert_info("Loading the shiny app ...")
  
  # find and launch the app
  appDir <- system.file("ShinyApps", shinyApp, package = "influential")
  shiny::runApp(appDir, display.mode = "normal")
}