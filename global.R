required_packages <- c("shiny", "bslib", "DT", "plotly", "jsonlite", "htmltools")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Install the following R packages before starting REVIVE: ",
       paste(missing_packages, collapse = ", "), ".")
}

library(shiny)
library(bslib)

source("R/helpers_dashboard.R")
source("R/mod_dashboard.R")

# This prepared dataset is shared read-only across sessions.
revive_data <- read_revive_dashboard_data("data/preprocessed/longitudinal_events.json")
