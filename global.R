required_packages <- c("shiny", "bslib", "DT", "plotly", "jsonlite", "htmltools", "readr", "dplyr")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Install the following R packages before starting REVIVE: ",
       paste(missing_packages, collapse = ", "), ".")
}

library(shiny)
library(bslib)

source("R/REVIVE_harmonization.R")
source("R/json_generation.R")
source("R/RDS_generation.R")
source("R/helpers_dashboard.R")
source("R/mod_dashboard.R")

# Generate once at startup and share the prepared dataset across sessions.
revive_data <- read_revive_dashboard_data(
  "data/Patients_characterisitics.csv",
  data_dir = "data",
  # Remove this exclusion after replacing the unreadable, 8-byte export.
  skip_files = "Flu_infections.csv"
)
