# REVIVE

An R Shiny dashboard for **Respiratory Viral Infection Vaccine Effectiveness in
Transplant Recipients**, using the sidebar, summary cards, and content cards from
the [Posit dashboard template](https://shiny.posit.co/py/templates/dashboard-tips/)
as a layout reference.

From the project root, run:

```r
shiny::runApp(".")
```

Required R packages: `shiny`, `bslib`, `DT`, `plotly`, `jsonlite`, and `htmltools`.
The app reads `data/preprocessed/longitudinal_events.json` once on startup.
After updating the prepared JSON, restart the app. Data preparation is documented
in `R/00_data_preprocess.qmd`.

`global.R`, `ui.R`, and `server.R` coordinate the application. The first module,
`R/mod_dashboard.R`, owns the filters, patient count, table, CSV download, and
swimmer plot. Shared data and plotting functions are in `R/helpers_dashboard.R`;
the Arctic reflection palette is defined in `www/custom.css`.

Filters apply together at the patient level. A category matches any selection
or report for that patient. The full age range means **All**, including missing
ages; a narrower range includes only known ages within that range. Reset restores
all six filters. SOTR describes **prior** transplant organs from `SOTRs`, not
overall transplant status. The table and CSV include every matching patient;
the CSV includes all matching rows, regardless of the table's current page.

The swimmer plot shows all dated events for matching patients. Infection uses
red (`#a00000`) circles and vaccination uses teal (`#298c8c`) squares. Patients
without dated events remain in the count, table, and export. Lines connect first
and last recorded events and do not imply observed follow-up time. Multiple
reports at the same coordinates can overlap. The chart scrolls vertically and
supports hover details, zoom, and PNG download.

Run dashboard checks from the project root:

```sh
Rscript tests/test_dashboard.R
```
