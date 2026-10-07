revive_filter_choices <- function(values, missing_label = "Not recorded") {
  if (!is.list(values)) values <- as.list(values)
  levels <- sort(unique(as.character(unlist(values, use.names = FALSE))))
  levels <- levels[!is.na(levels) & nzchar(levels)]
  choices <- c("All" = "__all__", stats::setNames(levels, levels))
  missing <- vapply(values, function(value) !length(value) || all(is.na(value)), logical(1))
  if (any(missing)) choices <- c(choices, stats::setNames("__missing__", missing_label))
  choices
}

filter_revive_patients <- function(patients, filters = list(), age = NULL, age_limits = NULL) {
  keep <- rep(TRUE, nrow(patients))
  for (field in c("study_group", "infection", "vaccination", "sex", "sotr", "immunosuppression")) {
    selected <- filters[[field]]
    if (!length(selected) || "__all__" %in% selected) next
    values <- if (is.list(patients[[field]])) patients[[field]] else as.list(patients[[field]])
    matches <- vapply(values, function(value) {
      missing <- !length(value) || all(is.na(value))
      ("__missing__" %in% selected && missing) || any(selected %in% value)
    }, logical(1))
    keep <- keep & matches
  }
  # The full slider range means All, including patients with unknown age.
  age_is_all <- is.null(age) || identical(as.numeric(age), as.numeric(age_limits))
  if (!age_is_all) {
    keep <- keep & !is.na(patients$age) & patients$age >= age[1] & patients$age <= age[2]
  }
  patients[keep, , drop = FALSE]
}

revive_patient_table <- function(patients) {
  collapse_labels <- function(values, empty) {
    vapply(values, function(value) if (length(value)) paste(value, collapse = "; ") else empty, character(1))
  }
  data.frame(
    "Patient ID" = patients$record_id,
    "Study group" = ifelse(is.na(patients$study_group), "Not recorded", patients$study_group),
    Age = patients$age,
    Sex = ifelse(is.na(patients$sex), "Not recorded", patients$sex),
    Race = collapse_labels(patients$race, "Not recorded"),
    SOTR = collapse_labels(patients$sotr, "No transplant organs recorded"),
    Immunosuppression = collapse_labels(patients$immunosuppression, "Not recorded"),
    Infection = collapse_labels(patients$infection, "No infection reported"),
    Vaccination = collapse_labels(patients$vaccination, "No vaccination reported"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

revive_export_columns <- function() {
  # Edit this list to choose which extra variables users may download.
  c("ig_plasma", "ig_plasma_3mo")
}

revive_patient_export <- function(patients, patient_characteristics, extra_columns = character()) {
  export <- revive_patient_table(patients)
  extra_columns <- intersect(extra_columns, revive_export_columns())
  # Select only the requested columns, then align them with the filtered IDs.
  extra <- patient_characteristics[extra_columns]
  rows <- match(patients$record_id, patient_characteristics$record_id)
  export[extra_columns] <- extra[rows, , drop = FALSE]
  export
}

revive_swimmer_plot <- function(events) {
  events <- events[!is.na(events$event_date), , drop = FALSE]
  if (!nrow(events)) return(NULL)
  groups <- split(events$event_date, events$record_id)
  lanes <- data.frame(
    record_id = names(groups),
    first = as.Date(vapply(groups, function(dates) as.numeric(min(dates)), numeric(1)), origin = "1970-01-01"),
    last = as.Date(vapply(groups, function(dates) as.numeric(max(dates)), numeric(1)), origin = "1970-01-01")
  )
  lanes <- lanes[order(lanes$first, lanes$record_id), ]
  lanes$lane <- seq_len(nrow(lanes))
  events$lane <- lanes$lane[match(events$record_id, lanes$record_id)]
  events$hover <- paste0(
    "Patient ", htmltools::htmlEscape(events$record_id), "<br>",
    htmltools::htmlEscape(events$event_label), "<br>", format(events$event_date, "%d %b %Y")
  )
  plot <- plotly::plot_ly(height = max(340, nrow(lanes) * 25 + 100)) |>
    plotly::add_segments(
      data = lanes, x = ~first, xend = ~last, y = ~lane, yend = ~lane,
      line = list(color = "#ACBCBF", width = 2), hoverinfo = "skip", showlegend = FALSE
    )
  for (kind in c("Infection", "Vaccination")) {
    subset <- events[events$event_type == kind, , drop = FALSE]
    if (!nrow(subset)) next
    plot <- plotly::add_trace(
      plot, data = subset, x = ~event_date, y = ~lane, text = ~hover,
      type = "scatter", mode = "markers", name = kind,
      marker = list(
        color = if (kind == "Infection") "#a00000" else "#298c8c",
        symbol = if (kind == "Infection") "circle" else "square",
        size = if (kind == "Infection") 11 else 10,
        line = list(color = "#ffffff", width = 1)
      ),
      hoverinfo = "text", showlegend = FALSE, inherit = FALSE
    )
  }
  date_range <- range(events$event_date) + c(-7, 7)
  plot |>
    plotly::layout(
      font = list(family = "Arial, sans-serif", color = "#243C4C", size = 12),
      paper_bgcolor = "#ffffff", plot_bgcolor = "#ffffff",
      margin = list(l = 100, r = 25, t = 20, b = 65),
      xaxis = list(title = "Event date", type = "date", range = as.character(date_range),
                   tickformat = if (diff(date_range) <= 90) "%d %b" else "%b %Y",
                   dtick = if (diff(date_range) > 90) "M1" else NULL, nticks = 6,
                   gridcolor = "#edf2f3", zeroline = FALSE),
      yaxis = list(title = "Patient ID", tickvals = lanes$lane, ticktext = lanes$record_id,
                   range = c(nrow(lanes) + 0.7, 0.3), fixedrange = TRUE,
                   gridcolor = "#f4f7f8", zeroline = FALSE),
      hovermode = "closest"
    ) |>
    plotly::config(
      displaylogo = FALSE, responsive = TRUE,
      modeBarButtonsToRemove = c("select2d", "lasso2d"),
      toImageButtonOptions = list(format = "png", filename = "REVIVE_swimmer_plot")
    )
}
