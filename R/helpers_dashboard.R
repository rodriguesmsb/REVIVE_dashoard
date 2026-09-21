read_revive_dashboard_data <- function(path) {
  if (!file.exists(path)) {
    stop("Prepared data are missing. Run R/00_data_preprocess.qmd before starting REVIVE.")
  }
  people <- tryCatch(
    jsonlite::read_json(path, simplifyVector = FALSE),
    error = function(e) stop("The prepared JSON could not be read. Rerun R/00_data_preprocess.qmd.")
  )
  required <- c("record_id", "Sex", "Age", "SOTRs", "Immunosuppression", "infections", "vaccinations")
  if (!is.list(people) || !length(people) || !all(vapply(people, function(person) {
    is.list(person) && all(required %in% names(person)) &&
      is.character(person$record_id) && length(person$record_id) == 1L &&
      !is.na(person$record_id) && nzchar(person$record_id)
  }, logical(1)))) {
    stop("The prepared JSON must contain participant attributes and event histories. Rerun preprocessing.")
  }
  labels <- function(value) {
    value <- as.character(unlist(value, use.names = FALSE))
    sort(unique(value[!is.na(value) & nzchar(value)]))
  }
  event_labels <- function(reports, field) {
    labels(lapply(reports, function(report) {
      if (is.null(report[[field]]) || !nzchar(report[[field]])) "Type not recorded" else report[[field]]
    }))
  }
  patients <- data.frame(
    record_id = vapply(people, `[[`, character(1), "record_id"),
    sex = vapply(people, function(person) {
      if (is.null(person$Sex)) NA_character_ else person$Sex
    }, character(1)),
    age = vapply(people, function(person) {
      if (is.null(person$Age)) NA_real_ else as.numeric(person$Age)
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
  if (anyDuplicated(patients$record_id)) stop("The prepared JSON contains duplicate participant IDs.")
  patients$infection <- lapply(people, function(person) event_labels(person$infections, "infection"))
  patients$vaccination <- lapply(people, function(person) event_labels(person$vaccinations, "vaccine_type"))
  patients$sotr <- lapply(people, function(person) labels(person$SOTRs))
  patients$immunosuppression <- lapply(people, function(person) labels(person$Immunosuppression))

  event_rows <- lapply(people, function(person) {
    rows <- lapply(c("infections", "vaccinations"), function(kind) {
      if (!length(person[[kind]])) return(NULL)
      is_infection <- kind == "infections"
      type_field <- if (is_infection) "infection" else "vaccine_type"
      date_field <- if (is_infection) "infection_date" else "vaccine_date"
      do.call(rbind, lapply(person[[kind]], function(report) {
        raw_date <- report[[date_field]]
        date <- if (is.null(raw_date)) as.Date(NA) else as.Date(raw_date, format = "%Y-%m-%d")
        data.frame(
          record_id = person$record_id,
          event_type = if (is_infection) "Infection" else "Vaccination",
          event_label = if (is.null(report[[type_field]])) "Type not recorded" else report[[type_field]],
          event_date = date, stringsAsFactors = FALSE
        )
      }))
    })
    do.call(rbind, rows)
  })
  events <- do.call(rbind, event_rows)
  if (is.null(events)) {
    events <- data.frame(record_id = character(), event_type = character(),
                         event_label = character(), event_date = as.Date(character()))
  }
  ages <- patients$age[is.finite(patients$age)]
  age_limits <- if (length(ages)) c(floor(min(ages)), ceiling(max(ages))) else c(0, 100)
  if (diff(age_limits) == 0) age_limits[2] <- age_limits[2] + 1
  list(patients = patients, events = events, age_limits = age_limits)
}

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
  for (field in c("infection", "vaccination", "sex", "sotr", "immunosuppression")) {
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
    Infection = collapse_labels(patients$infection, "No infection reported"),
    Vaccination = collapse_labels(patients$vaccination, "No vaccination reported"),
    Sex = ifelse(is.na(patients$sex), "Not recorded", patients$sex),
    SOTR = collapse_labels(patients$sotr, "No prior organs recorded"),
    Immunosuppression = collapse_labels(patients$immunosuppression, "Not recorded"),
    Age = patients$age, check.names = FALSE, stringsAsFactors = FALSE
  )
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
