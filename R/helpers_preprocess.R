# Labels transcribed from ref/RespiratoryViralInfectionVacci.pdf:
# printed page 120 (infections) and page 79 (vaccinations).
revive_event_dictionary <- function() {
  list(
    infection = c(
      `1` = "Flu infection",
      `2` = "COVID-19 infection",
      `3` = "RSV infection"
    ),
    vaccination = c(
      `1` = "Flu vaccine",
      `2` = "COVID-19 vaccine",
      `3` = "Respiratory syncytial virus (RSV) vaccine"
    )
  )
}

parse_revive_date <- function(value) {
  result <- rep(as.Date(NA), length(value))
  formats <- c("%Y-%m-%d", "%m/%d/%Y", "%m-%d-%Y")
  for (date_format in formats) {
    parsed <- suppressWarnings(as.Date(value, format = date_format))
    # Round-trip validation rejects impossible dates and trailing text.
    valid <- !is.na(value) & !is.na(parsed) &
      format(parsed, date_format) == value
    result[valid] <- parsed[valid]
  }
  result
}

read_revive_events <- function(path, source_file = basename(path)) {
  raw <- readr::read_csv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    na = "", trim_ws = TRUE, name_repair = "minimal",
    show_col_types = FALSE, progress = FALSE
  )
  if (nrow(readr::problems(raw)) > 0L || anyDuplicated(names(raw))) {
    stop("Check malformed rows or duplicate column names in ", source_file, ".")
  }
  infection_fields <- grep("^new_inf[1-9][0-9]*(_dt)?$", names(raw), value = TRUE)
  vaccination_fields <- grep(
    "^new([1-9][0-9]*)?_vax_(recvd|dt)$", names(raw), value = TRUE
  )
  type_fields <- unique(c(
    sub("_dt$", "", infection_fields),
    sub("_dt$", "_recvd", vaccination_fields)
  ))
  if (!length(type_fields)) {
    message("Skipping CSV without infection/vaccination fields: ", source_file)
    return(NULL)
  }
  required <- c("record_id", "redcap_event_name")
  if (!all(required %in% names(raw))) {
    stop("Export record_id and redcap_event_name with the event fields in ", source_file, ".")
  }
  for (column in names(raw)) {
    raw[[column]] <- trimws(raw[[column]])
    raw[[column]][raw[[column]] == ""] <- NA_character_
  }
  for (column in c("redcap_repeat_instrument", "redcap_repeat_instance")) {
    if (!column %in% names(raw)) raw[[column]] <- NA_character_
  }

  dictionary <- revive_event_dictionary()
  pieces <- lapply(type_fields, function(type_field) {
    is_infection <- startsWith(type_field, "new_inf")
    event_type <- if (is_infection) "infection" else "vaccination"
    date_field <- if (is_infection) {
      paste0(type_field, "_dt")
    } else {
      sub("_recvd$", "_dt", type_field)
    }
    if (!all(c(type_field, date_field) %in% names(raw))) {
      stop("Export both ", type_field, " and ", date_field, " in ", source_file, ".")
    }
    keep <- !is.na(raw[[type_field]]) | !is.na(raw[[date_field]])
    rows <- which(keep)
    if (!length(rows)) return(NULL)
    if (anyNA(raw$record_id[rows]) || anyNA(raw$redcap_event_name[rows])) {
      stop("An event is missing record_id or redcap_event_name in ", source_file, ".")
    }
    slot <- if (is_infection) {
      sub("^new_inf", "", type_field)
    } else {
      sub("^new([0-9]*)_vax_recvd$", "\\1", type_field)
    }
    if (slot == "") slot <- "1"
    code <- raw[[type_field]][rows]
    date_raw <- raw[[date_field]][rows]
    date <- parse_revive_date(date_raw)
    label <- unname(dictionary[[event_type]][code])
    issues <- vapply(seq_along(rows), function(i) {
      flags <- c(
        missing_type = is.na(code[i]),
        unknown_type_code = !is.na(code[i]) && is.na(label[i]),
        missing_date = is.na(date_raw[i]),
        invalid_date = !is.na(date_raw[i]) && is.na(date[i])
      )
      paste(names(flags)[flags], collapse = ";")
    }, character(1))
    data.frame(
      record_id = raw$record_id[rows],
      redcap_event_name = raw$redcap_event_name[rows],
      redcap_repeat_instrument = raw$redcap_repeat_instrument[rows],
      redcap_repeat_instance = raw$redcap_repeat_instance[rows],
      event_type = event_type, event_label = label, event_date = date,
      type_code = code, event_date_raw = date_raw,
      source_slot = as.integer(slot), source_file = source_file,
      source_row = rows, quality_issue = issues,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, pieces)
}

revive_participant_histories <- function(events) {
  groups <- split(seq_len(nrow(events)), events$record_id)
  unname(lapply(groups, function(rows) {
    participant <- events[rows, , drop = FALSE]
    infection_rows <- which(participant$event_type == "infection")
    vaccination_rows <- which(participant$event_type == "vaccination")
    list(
      record_id = participant$record_id[1],
      infections = lapply(infection_rows, function(i) {
        list(
          infection = participant$event_label[i],
          infection_date = as.character(participant$event_date[i])
        )
      }),
      vaccinations = lapply(vaccination_rows, function(i) {
        list(
          vaccine_type = participant$event_label[i],
          vaccine_date = as.character(participant$event_date[i])
        )
      })
    )
  }))
}

preprocess_revive_data <- function(data_dir = "data",
                                   output_dir = file.path(data_dir, "processed"),
                                   participant_data = NULL) {
  for (package in c("readr", "jsonlite")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Install the required R package: ", package, ".")
    }
  }
  if (!dir.exists(data_dir)) stop("Source data directory does not exist: ", data_dir)
  data_dir <- normalizePath(data_dir, mustWork = TRUE)
  paths <- sort(list.files(data_dir, pattern = "\\.csv$", full.names = TRUE,
                           recursive = TRUE, ignore.case = TRUE))
  relative_paths <- substring(paths, nchar(data_dir) + 2L)
  # Process source CSVs only; derived files never become inputs.
  output_path <- normalizePath(output_dir, mustWork = FALSE)
  keep <- !startsWith(relative_paths, "processed/") &
    !startsWith(relative_paths, "preprocessed/") &
    !startsWith(paths, paste0(output_path, "/"))
  paths <- paths[keep]
  relative_paths <- relative_paths[keep]
  if (!length(paths)) stop("No source CSV files found in ", data_dir, ".")
  pieces <- Map(read_revive_events, paths, relative_paths)
  events <- do.call(rbind, pieces)
  if (is.null(events) || !nrow(events)) {
    stop("No infection or vaccination reports found; existing outputs were left unchanged.")
  }

  report_key <- c(
    "record_id", "redcap_event_name", "redcap_repeat_instrument",
    "redcap_repeat_instance", "event_type", "source_slot"
  )
  # Identical copies of a REDCap report in overlapping exports count once.
  # Different monthly reports and repeat instances remain separate observations.
  report_values <- c(report_key, "type_code", "event_date_raw")
  events <- events[!duplicated(events[report_values]), , drop = FALSE]
  if (anyDuplicated(events[report_key])) {
    stop(
      "Conflicting exports contain different values for the same REDCap report. ",
      "Reconcile the source CSVs before rebuilding; existing outputs were left unchanged."
    )
  }
  events <- events[order(
    events$record_id, events$event_date, events$event_type,
    events$redcap_event_name, events$redcap_repeat_instrument,
    events$redcap_repeat_instance, events$source_slot, na.last = TRUE
  ), , drop = FALSE]
  rownames(events) <- NULL

  events$same_event_report_count <- 1L
  complete <- !is.na(events$event_label) & !is.na(events$event_date)
  if (any(complete)) {
    # JSON-encoded tuples avoid collisions from delimiters inside identifiers.
    keys <- vapply(which(complete), function(i) {
      as.character(jsonlite::toJSON(unname(as.list(events[i, c(
        "record_id", "event_type", "type_code", "event_date"
      )])), auto_unbox = TRUE, Date = "ISO8601"))
    }, character(1))
    counts <- table(keys)
    events$same_event_report_count[complete] <- as.integer(counts[keys])
  }
  events$possible_duplicate <- events$same_event_report_count > 1L

  participants <- revive_participant_histories(events)
  expected_ids <- sort(unique(events$record_id))
  if (!is.null(participant_data)) {
    participants <- join_revive_participants(participants, participant_data)
    expected_ids <- sort(union(
      expected_ids, vapply(participant_data, `[[`, character(1), "record_id")
    ))
  }
  json <- jsonlite::toJSON(
    participants, auto_unbox = TRUE, pretty = TRUE, na = "null", null = "null"
  )
  decoded <- jsonlite::fromJSON(json, simplifyVector = FALSE)
  decoded_ids <- vapply(decoded, `[[`, character(1), "record_id")
  event_count <- sum(vapply(decoded, function(person) {
    length(person$infections) + length(person$vaccinations)
  }, integer(1)))
  stopifnot(
    jsonlite::validate(json),
    identical(decoded_ids, expected_ids),
    event_count == nrow(events)
  )

  json_path <- file.path(output_dir, "longitudinal_events.json")
  if (file.exists(json_path)) {
    previous <- jsonlite::read_json(json_path, simplifyVector = FALSE)
    expected <- c("record_id", "infections", "vaccinations")
    allowed <- c(expected, if (!is.null(participant_data)) {
      c("Sex", "Age", "Race", "Immunosuppression", "SOTRs")
    })
    compatible <- is.list(previous) && all(vapply(previous, function(person) {
      is.list(person) && all(expected %in% names(person)) &&
        all(names(person) %in% allowed)
    }, logical(1)))
    if (!compatible) stop("Existing JSON uses a different schema; it was left unchanged.")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  # Validate before replacing the JSON; rename a completed file on the same disk.
  temporary_json <- tempfile("longitudinal_events_", tmpdir = output_dir, fileext = ".json")
  on.exit(unlink(temporary_json), add = TRUE)
  writeLines(json, temporary_json, useBytes = TRUE)
  if (!file.rename(temporary_json, json_path)) stop("Could not replace ", json_path, ".")
  saveRDS(events, file.path(output_dir, "longitudinal_events.rds"))
  message(
    "Prepared ", nrow(events), " reports for ", length(participants), " participants. ",
    "Reports with quality issues: ", sum(events$quality_issue != ""), "; ",
    "possible duplicate reports: ", sum(events$possible_duplicate), "."
  )
  invisible(events)
}
