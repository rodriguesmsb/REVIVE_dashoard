# Generate longitudinal event histories from the REDCap exports in data/.
# Dictionary: ref/Revive_dictionary.pdf, new_inf1 and new_vax_recvd.

read_revive_event_slots <- function(paths, type_fields, date_fields) {
  reports <- lapply(paths, function(path) {
    data <- readr::read_csv(
      path,
      col_types = readr::cols(.default = readr::col_character())
    )
    slots <- lapply(seq_along(type_fields), function(slot) {
      dplyr::transmute(
        data,
        record_id,
        event_code = .data[[type_fields[slot]]],
        event_date = .data[[date_fields[slot]]]
      )
    })
    dplyr::bind_rows(slots)
  })
  empty_reports <- data.frame(
    record_id = character(), event_code = character(), event_date = character()
  )
  dplyr::bind_rows(empty_reports, reports) |>
    dplyr::filter(!is.na(event_code) | !is.na(event_date)) |>
    dplyr::mutate(event_date = as.Date(event_date, format = "%Y-%m-%d"))
}

generate_revive_json <- function(data_dir = "data", output_dir = "data/preprocessed",
                                 record_ids = character(), skip_files = character()) {
  paths <- list.files(data_dir, full.names = TRUE)
  csv_files <- paths[endsWith(tolower(paths), ".csv")]
  excluded <- basename(csv_files) %in% skip_files
  if (any(excluded)) {
    warning("Skipping event exports: ", paste(basename(csv_files[excluded]), collapse = ", "),
            call. = FALSE)
  }
  csv_files <- csv_files[!excluded]
  filenames <- tolower(basename(csv_files))
  infection_files <- csv_files[grepl("infection", filenames, fixed = TRUE)]
  vaccine_files <- csv_files[grepl("vaccin", filenames, fixed = TRUE)]

  infections <- read_revive_event_slots(
    infection_files,
    type_fields = c("new_inf1", "new_inf2", "new_inf3", "new_inf4", "new_inf5", "new_inf6"),
    date_fields = c("new_inf1_dt", "new_inf2_dt", "new_inf3_dt", "new_inf4_dt", "new_inf5_dt", "new_inf6_dt")
  ) |>
    dplyr::transmute(
      record_id,
      infection = dplyr::case_when(
        event_code == 1 ~ "Flu infection",
        event_code == 2 ~ "COVID-19 infection",
        event_code == 3 ~ "RSV infection",
        .default = NA_character_
      ),
      infection_date = as.character(event_date)
    )

  vaccinations <- read_revive_event_slots(
    vaccine_files,
    type_fields = c("new_vax_recvd", "new2_vax_recvd", "new3_vax_recvd", "new4_vax_recvd", "new5_vax_recvd", "new6_vax_recvd"),
    date_fields = c("new_vax_dt", "new2_vax_dt", "new3_vax_dt", "new4_vax_dt", "new5_vax_dt", "new6_vax_dt")
  ) |>
    dplyr::transmute(
      record_id,
      vaccine_type = dplyr::case_when(
        event_code == 1 ~ "Flu vaccine",
        event_code == 2 ~ "COVID-19 vaccine",
        event_code == 3 ~ "Respiratory syncytial virus (RSV) vaccine",
        .default = NA_character_
      ),
      vaccine_date = as.character(event_date)
    )

  # Keep every reported slot. No event is discarded because another report
  # has the same date. Participants without reports have empty arrays.
  record_ids <- unique(c(record_ids, infections$record_id, vaccinations$record_id))
  people <- lapply(record_ids, function(id) {
    patient_infections <- infections[infections$record_id == id, ]
    patient_vaccinations <- vaccinations[vaccinations$record_id == id, ]
    list(
      record_id = id,
      infections = lapply(seq_len(nrow(patient_infections)), function(row) {
        list(infection = patient_infections$infection[row],
             infection_date = patient_infections$infection_date[row])
      }),
      vaccinations = lapply(seq_len(nrow(patient_vaccinations)), function(row) {
        list(vaccine_type = patient_vaccinations$vaccine_type[row],
             vaccine_date = patient_vaccinations$vaccine_date[row])
      })
    )
  })

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(
    people, file.path(output_dir, "longitudinal_events.json"),
    auto_unbox = TRUE, pretty = TRUE, na = "null", null = "null"
  )
  people
}
