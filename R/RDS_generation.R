# Read source data, harmonize patients, and prepare dashboard inputs.
# Dictionary conversions belong in REVIVE_harmonization.R.

generate_revive_patients <- function(patients_path, output_dir = "data/preprocessed") {
  patients <- readr::read_csv(
    "data/Patients_characterisitics.csv",
    col_types = readr::cols(.default = readr::col_character())
  )
  patients_harmonized <- harmonize_patients(patients)

  # remove patients with no enrrolment data
  patients_harmonized <- patients_harmonized |> 
    dplyr::filter(!is.na(enrollment_dt))

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(
    patients_harmonized,
    file.path(output_dir, "patients_characteristics_harmonized.rds")
  )
  readr::write_csv(
    patients_harmonized,
    file.path(output_dir, "patients_characteristics_harmonized.csv"),
    na = ""
  )

  patients_harmonized
}

revive_selected_labels <- function(patients, fields) {
  selections <- patients[unname(fields)]
  lapply(seq_len(nrow(patients)), function(row) {
    values <- unlist(selections[row, ], use.names = FALSE)
    names(fields)[which(values == "Yes")]
  })
}

generate_revive_dashboard_data <- function(patients_harmonized, events) {
  race_fields <- c(
    White = "race_white",
    "Black or African American" = "race_black_african_american",
    Asian = "race_asian",
    "American Indian or Alaska Native" = "race_american_indian_alaska_native",
    "Native Hawaiian or Other Pacific Islander" = "race_native_hawaiian_other_pacific_islander",
    Other = "race_other",
    "Prefer not to answer" = "race_prefer_not_to_answer"
  )
  organ_fields <- c(
    Kidney = "kidney_transplant",
    Liver = "liver_transplant",
    Lung = "lung_transplant",
    Heart = "heart_transplant",
    Pancreas = "pancreas_transplant",
    Intestine = "intestine_transplant",
    "Islet cell" = "islet_cell_transplant",
    Other = "other_transplant"
  )
  medication_fields <- c(
    None = "no_immunosuppression",
    Azathioprine = "azathioprine",
    Belatacept = "belatacept",
    Cyclosporine = "cyclosporine",
    Tacrolimus = "tacrolimus",
    Everolimus = "everolimus",
    Sirolimus = "sirolimus",
    "Mycophenolate mofetil or Mycophenolic acid" = "mycophenolate",
    "Prednisone or Methylprednisolone" = "prednisone_methylprednisolone",
    "Other immunosuppression" = "other_immunosuppression"
  )

  patients <- data.frame(
    record_id = patients_harmonized$record_id,
    study_group = patients_harmonized$study_group,
    sex = patients_harmonized$sex_at_birth,
    age = patients_harmonized$age,
    stringsAsFactors = FALSE
  )
  patients$race <- revive_selected_labels(patients_harmonized, race_fields)
  patients$sotr <- revive_selected_labels(patients_harmonized, organ_fields)
  patients$immunosuppression <- revive_selected_labels(patients_harmonized, medication_fields)

  infections <- events[events$event_type == "Infection", ]
  vaccinations <- events[events$event_type == "Vaccination", ]
  infection_labels <- split(infections$event_label, infections$record_id)
  vaccination_labels <- split(vaccinations$event_label, vaccinations$record_id)
  patients$infection <- lapply(patients$record_id, function(id) {
    sort(unique(as.character(infection_labels[[id]])))
  })
  patients$vaccination <- lapply(patients$record_id, function(id) {
    sort(unique(as.character(vaccination_labels[[id]])))
  })

  # A usable slider also needs limits for an empty cohort or a single age.
  ages <- patients$age[is.finite(patients$age)]
  age_limits <- if (length(ages)) c(floor(min(ages)), ceiling(max(ages))) else c(0, 100)
  age_limits[2] <- max(age_limits[2], age_limits[1] + 1)

  list(
    patients = patients,
    patient_characteristics = patients_harmonized,
    events = events,
    age_limits = age_limits
  )
}

read_revive_dashboard_data <- function(patients_path, data_dir = "data",
                                      output_dir = "data/preprocessed",
                                      skip_files = character()) {
  patients_harmonized <- generate_revive_patients(patients_path, output_dir)
  people <- generate_revive_json(
    data_dir, output_dir,
    record_ids = patients_harmonized$record_id,
    skip_files = skip_files
  )

  # The JSON supplies event histories; patient attributes come from the
  # source CSV through harmonize_patients(). Empty event arrays contribute no rows.
  infection_rows <- lapply(people, function(person) {
    lapply(person$infections, function(report) {
      data.frame(
        record_id = person$record_id,
        event_type = "Infection",
        event_label = as.character(report$infection)[1],
        event_date = as.Date(as.character(report$infection_date)[1], format = "%Y-%m-%d")
      )
    })
  })
  vaccination_rows <- lapply(people, function(person) {
    lapply(person$vaccinations, function(report) {
      data.frame(
        record_id = person$record_id,
        event_type = "Vaccination",
        event_label = as.character(report$vaccine_type)[1],
        event_date = as.Date(as.character(report$vaccine_date)[1], format = "%Y-%m-%d")
      )
    })
  })
  empty_events <- data.frame(
    record_id = character(), event_type = character(),
    event_label = character(), event_date = as.Date(character())
  )
  event_rows <- unlist(c(infection_rows, vaccination_rows), recursive = FALSE)
  events <- do.call(rbind, c(list(empty_events), event_rows))
  events$event_label[is.na(events$event_label) | events$event_label == ""] <- "Type not recorded"

  generate_revive_dashboard_data(patients_harmonized, events)
}
