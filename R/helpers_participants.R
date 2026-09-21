revive_participant_dictionary <- function() {
  # ref/RespiratoryViralInfectionVacci.pdf: sex/race p. 19, medications p. 23.
  # The embedded prior_tx_organ field omits its choices on pp. 21-22.
  # Reuse the organ code list printed for organ_tx_28d_type on p. 135;
  # prior_tx_organ branching rules independently confirm 1 = Kidney, 2 = Liver.
  list(
    sex = c(`1` = "Male", `2` = "Female", `97` = "Prefer not to answer"),
    race = c(
      `1` = "White", `2` = "Black or African American", `3` = "Asian",
      `4` = "American Indian or Alaska Native",
      `5` = "Native Hawaiian or Other Pacific Islander",
      `88` = "Other, specify in the text box below:", `97` = "Prefer not to answer"
    ),
    is_meds_elig = c(
      `0` = "None",
      `1` = "Azathioprine",
      `2` = "Belatacept",
      `3` = "Cyclosporine",
      `4` = "Tacrolimus",
      `5` = "Everolimus",
      `6` = "Sirolimus",
      `7` = "MMF",
      `8` = "Steroids",
      `88` = "Other"
    ),
    prior_tx_organ = c(
      `1` = "Kidney", `2` = "Liver", `3` = "Lung", `4` = "Heart",
      `5` = "Pancreas", `6` = "Intestine",
      `7` = "Islet",
      `88` = "Other"
    )
  )
}

read_revive_participants <- function(path) {
  dictionary <- revive_participant_dictionary()
  scalar_fields <- c("record_id", "sex", "age_calc_consent")
  checkbox_fields <- lapply(c("race", "is_meds_elig", "prior_tx_organ"), function(field) {
    paste0(field, "___", names(dictionary[[field]]))
  })
  names(checkbox_fields) <- c("race", "is_meds_elig", "prior_tx_organ")
  # Read only the requested attributes; preserve string IDs and checkbox zeros.
  raw <- readr::read_csv(
    path,
    col_select = tidyselect::all_of(c(scalar_fields, unlist(checkbox_fields, use.names = FALSE))),
    col_types = readr::cols(.default = readr::col_character()),
    na = "", trim_ws = TRUE, name_repair = "minimal",
    show_col_types = FALSE, progress = FALSE
  )
  if (nrow(readr::problems(raw)) > 0L || anyDuplicated(names(raw))) {
    stop("Check malformed rows or duplicate columns in the participant export.")
  }
  for (field in names(raw)) {
    raw[[field]] <- trimws(raw[[field]])
    raw[[field]][raw[[field]] == ""] <- NA_character_
  }
  if (anyNA(raw$record_id)) stop("The participant export contains a missing record_id.")

  coalesce_attribute <- function(values, field) {
    values <- unique(values[!is.na(values)])
    if (length(values) > 1L) {
      stop("Conflicting nonmissing values for ", field,
           " within a participant; reconcile the export before joining.")
    }
    if (length(values)) values[[1]] else NA_character_
  }
  label_checkboxes <- function(person, field) {
    values <- vapply(checkbox_fields[[field]], function(column) {
      coalesce_attribute(person[[column]], column)
    }, character(1))
    if (any(!is.na(values) & !values %in% c("0", "1"))) {
      stop("Checkbox values must be 0, 1, or blank for ", field, ".")
    }
    if (all(is.na(values))) return(NULL)
    # Lists retain JSON arrays even when exactly one option was selected.
    unname(as.list(dictionary[[field]][which(values == "1")]))
  }

  groups <- split(seq_len(nrow(raw)), raw$record_id)
  unname(lapply(groups, function(rows) {
    person <- raw[rows, , drop = FALSE]
    sex <- coalesce_attribute(person$sex, "sex")
    if (!is.na(sex) && !sex %in% names(dictionary$sex)) {
      stop("Unknown sex code in the participant export; update the dictionary.")
    }
    age_raw <- coalesce_attribute(person$age_calc_consent, "age_calc_consent")
    age <- suppressWarnings(as.numeric(age_raw))
    if (!is.na(age_raw) && (is.na(age) || !is.finite(age) || age < 0)) {
      stop("age_calc_consent must be a nonnegative number or blank.")
    }
    list(
      record_id = person$record_id[[1]],
      Sex = unname(dictionary$sex[sex]),
      Age = age,
      Race = label_checkboxes(person, "race"),
      Immunosuppression = label_checkboxes(person, "is_meds_elig"),
      SOTRs = label_checkboxes(person, "prior_tx_organ")
    )
  }))
}

join_revive_participants <- function(histories, participant_data) {
  history_ids <- vapply(histories, `[[`, character(1), "record_id")
  participant_ids <- vapply(participant_data, `[[`, character(1), "record_id")
  if (anyNA(c(history_ids, participant_ids)) ||
      any(!nzchar(c(history_ids, participant_ids))) ||
      anyDuplicated(history_ids) || anyDuplicated(participant_ids)) {
    stop("Each input to the participant join must have one nonmissing record_id per participant.")
  }
  ids <- sort(union(history_ids, participant_ids))
  history_index <- match(ids, history_ids)
  participant_index <- match(ids, participant_ids)
  lapply(seq_along(ids), function(i) {
    person <- if (is.na(history_index[i])) {
      list(record_id = ids[i], infections = list(), vaccinations = list())
    } else {
      histories[[history_index[i]]]
    }
    attributes <- if (is.na(participant_index[i])) {
      list(Sex = NA_character_, Age = NA_real_, Race = NULL,
           Immunosuppression = NULL, SOTRs = NULL)
    } else {
      participant_data[[participant_index[i]]]
    }
    # Single-bracket assignment preserves explicit NULL fields in the JSON.
    person[names(attributes)] <- attributes
    person
  })
}
