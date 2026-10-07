# Harmonize the REVIVE patient characteristics export.
# Follow the explicit variable conversions used in u01l_harmonization.R.
# Original columns and REDCap identifiers are retained.
# Dictionary: ref/Revive_dictionary.pdf (October 7, 2026).

harmonize_patients <- function(patients) {
  patients |>
    dplyr::mutate(
      record_id = as.character(record_id),
      age = as.numeric(age_calc_consent),

      # Consent group: dictionary field 24, page 5.
      study_group = dplyr::case_when(
        group_consent_select == 1 ~ "Standard Cohort",
        group_consent_select == 2 ~ "Comprehensive Cohort",
        .default = NA_character_
      ),

      # Demographics: dictionary fields 61, 64, and 66, pages 9-10.
      sex_at_birth = dplyr::case_when(
        sex == 1 ~ "Male",
        sex == 2 ~ "Female",
        sex == 97 ~ "Prefer not to answer",
        .default = NA_character_
      ),
      hispanic_latino = dplyr::case_when(
        ethnicity == 1 ~ "Yes",
        ethnicity == 2 ~ "No",
        ethnicity == 97 ~ "Prefer not to answer",
        .default = NA_character_
      ),
      # Separate indicators preserve every selected racial identity.
      race_white = dplyr::case_when(
        race___1 == 1 ~ "Yes",
        race___1 == 0 ~ "No",
        .default = NA_character_
      ),
      race_black_african_american = dplyr::case_when(
        race___2 == 1 ~ "Yes",
        race___2 == 0 ~ "No",
        .default = NA_character_
      ),
      race_asian = dplyr::case_when(
        race___3 == 1 ~ "Yes",
        race___3 == 0 ~ "No",
        .default = NA_character_
      ),
      race_american_indian_alaska_native = dplyr::case_when(
        race___4 == 1 ~ "Yes",
        race___4 == 0 ~ "No",
        .default = NA_character_
      ),
      race_native_hawaiian_other_pacific_islander = dplyr::case_when(
        race___5 == 1 ~ "Yes",
        race___5 == 0 ~ "No",
        .default = NA_character_
      ),
      race_other = dplyr::case_when(
        race___88 == 1 ~ "Yes",
        race___88 == 0 ~ "No",
        .default = NA_character_
      ),
      race_prefer_not_to_answer = dplyr::case_when(
        race___97 == 1 ~ "Yes",
        race___97 == 0 ~ "No",
        .default = NA_character_
      ),

      # Transplant recipient status: dictionary field 69, pages 10-11.
      transplant_recipient = dplyr::case_when(
        tx_status_yn == 1 ~ "Yes",
        tx_status_yn == 0 ~ "No",
        .default = NA_character_
      ),

      # Most recent transplant organs: dictionary field 73, page 11.
      # Separate indicators preserve multi-organ transplants.
      kidney_transplant = dplyr::case_when(
        recent_tx_organ___1 == 1 ~ "Yes",
        recent_tx_organ___1 == 0 ~ "No",
        .default = NA_character_
      ),
      liver_transplant = dplyr::case_when(
        recent_tx_organ___2 == 1 ~ "Yes",
        recent_tx_organ___2 == 0 ~ "No",
        .default = NA_character_
      ),
      lung_transplant = dplyr::case_when(
        recent_tx_organ___3 == 1 ~ "Yes",
        recent_tx_organ___3 == 0 ~ "No",
        .default = NA_character_
      ),
      heart_transplant = dplyr::case_when(
        recent_tx_organ___4 == 1 ~ "Yes",
        recent_tx_organ___4 == 0 ~ "No",
        .default = NA_character_
      ),
      pancreas_transplant = dplyr::case_when(
        recent_tx_organ___5 == 1 ~ "Yes",
        recent_tx_organ___5 == 0 ~ "No",
        .default = NA_character_
      ),
      intestine_transplant = dplyr::case_when(
        recent_tx_organ___6 == 1 ~ "Yes",
        recent_tx_organ___6 == 0 ~ "No",
        .default = NA_character_
      ),
      islet_cell_transplant = dplyr::case_when(
        recent_tx_organ___7 == 1 ~ "Yes",
        recent_tx_organ___7 == 0 ~ "No",
        .default = NA_character_
      ),
      other_transplant = dplyr::case_when(
        recent_tx_organ___88 == 1 ~ "Yes",
        recent_tx_organ___88 == 0 ~ "No",
        .default = NA_character_
      ),
      other_transplant_description = dplyr::na_if(trimws(recent_tx_organ_oth), ""),

      # Current immunosuppression: dictionary field 92, page 13.
      no_immunosuppression = dplyr::case_when(
        is_meds_elig___0 == 1 ~ "Yes",
        is_meds_elig___0 == 0 ~ "No",
        .default = NA_character_
      ),
      azathioprine = dplyr::case_when(
        is_meds_elig___1 == 1 ~ "Yes",
        is_meds_elig___1 == 0 ~ "No",
        .default = NA_character_
      ),
      belatacept = dplyr::case_when(
        is_meds_elig___2 == 1 ~ "Yes",
        is_meds_elig___2 == 0 ~ "No",
        .default = NA_character_
      ),
      cyclosporine = dplyr::case_when(
        is_meds_elig___3 == 1 ~ "Yes",
        is_meds_elig___3 == 0 ~ "No",
        .default = NA_character_
      ),
      tacrolimus = dplyr::case_when(
        is_meds_elig___4 == 1 ~ "Yes",
        is_meds_elig___4 == 0 ~ "No",
        .default = NA_character_
      ),
      everolimus = dplyr::case_when(
        is_meds_elig___5 == 1 ~ "Yes",
        is_meds_elig___5 == 0 ~ "No",
        .default = NA_character_
      ),
      sirolimus = dplyr::case_when(
        is_meds_elig___6 == 1 ~ "Yes",
        is_meds_elig___6 == 0 ~ "No",
        .default = NA_character_
      ),
      mycophenolate = dplyr::case_when(
        is_meds_elig___7 == 1 ~ "Yes",
        is_meds_elig___7 == 0 ~ "No",
        .default = NA_character_
      ),
      prednisone_methylprednisolone = dplyr::case_when(
        is_meds_elig___8 == 1 ~ "Yes",
        is_meds_elig___8 == 0 ~ "No",
        .default = NA_character_
      ),
      other_immunosuppression = dplyr::case_when(
        is_meds_elig___88 == 1 ~ "Yes",
        is_meds_elig___88 == 0 ~ "No",
        .default = NA_character_
      ),
      other_immunosuppression_description = dplyr::na_if(trimws(is_meds_elig_oth), ""),
      autoimmune_medication = dplyr::na_if(trimws(autoimm_med_name), ""),
      ig_plasma = dplyr::if_else(ig_plasma == 1, "Yes", "NO", NA_character_),
      ig_plasma_3mo = dplyr::if_else(ig_plasma_3mo == 1, "Yes", "NO", NA_character_),
    )
}
