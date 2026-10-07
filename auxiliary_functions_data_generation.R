# load necessary files
# functions to load HOPE data


generate_liver_data <- function() {
  print("Load necessary files")

  source("//win.ad.jhu.edu/Users$/HOME/Data_Harmonization/load_hope_data.R")
  # load functions to harmonize data sets
  source("//win.ad.jhu.edu/Users$/HOME/Data_Harmonization/harmonize_data.R")

  donors <- load_donors_data()
  hope_liver <- load_hope_data("u01l")

  # load auxiliary files
  # load malignancy data
  malignacy_cases <- read_csv("data/Malignancies.csv")

  # load rejection biopsy results
  # rejection_cases <- read_csv("../data/2026_02_10_Liver_rejection.csv")

  # load centers code
  centers_cod <- read_csv("data/Centers_and_OPO_codes.csv")

  # load primary indication for transplant
  primary_indication <- read_csv("data/Adjudicated_primary_indication_for_LT.csv")

  # load centers names
  centers_name <- read_csv("data/Centers_names.csv")

  # read key pairs data
  key_pairs <- read_csv("data/Recipients_ids.csv")

  # harmonize data set
  # hope_liver <- harmonize_data(hope_liver, dataset = "u01k")

  print("Starting data harmonization process")

  hope_liver <- hope_liver |>
    mutate(group = case_when(
      organ_hiv_status == 1 ~ "HIV D+/R+",
      organ_hiv_status %in% c(2, 3) ~ "HIV D-/R+"
    ))
  
  # get a list of transplanted patients
  tranplanted_ids <- hope_liver |>
    filter(!is.na(dot)) |>
    pull(studyid) |>
    unique()

  # filter hope_liver to only transplanted patients
  hope_liver <- hope_liver |>
    filter(studyid %in% tranplanted_ids)
  

  # get infection data
  infection_data <- infection_long(
    data = hope_liver,
    events = "infection[0-9]{1,2}_date",
    origin = dot
  )

  # check if infection occured before or after transplant
  infection_data <- infection_data |>
    mutate(
      after_trans =
        case_when(
          infection_date >= dot ~ "Yes",
          TRUE ~ "No"
        )
    ) |> 
    filter(studyid %in% tranplanted_ids)
  
  # create a data for sae events
  sae_data <-  hope_liver |> 
    select(studyid, dot, group,sae_date, starts_with("sae_criteria")) |> 
    filter(studyid %in% tranplanted_ids) |> 
    group_by(studyid) |>
    mutate(group = max(group, na.rm = TRUE),
          dot = max(dot, na.rm = TRUE)
    ) |>
    ungroup() |> 
    filter(if_any(starts_with("sae_criteria"), ~ . == 1)) |> 
    pivot_longer(cols = starts_with("sae_criteria"), 
                 names_to = "sae_criteria", 
                 values_to = "sae_criteria_value") |> 
    filter(sae_criteria_value == 1)

  
  # ordered sae rank
  sae_data <- sae_data |> 
    mutate(sae_criteria = case_when(
      sae_criteria == "sae_criteria___1" ~ "Death",
      sae_criteria == "sae_criteria___2" ~ "Life-threatening",
      sae_criteria == "sae_criteria___3" ~ "Hospitalization",
      .default = "Other"
    ))
  
  # add graf loss to sae
  graft_loss_sae <- hope_liver |> 
    select(studyid, dot, group, graft_loss_date_l) |>
    group_by(studyid) |>
    summarise(
      dot = max(dot, na.rm = TRUE),
      group = max(group, na.rm = TRUE),
      graft_loss_date_l = max(graft_loss_date_l, na.rm = TRUE)
    ) |> 
    filter(!is.infinite(graft_loss_date_l)) |> 
    rename(sae_date = graft_loss_date_l) |>
    mutate(sae_criteria_value = 1,
           sae_criteria = "Graft loss"
    )
  sae_data <- bind_rows(sae_data, graft_loss_sae)
  
  
  #defining a list with time points for longitudianal analysis
  #time points will be get at each 6 monhts (26 weeks)
  time_points <- c("week_26_arm_1","week_52_year_1_arm_1", "week_78_arm_1",
                   "week_104_year_2_arm_1","week_130_arm_1", 
                   "week_156_year_3_arm_1", "week_182_arm_1", "week_208_year_4_arm_1")
  
  # create longitudinal data for fib4 index at each time point
  fib4_index <- hope_liver |> 
    filter(redcap_event_name %in% time_points) |> 
    select(studyid, redcap_event_name, ast_fu, ast_unit_fu, alt_fu, alt_unit_fu, bilirubin_fu, 
           platelet_fu, platelet_unit_fu, creatinine_fu, creatinine_unit_fu) |> 
    group_by(studyid,redcap_event_name) |>
    summarise_all(max, na.rm = TRUE)

  # define race based on race and ethinicity in black white non hispanic, white hispanic and others
  hope_liver <- hope_liver |>
    mutate(
      race_combined =
      # check about the definition of race and ethnicity
        case_when(
          race___2 == 1 ~ "Black or African American",
          race___1 == 1 ~ "White",
          race___3 == 1 ~ "Other",
          race___4 == 1 ~ "Other",
          race___5 == 1 ~ "Other",
          race___0 == 1 ~ "Other",
          .default = NA_character_
        ),
      ethnicity = if_else(ethnicity == 1, "Hispanic or Latino",
        "Not Hispanic", NA_character_
      ),
      gender = if_else(gender == 0, "male", "female", NA_character_)
    )

  # recode hbv_core_ab and hcv_ab
  hope_liver <- hope_liver |>
    mutate(
      # Hepatitis B core ab positive
      hep_b_antibody =
        case_when(
          hbv_core_ab == 1 ~ "Positive",
          hbv_core_ab == 0 ~ "Negative",
          .default =  NA
        ),

      # hepatitis B surface antigen
      hep_b_surface_ag =
        case_when(
          hbv_sag == 1 ~ "Positive",
          hbv_sag == 0 ~ "Negative",
          .default = NA
        ),
      hep_b_surface_ab =
        case_when(
          hbv_sab == 1 ~ "Positive",
          hbv_sab == 0 ~ "Negative",
          .default = NA
        ),

      # Hepatitis B IGG or total core antibody
      hep_b_igg_core =
        case_when(
          hbv_core_ab == 1 ~ "Positive",
          hbv_core_ab == 0 ~ "Negative",
          .default = NA
        ),

      # hepatitis C antibody
      hep_c_antibody =
        case_when(
          hcv_ab == 1 ~ "Positive",
          hcv_ab == 0 ~ "Negative",
          .default = NA
        ),

      # hepatitis C PCR
      hep_c_pcr =
        case_when(
          hcv_pcr == 1 ~ "Positive",
          hcv_pcr == 0 ~ "Negative",
          .default = NA
        ),

      # CMV IgG Ab
      cmv_igg_pos =
        case_when(
          cmv == 1 ~ "Positive",
          cmv == 0 ~ "Negative",
          .default = NA
        ),

      # donor CMV status
      donor_cmv_result =
        case_when(
          donorsero_cmvab == 1 ~ "Positive",
          donorsero_cmvab == 0 ~ "Negative",
          .default = NA
        ),

      # blood_type
      blood_type = case_when(
        blood_type == 1 ~ "A",
        blood_type == 2 ~ "B",
        blood_type == 3 ~ "AB",
        blood_type == 4 ~ "O",
      ),
      # is patient a slk
      is_slk = case_when(
        donor_organs == 1 ~ "No",
        donor_organs == 2 ~ "Yes"
      ),
      # define diabetes
      diabetes = case_when(
        endocrine___1 == 1 ~ "Yes",
        endocrine___1 == 0 ~ "No",
        .default = "No"
      ),
      # define hypertension
      hypertension = case_when(
        renal___8 == 1 | cardiovascular___7 == 1 ~ "Yes",
        renal___8 == 0 & cardiovascular___7 == 0 ~ "No",
        .default = "No"
      )
    )


  # get viral load at day 0
  hope_liver <- hope_liver |>
    mutate(
      viral_load_transp =
        case_when(
          str_detect(str_to_lower(hiv_pcr), "detected") ~ 0,
          str_detect(str_to_lower(hiv_pcr), "undet") ~ 0,
          str_detect(str_to_lower(hiv_pcr), "non") ~ 0,
          .default = as.numeric(
            str_extract_all(
              hiv_pcr,
              "\\d+\\.?\\d*",
              simplify = TRUE
            )
          )
        ),
      viral_load_fu =
        case_when(
          str_detect(str_to_lower(hiv_pcr), "detected") ~ 0,
          str_detect(str_to_lower(hiv_pcr), "undet") ~ 0,
          str_detect(str_to_lower(hiv_pcr), "non") ~ 0,
          .default = as.numeric(
            str_extract_all(
              hiv_pcr,
              "\\d+\\.?\\d*",
              simplify = TRUE
            )
          )
        )
    )

  # get years with HIV
  hope_liver <- hope_liver |>
    mutate(
      year_hiv_diag = get_year(first_positive_ab_y),
      year_dot = year(dot)
    )


  # create a list of others cause of HIV
  other_cause_of_hiv <- c(
    "hiv_risk_category___6", "hiv_risk_category___3",
    "hiv_risk_category___5", "hiv_risk_category___0"
  )

  # define hiv acquisition risk
  hope_liver <- hope_liver |>
    mutate(
      msm = if_else(hiv_risk_category___1 == 1, "Yes", "No"),
      hetero_sex_risk = if_else(hiv_risk_category___2 == 1, "Yes", "No"),
      idu = if_else(hiv_risk_category___4 == 1, "Yes", "No"),
      other = case_when(
        any(across(other_cause_of_hiv) == 1) ~ "Yes",
        .default = "No"
      )
    )


  # create a binary varible for viral load at transplant < 200 ml
  hope_liver <- hope_liver |>
    mutate(viral_load_transp_binary = if_else(viral_load_transp < 200, "<= 200", ">200"))

  # check most common art classes
  hope_liver |>
    select(starts_with("art_class___")) |>
    summarise(across(everything(), ~ sum(., na.rm = TRUE))) |>
    t()


  # check if recipient get insti or nrrti treatment
  hope_liver <- hope_liver |>
    mutate(
      INSTI =
        case_when(
          art_class___1 == 1 ~ "Yes",
          .default = "No"
        ),
      NNRTI =
        case_when(
          art_class___2 == 1 ~ "Yes",
          .default = "No"
        ),
      pi_cobicistat =
        case_when(
          art_class___3 == 1 ~ "Yes",
          .default = "No"
        ),
    )

  # define primary indication for transplant
  hope_liver <- hope_liver |>
    left_join(primary_indication |> select(studyid, Adjudicated_primary_ind_for_Table_1),
      by = "studyid"
    ) |>
    rename(primary_ind_for_transp = Adjudicated_primary_ind_for_Table_1)


  # recode imunossupression induction
  hope_liver <- hope_liver |>
    mutate(
      steroids =
        case_when(
          induction_med___1 == 1 ~ "Yes",
          .default = "No"
        ),
      basiliximab =
        case_when(
          induction_med___3 == 1 ~ "Yes",
          .default = "No"
        )
    )

  # limphocyte depleting therapy are ATG, ATGAM and alemtuzumab
  hope_liver <- hope_liver |>
    rowwise() |>
    mutate(
      lymphocite_dep_the =
        case_when(
          induction_med___2 == 1 | induction_med___4 == 1 | induction_med___5 == 1 |
            induction_med___7 == 1 ~ "Yes",
          .default = "No"
        )
    )

  # check patients that received lymphocyte depleting therapy
  hope_liver |>
    group_by(studyid) |>
    summarise(
      lymphocite_dep_the = max(lymphocite_dep_the, na.rm = TRUE),
      group = max(group, na.rm = TRUE)
    ) |>
    group_by(group, lymphocite_dep_the) |>
    summarise(n = n())


  # recode maintenance of immunosuppression
  hope_liver <- hope_liver |>
    mutate(
      predinisone =
        case_when(
          is_central___5 == 1 ~ "Yes",
          .default = "No"
        ),
      mycophenolate =
        case_when(
          is_central___4 == 1 ~ "Yes",
          .default = "No"
        ),
      tacrolimus =
        case_when(
          is_central___1 == 1 ~ "Yes",
          .default = "No"
        ),
      cd4_percent = as.numeric(str_extract_all(cd4_percent, "\\d+\\.?\\d*"))
    )


  # OI tab
  list_of_oi <- c(
    "bl_aspergillosis", "bl_bartonella", "bl_candi_lungs", "bl_cmv", "bl_candi_esophagus",
    "bl_cervical_cancer", "bl_coccidioidomycosis", "bl_cryptococcosis",
    "bl_cryptosporidiosis", "bl_cryptosporidiosis", "bl_encephalopathy",
    "bl_herpes_simplex", "bl_hhv6", "bl_histoplasmosis", "bl_isosporiasis",
    "bl_kaposi", "bl_lymphoid", "bl_lymphoid", "bl_mucormycosis",
    "bl_mucormycosis", "bl_mucormycosis", "bl_myco_tuberculosis", "bl_myco_other",
    "bl_nocardia", "bl_pjp", "bl_pml", "bl_salmonella",
    "bl_toxoplasmosis", "bl_varicella_zoster", "bl_wasting_syndrome"
  )


  # check for history of previous OI
  hope_liver <- hope_liver |>
    rowwise() |>
    mutate(
      history_Oi =
        case_when(
          any(c_across(list_of_oi)) == 1 ~ "Yes",
          .default = "No"
        )
    )


  # calculating allcoation MELD and Biologic MELD
  hope_liver <- hope_liver |>
    # Create allo_meld_day0 based on meld_exception_yn_d0 & missing values
    mutate(
      allo_meld = ifelse(meld_exception_yn_d0 == 0 & is.na(meld_exception_d0),
        meld_day0, NA
      ),

      # Copy meld_exception_d0 to meld_excp_d0
      meld_excp_d0 = meld_exception_d0,

      # Replace meld_excp_d0 with empty string if it matches "Status 1a"
      meld_excp_d0 = ifelse(str_detect(meld_excp_d0, "Status 1a"), "",
        meld_excp_d0
      ),

      # Convert meld_excp_d0 to numeric, replacing non-numeric values with NA
      meld_excp_d0 = as.numeric(meld_excp_d0)
    ) |>
    # Replace allo_meld_day0 based on meld_exception_yn_d0 & non-missing meld_exception_d0
    mutate(
      allo_meld = ifelse(meld_exception_yn_d0 == 1 & !is.na(meld_exception_d0),
        meld_excp_d0, allo_meld
      )
    )
  hope_liver <- hope_liver |>
    mutate(
      allo_meld = 
        ifelse(studyid == "L201-849R", 28, allo_meld)
    )

  hope_liver <- hope_liver |>
    mutate(
      primary_ind_for_transp =
        case_when(
          primary_ind_for_transp == "Acute hepatic necrosis due to drug" ~ "Other",
          primary_ind_for_transp == "Cryptogenic cirrhosis" ~ "Other",
          primary_ind_for_transp == "Autoimmune" ~ "Autoimmune/PBC/PSC",
          primary_ind_for_transp == "Primary biliary" ~ "Autoimmune/PBC/PSC",
          primary_ind_for_transp == "Primary biliary cirrhosis" ~ "Autoimmune/PBC/PSC",
          primary_ind_for_transp == "Primary sclerosing cholangitis" ~ "Autoimmune/PBC/PSC",
          primary_ind_for_transp == "Alcoholic liver disease" ~ "Alcoholic liver disease",
          primary_ind_for_transp == "Acute alcoholic hepatitis" ~ "Alcoholic liver disease",
          .default = primary_ind_for_transp
        )
    )
  hope_liver <- hope_liver |>
    group_by(studyid) |>
    select(
      # socio demographics
      group, tx_age, gender, race_combined, ethnicity, blood_type, hypertension, diabetes,
      primary_ind_for_transp, meld_day0, allo_meld, is_slk,


      # hiv and co-infections
      year_hiv_diag, year_dot, history_Oi, msm, hetero_sex_risk, idu, other, cd4_absolute,
      cd4_percent, viral_load_transp_binary, hep_b_igg_core, hep_b_surface_ag,
      hep_b_surface_ab, hep_c_antibody, hep_c_pcr, cmv_igg_pos, donor_cmv_result, INSTI,
      pi_cobicistat, cold_ischemia_liver,

      # immunossupression induction
      steroids, basiliximab, predinisone, mycophenolate, tacrolimus, lymphocite_dep_the,
      donor_unos_id_liver,

      # dates
      date_death, graft_loss_date_l, dot, sae_date, infection1_date, rejection_date_liver,
      rejection_biopsy_liver, date_last_visit, date_remove, date_visit, date_test1,

      # sae
      starts_with("sae_criteria"), infection1_aids, viral_breakthrough
    ) |>
    summarise(
      max_sae_date = max(na.omit(sae_date)),
      date_visit = max(na.omit(date_visit)),
      across(where(is.numeric), \(x) max(na.omit(x))),
      across(where(is.character), \(x) max(na.omit(x))),
      across(where(is.Date), \(x) min(na.omit(x)))
    ) |>
    mutate(
      years_with_hiv = year_dot - year_hiv_diag,
      cmv_mismatch = case_when(
        cmv_igg_pos == "Negative" & donor_cmv_result == "Positive" ~ "D+/R-",
        cmv_igg_pos == "Negative" & donor_cmv_result == "Negative" ~ "D-/R-",
        .default = "No"
      ),
      primary_ind_for_transp = if_else(
        is.na(primary_ind_for_transp), "Other", primary_ind_for_transp
      )
    )

  # change all infinite values to NA
  hope_liver <- hope_liver |>
    mutate(
      across(where(is.numeric), ~ if_else(is.infinite(.x), NA_real_, .x)),
      across(where(is.Date), ~ if_else(is.infinite(.x), NA_Date_, .x))
    )


  # select donors that were used in the transplantation
  donors_ids <- hope_liver |>
    filter(studyid %in% tranplanted_ids) |>
    pull(donor_unos_id_liver) |>
    na.omit()

  # filter donors for those used in transplantation
  donors <- donors |>
    filter(unos_donor_id %in% donors_ids)


  # classify donors 2 is FP
  donors <- donors |>
    mutate(
      group =
        case_when(
          donor_hiv_status == 1 ~ "HIV D+",
          donor_hiv_status %in% c(2, 3, 4) ~ "HIV D-"
        ),
      age = as.numeric(age)
    )

  # recode race
  donors <- donors |>
    mutate(
      race_combined =
        case_when(
          race___0 == 1 ~ "White",
          race___1 == 1 ~ "Black or African American",
          race___2 == 1 ~ "Other",
          race___3 == 1 ~ "Other",
          race___4 == 1 ~ "Other",
          race___5 == 1 ~ "Other",
          .default = NA_character_
        ),
      ethnicity = case_when(
        hisp == 1 ~ "Hispanic or Latino",
        hisp == 0 ~ "Not Hispanic",
        .default = NA_character_
      )
    )

  # recode hepb
  donors <- donors |>
    mutate(
      hep_b_antibody =
        case_when(
          anti_hbcab == 1 ~ "Positive",
          anti_hbcab == 0 ~ "Negative",
          .default = NA
        )
    )
  # recode hep c
  donors <- donors |>
    mutate(
      hep_c_antibody =
        case_when(
          anti_hcv == 1 ~ "Positive",
          anti_hcv == 0 ~ "Negative",
          .default = NA
        ),
      cmv =
        case_when(
          anti_cmv == 1 ~ "Positive",
          anti_cmv == 0 ~ "Negative",
          .default = NA
        ),
      hbv_nat = case_when(
        hbv_nat == 1 ~ "Positive",
        hbv_nat == 0 ~ "Negative",
        .default = NA
      )
    )

  # get hep c nat
  donors <- donors |>
    mutate(
      hep_c_nat =
        case_when(
          hcv_nat == 1 ~ "Positive",
          hcv_nat == 0 ~ "Negative",
          .default = NA
        )
    )


  # code donors for FP tests
  donors <- donors |>
    mutate(
      fp =
        case_when(
          donor_hiv_status == 2 ~ "Yes",
          .default = "No",
        ),
      false_positive_category =
        case_when(
          false_positive_category == 1 ~ "Ab Positive",
          false_positive_category == 2 ~ "NAT Positive",
          .default = NA,
        ),
      diabetes =
        case_when(
          diabetes_hist == 0 ~ "Yes",
          diabetes_hist == 1 ~ "No",
          .default = NA,
        ),
      cancer = case_when(
        cancer_hist == 0 ~ "Yes",
        cancer_hist == 1 ~ "No",
      ),
      hypertension =
        case_when(
          htn_hist == 0 ~ "Yes",
          htn_hist == 1 ~ "No",
          .default = NA,
        ),
      dcd =
        case_when(
          dcd_yn == 0 ~ "Yes",
          dcd_yn == 1 ~ "No",
          .default = NA
        ),
      new_hiv_diag =
        case_when(
          hiv_status_discovery_time == 1 ~ "Yes",
          hiv_status_discovery_time == 0 ~ "No",
          .default = NA,
        ),
      on_art =
        case_when(
          art_yn == 1 ~ "Yes",
          art_yn == 0 ~ "No",
          .default = NA,
        ),
      hiv_ab_test =
        case_when(
          anti_hiv_i == 1 ~ "Positive",
          anti_hiv_i == 0 ~ "Negative",
          .default = NA,
        ),
      hiv_nat_test =
        case_when(
          hiv_nat == 1 ~ "Positive",
          hiv_nat == 0 ~ "Negative",
          .default = NA,
        ),
      # blood_type
      blood_type = case_when(
        blood_type == 0 ~ "A",
        blood_type == 1 ~ "B",
        blood_type == 2 ~ "AB",
        blood_type == 3 ~ "O",
      ),
    )


  # create a list of others cause of HIV
  other_cause_of_hiv_donors <- c(
    "hiv_transmission___5", "hiv_transmission___7",
    "hiv_transmission___8", "hiv_transmission___9"
  )

  # define hiv donors risk
  donors <- donors |>
    rowwise() |>
    mutate(
      msm = if_else(hiv_transmission___1 == 1, "Yes", "No"),
      hetero_sex_risk = if_else(hiv_transmission___10 == 1, "Yes", "No"),
      idu = if_else(hiv_transmission___4 == 1, "Yes", "No"),
      other = case_when(
        any(across(other_cause_of_hiv_donors) == 1) ~ "Yes",
        .default = "No"
      )
    )

  # standardize ND in hiv_pcr to be <40
  donors <- donors |>
    mutate(
      hiv_pcr =
        case_when(
          str_detect(str_to_lower(hiv_pcr), "nd") ~ "<40",
          .default = hiv_pcr
        )
    )

  # get hwo many patients on ART had less than 400 copies of virus per mL
  donors <- donors |>
    mutate(
      hiv_pcr_parsed = as.numeric(parse_number(hiv_pcr)),
      cd4_percentage = as.numeric(parse_number(cd4_percentage))
    ) |>
    mutate(
      vl_below_1000 =
        case_when(
          on_art == "Yes" & hiv_pcr_parsed < 1000 ~ "Yes",
          on_art == "Yes" & hiv_pcr_parsed >= 1000 ~ "No",
          .default = NA
        ),
      vl_below_400 =
        case_when(
          on_art == "Yes" & hiv_pcr_parsed < 400 ~ "Yes",
          on_art == "Yes" & hiv_pcr_parsed >= 400 ~ "No",
          .default = NA
        )
    )


  # ot those of art get viral load
  donors <- donors |>
    mutate(
      vl_off_art =
        case_when(
          on_art == "No" ~ hiv_pcr_parsed,
          .default = NA
        )
    )
  donors_dri <- donors |>
    mutate(
      age = as.numeric(age),
      cod =
        case_when(
          death_cause == 0 ~ "anoxia",
          death_cause == 1 ~ "cva",
          death_cause == 2 ~ "trauma",
          .default = "other"
        ),
      eth =
        case_when(
          race_combined == "Black" ~ "black",
          race_combined %in% c("White Non hispanic", "White Hispanic") ~ "white",
          .default = "other"
        ),
      dcd =
        case_when(
          dcd_yn == 0 ~ 1,
          .default = 0
        ),
      height = as.numeric(height_cm),
      split = 0
    )

  # get cold schemia time from hope_liver
  organ_cit <- hope_liver |>
    filter(studyid %in% tranplanted_ids) |>
    group_by(studyid) |>
    summarise(
      cold_ischemia_liver = max(cold_ischemia_liver, na.rm = TRUE),
      donor_unos_id_liver = max(donor_unos_id_liver, na.rm = TRUE)
    )

  # add CIT to donors_dri
  donors_dri <- donors_dri |>
    left_join(organ_cit |> select(donor_unos_id_liver, cold_ischemia_liver),
      by = c("unos_donor_id" = "donor_unos_id_liver")
    ) |>
    rename(cit = cold_ischemia_liver)


  # add center name to donors_dri
  donors_dri <- donors_dri |>
    left_join(centers_cod |> select(center_location, center_cod),
      by = c("tx_ctr_liver" = "center_cod")
    )

  # add opo location
  donors_dri <- donors_dri |>
    left_join(centers_cod |> select(opo_cod, opo_city_location),
      by = c("opo" = "opo_cod")
    )

  # split center location into city and sate using ,
  donors_dri <- donors_dri |>
    mutate(
      center_city = str_trim(str_split_fixed(center_location, ",", 2)[, 1]),
      center_state = str_trim(str_split_fixed(center_location, ",", 2)[, 2]),
      opo_city = str_trim(str_split_fixed(opo_city_location, ",", 2)[, 1]),
      opo_state = str_trim(str_split_fixed(opo_city_location, ",", 2)[, 2])
    )


  # getiing sahre for each donor
  donors_dri <- donors_dri |>
    mutate(
      share = case_when(
        center_city == opo_city ~ "local",
        center_city != opo_city & center_state == opo_state ~ "regional",
        center_state != opo_state ~ "national",
        .default = "national"
      )
    )

  # apply function to compute LDRI
  donors_dri <- donors_dri |>
    rowwise() |>
    mutate(
      ldri = liver_dri_feng2006(
        age = age, cod = cod, eth = eth,
        dcd = dcd, split = split, share = share,
        cit = cit / 60, height = height
      )
    )

  # add ldri to donors
  donors <- donors |>
    left_join(donors_dri |> select(studyid, ldri, cod), by = "studyid")

  # create outcomes data
  outcomes_data <- hope_liver |>
    select(
      studyid, group, date_death, graft_loss_date_l, dot, sae_date,
      infection1_date, starts_with("sae_criteria"), rejection_date_liver,
      rejection_biopsy_liver, viral_breakthrough, date_last_visit, date_remove,
      date_visit, infection1_aids,donor_unos_id_liver, allo_meld, cold_ischemia_liver,
      hep_c_antibody, lymphocite_dep_the, date_test1, blood_type
    ) |> 
    rename(blood_type_recipient = blood_type)

  # adjucate infection according to Tao's list
  outcomes_data <- outcomes_data |>
    mutate(
      infection1_date =
        case_when(
          studyid == "L201-655R" ~ ymd("2019-12-02"),
          studyid == "L201-716R" ~ ymd("2023-11-14"),
          studyid == "L201-991R" ~ ymd("2023-12-26"),
          studyid == "L201-407R" ~ ymd("2024-03-01"),
          studyid == "L201-738R" ~ ymd("2021-03-02"),
          studyid == "L212-021R" ~ ymd("2022-12-19"),
          .default = NA
        )
    )

  # fix time contr for recipient L212-021
  outcomes_data <- outcomes_data |>
    rowwise() |>
    mutate(date_test1 = if_else(studyid == "L212-021R", ymd("2021-12-21"), date_test1))


  # fix date last visit and add censoring date as 01 december 2025
  outcomes_data <- outcomes_data |>
    rowwise() |>
    mutate(
      date_last_visit = if_else(is.na(date_last_visit), date_visit, date_last_visit),
      censoring_date = ymd("2025-12-01")
    )

  # fix sae date for patient L201-407R
  outcomes_data <- outcomes_data |>
    rowwise() |>
    mutate(sae_date = if_else(studyid == "L201-407R", ymd("2024-03-01"), sae_date))


  # define date list
  date_list <- c(
    "date_death", "graft_loss_date_l", "sae_date", "infection1_date",
    "date_test1"
  )

  date_list_sensitivity <- c(
    "date_death", "graft_loss_date_l", "infection1_date",
    "date_test1"
  )

  date_to_fix <- c("sae_date", "infection1_date", "max_sae_date")

  end_date <- c("date_last_visit", "date_remove", "date_visit", "censoring_date")

  # fix last visit date of patient L223-817R
  outcomes_data <- outcomes_data |>
    rowwise() |>
    mutate(
      date_last_visit = if_else(studyid == "L223-817R", ymd("2023-12-22"), date_last_visit),
      date_visit = if_else(studyid == "L223-817R", ymd("2023-12-22"), date_visit)
    )


  # create a composite outcome
  outcomes_data <- outcomes_data |>
    rowwise() |>
    mutate(
      composite =
        case_when(
          any(across(date_list) >= 1) ~ 1,
          .default = 0
        ),
      composite_sentitivity =
        case_when(
          any(across(date_list_sensitivity) >= 1) ~ 1,
          .default = 0
        )
    )


  # get date of composite
  outcomes_data <- outcomes_data |>
    mutate(
      time_to_composite =
        case_when(
          composite == 1 ~ min(c_across(date_list), na.rm = TRUE),
          .default = NA
        ),
      time_to_composite_sentitivity =
        case_when(
          composite_sentitivity == 1 ~ min(c_across(date_list_sensitivity), na.rm = TRUE),
          .default = NA
        )
    )


  # compute time contribution in months
  outcomes_data <- outcomes_data |>
    rowwise() |>
    mutate(
      time_contr_comp =
        case_when(
          composite == 1 ~
            as.numeric(time_to_composite - dot,
              units = "days"
            ),
          composite == 0 ~
            as.numeric(
              min(
                c_across(end_date),
                na.rm = TRUE
              ) - dot,
              units = "days"
            )
        ),
      time_contr_comp_sens =
        case_when(
          composite_sentitivity == 1 ~
            as.numeric(time_to_composite_sentitivity - dot,
              units = "days"
            ),
          composite_sentitivity == 0 ~
            as.numeric(
              min(
                c_across(end_date),
                na.rm = TRUE
              ) - dot,
              units = "days"
            )
        )
    ) |>
    mutate(
      time_contr_comp = as.integer((time_contr_comp / 31.44)),
      time_contr_comp_sens = as.integer((time_contr_comp_sens / 31.44)),
      group = factor(group, levels = c("HIV D-/R+", "HIV D+/R+"))
    )

  # give .5 if time_contr is equal 0
  outcomes_data <- outcomes_data |>
    mutate(
      time_contr_comp = if_else(as.integer(time_contr_comp) == 0, 0.5,
        time_contr_comp
      ),
      time_contr_comp_sens = if_else(as.integer(time_contr_comp_sens) == 0, 0.5,
        time_contr_comp_sens
      )
    )
  # compute median time contribution
  outcomes_data <- outcomes_data |>
    mutate(time_contr = as.numeric(date_last_visit - dot), units = "days")

  # convert time contribution to years
  outcomes_data <- outcomes_data |>
    mutate(time_contr = as.integer(time_contr / 365.25))

  # add donors information to outcomes data
  outcomes_data <- outcomes_data |>
    left_join(donors |> rename(donor_unos_id_liver = unos_donor_id) |> 
      select(donor_unos_id_liver, sex, race_combined, ethnicity,
             blood_type, cmv, cod)) |> 
    rename(blood_type_donor = blood_type, race_combined_donor = race_combined, 
           ethnicity_donor = ethnicity, donor_sex = sex, donor_cmv = cmv)
  
  # define time contribution for death 
  outcomes_data <- outcomes_data |> 
  
  #define death event
  mutate(
    death = 
      case_when(
        !is.na(date_death) ~ 1,
        .default = 0
      )
  ) |>
  
  #compute time contribution for death
  mutate(
    time_contr_death =
               case_when(
                   death == 1 ~ as.numeric(date_death - dot, 
                           units = "days"),
                   death == 0 ~ as.numeric(min(c_across(end_date), 
                         na.rm = TRUE) - dot,
                       units = "days"))) |>
  #change death to month
  mutate(time_contr_death = as.integer((time_contr_death/31.44))) |>
  # add 0.5 for time = 0
  mutate(time_contr_death = 
           if_else(as.integer(time_contr_death) == 0, 0.5,time_contr_death))
  
  outcomes_data <- outcomes_data |> 
  
    #define death or hiv breaktrough
    mutate(death_or_hivB = 
            case_when(
              !is.na(date_death) | viral_breakthrough == 1 ~ 1,
              .default = 0))

  #define time contribution to death of HIV breakthrough
  outcomes_data <- outcomes_data |> 
    mutate(
      time_contr_death_or_hivB =
        case_when(
          death_or_hivB == 1 ~ as.numeric(
            min(c_across(c("date_death", "date_test1")), na.rm = TRUE) - dot,
            units = "days"),
        death_or_hivB == 0 ~ as.numeric(
          min(c_across(end_date), na.rm = TRUE) - dot,
          units = "days")
      )
  ) |>
  mutate(time_contr_death_or_hivB = as.integer((time_contr_death_or_hivB/31.44))) |>
  mutate(time_contr_death_or_hivB = 
           if_else(as.integer(time_contr_death_or_hivB) == 0, 0.5,
                   time_contr_death_or_hivB))
  
  
  #define death or hiv breaktrough or graft loss
  outcomes_data <- outcomes_data |> 
    mutate(death_hivB_graftloss = 
           case_when(
             !is.na(date_death) | viral_breakthrough == 1 | !is.na(graft_loss_date_l) ~ 1,
             .default = 0
           ))

  #define time contribution
  outcomes_data <- outcomes_data |> 
    mutate(
      time_contr_death_hivB_graftloss =
        case_when(
          death_hivB_graftloss == 1 ~ as.numeric(
            min(
              c_across(c("date_death", "date_test1", "graft_loss_date_l")), 
              na.rm = TRUE) - dot,
            units = "days"),
          death_hivB_graftloss == 0 ~ as.numeric(
            min(c_across(end_date), na.rm = TRUE) - dot,units = "days"))) |> 
    mutate(time_contr_death_hivB_graftloss = 
      as.integer((time_contr_death_hivB_graftloss/31.44))) |> 
    mutate(time_contr_death_hivB_graftloss = 
            if_else(as.integer(time_contr_death_hivB_graftloss) == 0, 0.5,
                    time_contr_death_hivB_graftloss))
  
  # define death, hiv breakthrough, graft loss, or SAE  
  outcomes_data <- outcomes_data |> 
    mutate(death_hivB_graftloss_sae = 
            case_when(
              !is.na(date_death) | viral_breakthrough == 1 | 
                !is.na(graft_loss_date_l) | !is.na(sae_date) ~ 1,
              .default = 0
            ))
  
  #define time contribution
  outcomes_data <- outcomes_data |> 
    mutate(
      time_contr_death_hivB_graftloss_sae =
        case_when(
          death_hivB_graftloss_sae == 1 ~ as.numeric(
            min(
              c_across(c("date_death", "date_test1", "graft_loss_date_l",
                        "sae_date")), 
              na.rm = TRUE) - dot,
            units = "days"),
          death_hivB_graftloss_sae == 0 ~ as.numeric(
            min(c_across(end_date), na.rm = TRUE) - dot,units = "days"))) |> 
    mutate(time_contr_death_hivB_graftloss_sae = as.integer((time_contr_death_hivB_graftloss_sae/31.44))) |>
    mutate(time_contr_death_hivB_graftloss_sae = 
            if_else(as.integer(time_contr_death_hivB_graftloss_sae) == 0, 0.5,
                    time_contr_death_hivB_graftloss_sae))
  
  
  # define outcome death, hiv breakthrough, graft loss, SAE, or infection
  # this is the same as the composite outcome since there is no HIV virologic failure
  outcomes_data <- outcomes_data |> 
    mutate(death_hivB_graftloss_sae_infection = composite,
           time_contr_death_hivB_graftloss_sae_infection = time_contr_comp
            )
  # create rejection
  outcomes_data <- outcomes_data |> 
    mutate(rejection = 
            case_when(
              !is.na(rejection_date_liver) ~ 1,
              .default = 0
           ))
  #define time contribution
  outcomes_data <- outcomes_data |> 
    mutate(time_contr_rejection = case_when(
            rejection == 1 ~ as.numeric(rejection_date_liver - dot, units = "days")/31.44,
            rejection == 0 ~ as.numeric(min(c_across(end_date), 
                          na.rm = TRUE) - dot, units = "days")/31.44))
  
  #iunclude acgf in sensitivity analysis
  outcomes_data <- outcomes_data |> 
    #define all cause of gratf loss
    mutate(acgf = 
            case_when(
              !is.na(date_death)| !is.na(graft_loss_date_l) ~ 1,
              .default = 0
            ))

  outcomes_data <- outcomes_data |> 
    mutate(acgf_date = pmin(date_death, graft_loss_date_l, na.rm = TRUE))


  #define time contribution
  outcomes_data <- outcomes_data |> 
    mutate(time_contr_acgf = case_when(
            acgf == 1 ~ as.numeric(acgf_date - dot, units = "days"),
            acgf == 0 ~ as.numeric(min(c_across(end_date), 
                          na.rm = TRUE) - dot,
                        units = "days")))
            
  #change to months
  outcomes_data <- outcomes_data |> 
    mutate(time_contr_acgf = as.integer((time_contr_acgf/31.44))) |>
    #add 0.5 for time = 0
    mutate(time_contr_acgf = 
            if_else(as.integer(time_contr_acgf) == 0, 0.5,
                    time_contr_acgf)) |>
    #censoring at 48 monhts
    mutate(time_contr_acgf = if_else(time_contr_acgf > 48, 48, time_contr_acgf),
          acgf = if_else(time_contr_acgf == 48 & acgf == 0, 0, acgf))
  
  # add extra variables to fib4 data
  #add slk information
  fib4_index <- fib4_index |> 
    left_join(hope_liver |> select(studyid, is_slk, tx_age, gender, group), by = "studyid")

  # add dot and end date to fib4 data
  fib4_index <- fib4_index |> 
    left_join(outcomes_data |> select(studyid, dot, date_remove, date_death, graft_loss_date_l, censoring_date), by = "studyid")


  #correct age according to time point
  fib4_index <- fib4_index %>%
    mutate(age_at_timepoint = case_when(
      redcap_event_name == "week_26_arm_1" ~ tx_age + 0.5,
      redcap_event_name == "week_52_year_1_arm_1" ~ tx_age + 1,
      redcap_event_name == "week_78_arm_1" ~ tx_age + 1.5,
      redcap_event_name == "week_104_year_2_arm_1" ~ tx_age + 2,
      redcap_event_name == "week_130_arm_1" ~ tx_age + 2.5,
      redcap_event_name == "week_156_year_3_arm_1" ~ tx_age + 3,
      redcap_event_name == "week_182_arm_1" ~ tx_age + 3.5,
      redcap_event_name == "week_208_year_4_arm_1" ~ tx_age + 4,
      TRUE ~ tx_age
    ))
  
  #change values in ast alt and paltet to numeric
  fib4_index <- fib4_index %>%
    mutate(
      ast_fu = as.numeric(str_extract_all(ast_fu,"\\d+\\.?\\d*", simplify = TRUE)),
      alt_fu = as.numeric(str_extract_all(alt_fu,"\\d+\\.?\\d*", simplify = TRUE)),
      bilirubin_fu = as.numeric(str_extract_all(bilirubin_fu,"\\d+\\.?\\d*", simplify = TRUE))
    )
  #label units
  fib4_index <- fib4_index %>%
    mutate(
      alt_unit_fu = case_when(
        alt_unit_fu == 1 ~ "U/L",
        alt_unit_fu == 2 ~ "mU/mL",
        TRUE ~ NA_character_
      ),
      ast_unit_fu = case_when(
        ast_unit_fu == 1 ~ "U/L",
        ast_unit_fu == 2 ~ "mU/mL",
        TRUE ~ NA_character_
      ),
      platelet_unit_fu = case_when(
        platelet_unit_fu == 1 ~ "K/cu mm",
        platelet_unit_fu == 2 ~ "10^9/L",
        platelet_unit_fu == 3 ~ "10^3/mcL",
        TRUE ~ NA_character_
      )
    )
  
  #convert units if needed to compute fib-4
  fib4_index <- fib4_index %>%
    mutate(
      #ast should be in U/L
      ast_converted = case_when(
        ast_unit_fu == "U/L" ~ ast_fu,
        ast_unit_fu == "mU/mL" ~ ast_fu * 1,
        TRUE ~ NA_real_
      ),
      #alt should be in U/L
      alt_converted = case_when(
        alt_unit_fu == "U/L" ~ alt_fu,
        alt_unit_fu == "mU/mL" ~ alt_fu * 1,
        TRUE ~ NA_real_
      ),
      #platelet should be 10^3/microL
      platelet_converted = case_when(
        platelet_unit_fu == "K/cu mm" ~ platelet_fu,
        platelet_unit_fu == "10^9/L" ~ platelet_fu,
        platelet_unit_fu == "10^3/mcL" ~ platelet_fu,
        TRUE ~ NA_real_
      )
    )
  
  #compute fib-4 index
  fib4_index <- fib4_index %>%
    mutate(
      fib4 = case_when(
        !is.na(ast_converted) & !is.na(alt_converted) & !is.na(bilirubin_fu) & !is.na(platelet_converted) & !is.na(age_at_timepoint) ~
          (age_at_timepoint * ast_converted) / (platelet_converted * sqrt(alt_converted)),
        TRUE ~ NA_real_
      )
    )
  
  #order timeppoints for model
  fib4_index <- fib4_index %>%
    mutate(
      time_points = factor(
        redcap_event_name, 
        levels = c("week_26_arm_1", "week_52_year_1_arm_1","week_78_arm_1",
                  "week_104_year_2_arm_1","week_130_arm_1", "week_156_year_3_arm_1",
                  "week_182_arm_1", "week_208_year_4_arm_1"),
      labels = c("6 mo","12 mo","18 mo","24 mo","30 mo","36 mo",
                "42 mo","48 mo"))
    )

  return(
    list(
      demographics_info_lr = hope_liver,
      demographics_info_don = donors,
      outcomes_data = outcomes_data,
      infection_data = infection_data,
      sae_data = sae_data,
      biomarkers_data = fib4_index
    )
  )
}

# create a function that map dead
label_dead <- function(studyid) {
  case_when(
    studyid == "L201-484R" ~ "Cardiac Thrombus",
    studyid == "L210-772R" ~ "Cardiac Arrest",
    studyid == "L223-617R" ~ "Unknown - subject passed at OS LTACH; no records available. See AE log #1.",
    studyid == "L201-094R" ~ "Unknown",
    studyid == "L220-967R" ~ "Unknown",
    studyid == "L201-738R" ~ "Aspiration",
    studyid == "L201-655R" ~ "Metastatic Hepatocellular Carcinoma",
    studyid == "L201-547R" ~ "Recurrent HCC",
    studyid == "L215-781R" ~ "Stroke",
    studyid == "L203-399R" ~ "Unknown",
    .default = NA
  )
}


# create a function to define site of infection
label_condition <- function(code) {
  case_when(
    code == 1 ~ "CNS",
    code == 2 ~ "Eye",
    code == 3 ~ "Ear/Nose/Sinus",
    code == 4 ~ "Mouth/Neck",
    code == 5 ~ "Pulmonary",
    code == 6 ~ "Cardiovascular",
    code == 7 ~ "Intraabdominal",
    code == 8 ~ "GU tract",
    code == 9 ~ "Skeletal/joint",
    code == 10 ~ "Skin",
    code == 11 ~ "Lymph",
    code == 12 ~ "Disseminated infection",
    code == 13 ~ "Fever, unknown origin",
    code == 14 ~ "Neutropenic fever",
    code == 16 ~ "GI tract, excluding anal/rectal",
    code == 20 ~ "Anal/Rectal",
    code == 18 ~ "Blood stream",
    code == 0 ~ "Unknown",
    .default = NA
  )
}


# create a function to define type of infection
label_infection_type <- function(code) {
  case_when(
    code == 1 ~ "Bacterial",
    code == 2 ~ "Parasitic",
    code == 3 ~ "Viral",
    code == 4 ~ "Fungal",
    code == 5 ~ "Unknown",
    .default = NA
  )
}

# create a function to map fungal infection
# Create the function that maps the codes to fungal labels
label_fungal_species <- function(code) {
  case_when(
    code == 1 ~ "Absidia spp.",
    code == 2 ~ "Alternaria spp.",
    code == 3 ~ "Aspergillus fumigatus",
    code == 4 ~ "Aspergillus nidulans",
    code == 5 ~ "Aspergillus niger",
    code == 6 ~ "Aspergillus species, other",
    code == 7 ~ "Aspergillus terreus",
    code == 8 ~ "Blastomyces dermatitidis",
    code == 9 ~ "Candida albicans",
    code == 10 ~ "Candida dubliniensis",
    code == 11 ~ "Candida glabrata",
    code == 12 ~ "Candida guilliermondii",
    code == 13 ~ "Candida kefyr",
    code == 14 ~ "Candida krusei",
    code == 15 ~ "Candida lusitaniae",
    code == 16 ~ "Candida parapsilosis",
    code == 17 ~ "Candida species, other",
    code == 18 ~ "Candida tropicalis",
    code == 19 ~ "Coccidioides immitis / posadasii",
    code == 20 ~ "Cryptococcus gattii",
    code == 21 ~ "Cryptococcus neoformans",
    code == 22 ~ "Cryptococcus other",
    code == 23 ~ "Cunninghamella spp.",
    code == 24 ~ "Dermatophytes, NOS",
    code == 25 ~ "Fusarium spp.",
    code == 26 ~ "Histoplasma capsulatum",
    code == 27 ~ "Mucor spp.",
    code == 28 ~ "Paracoccidioides",
    code == 29 ~ "Pneumocystis jirovecii",
    code == 30 ~ "Rhizomucor spp.",
    code == 31 ~ "Rhizopus spp.",
    code == 32 ~ "Scedosporium apiospermum / Pseudallescheria boydii",
    code == 33 ~ "Scedosporium prolificans",
    code == 34 ~ "Sporothrix schenckii",
    code == 35 ~ "Zygomycetes, NOS",
    code == 36 ~ "Other",
    code == 0 ~ "None",
    .default = NA
  )
}

# create a function to map bacterias
label_bacterial_species <- function(code) {
  case_when(
    code == 1 ~ "Acinetobacter baumannii",
    code == 2 ~ "Actinomyces",
    code == 3 ~ "Aeromonas",
    code == 4 ~ "Bacillus anthracis",
    code == 5 ~ "Bacillus species",
    code == 6 ~ "Bacteroides fragilis",
    code == 7 ~ "Bacteroides species",
    code == 8 ~ "Bartonella species",
    code == 9 ~ "Bordetella species",
    code == 10 ~ "Borrelia burgdorferi",
    code == 11 ~ "Borrelia species",
    code == 12 ~ "Brucella Species",
    code == 13 ~ "Burkholderia cepacia",
    code == 14 ~ "Burkholderia mallei",
    code == 15 ~ "Burkholderia pseudomallei",
    code == 16 ~ "Campylobacter and related species",
    code == 17 ~ "Campylobacter jejuni",
    code == 18 ~ "Capnocytophaga canimorsus",
    code == 19 ~ "Chlamydia trachomatis",
    code == 20 ~ "Chlamydophila pneumonia",
    code == 21 ~ "Chlamydophila psittaci",
    code == 22 ~ "Citrobacter species",
    code == 23 ~ "Clostridium botulinum",
    code == 24 ~ "Clostridium difficile",
    code == 25 ~ "Clostridium species",
    code == 26 ~ "Clostridium tetani (Tetanus)",
    code == 27 ~ "Corynebacterium diphtheriae",
    code == 28 ~ "Corynebacterium species",
    code == 29 ~ "Coxiella burnetii",
    code == 30 ~ "Ehrlichia species",
    code == 31 ~ "Eikenella corrodens",
    code == 32 ~ "Enterobacter species",
    code == 33 ~ "Enterococcus",
    code == 34 ~ "Enterococcus faecalis",
    code == 35 ~ "Enterococcus faecium",
    code == 36 ~ "Erysipelothrix rhusiopathiae",
    code == 37 ~ "Escherichia coli",
    code == 38 ~ "Francisella tularensis",
    code == 39 ~ "Haemophilus ducreyi (Chancroid)",
    code == 40 ~ "Haemophilus influenzae",
    code == 41 ~ "Helicobacter cinaedi and related species",
    code == 42 ~ "Helicobacter pylori",
    code == 43 ~ "Klebsiella oxytoca",
    code == 44 ~ "Klebsiella pneumoniae",
    code == 45 ~ "Klebsiella species",
    code == 46 ~ "Lactobacillus",
    code == 47 ~ "Legionella pneumophila",
    code == 48 ~ "Legionella species",
    code == 49 ~ "Leptospira interrogans",
    code == 50 ~ "Listeria monocytogenes",
    code == 51 ~ "Lymphogranuloma venereum (LGV)",
    code == 52 ~ "Mixed flora",
    code == 53 ~ "Moraxella catarrhalis",
    code == 54 ~ "Morganella",
    code == 55 ~ "Mycobacterium abscessus",
    code == 56 ~ "Mycobacterium avium-complex (MAC, MAI, non-HIV)",
    code == 57 ~ "Mycobacterium chelonae",
    code == 58 ~ "Mycobacterium fortuitum",
    code == 59 ~ "Mycobacterium gordonae",
    code == 60 ~ "Mycobacterium kansasii",
    code == 61 ~ "Mycobacterium leprae",
    code == 62 ~ "Mycobacterium marinum",
    code == 63 ~ "Mycobacterium scrofulaceum",
    code == 64 ~ "Mycobacterium tuberculosis",
    code == 65 ~ "Mycobacterium ulcerans",
    code == 66 ~ "Mycobacterium xenopi",
    code == 67 ~ "Mycobacterium, other atypical",
    code == 68 ~ "Mycoplasma pneumoniae (Antibiotic Guide)",
    code == 69 ~ "Mycobacteria, other atypical",
    code == 70 ~ "Neisseria gonorrhoeae",
    code == 71 ~ "Neisseria meningitidis",
    code == 72 ~ "Nocardia farcinica",
    code == 73 ~ "Nocardia nova",
    code == 74 ~ "Nocardia asteroides",
    code == 75 ~ "Nocardia braziliensis",
    code == 76 ~ "Nocardia spp, other",
    code == 77 ~ "Pasteurella multocida",
    code == 78 ~ "Peptostreptococcus/Peptococcus",
    code == 79 ~ "Plesiomonas",
    code == 80 ~ "Propionibacterium species",
    code == 81 ~ "Proteus species",
    code == 82 ~ "Providencia",
    code == 83 ~ "Pseudomonas aeruginosa",
    code == 84 ~ "Pseudomonas spp.",
    code == 85 ~ "Rhodococcus equi",
    code == 86 ~ "Rickettsia rickettsii",
    code == 87 ~ "Rickettsia species",
    code == 88 ~ "Salmonella species",
    code == 89 ~ "Serratia species",
    code == 90 ~ "Shigella dysenteriae",
    code == 91 ~ "Shigella species",
    code == 92 ~ "Staphylococci, coagulase negative",
    code == 93 ~ "Staphylococcus aureus",
    code == 94 ~ "Stenotrophomonas maltophilia",
    code == 95 ~ "Streptobacillus moniliformis",
    code == 96 ~ "Streptococcus pneumoniae",
    code == 97 ~ "Streptococcus pyogenes (Group A)",
    code == 98 ~ "Streptococcus species",
    code == 99 ~ "Treponema pallidum (syphilis)",
    code == 100 ~ "Tropheryma whipplei",
    code == 101 ~ "Vibrio cholerae",
    code == 102 ~ "Vibrio species (non-cholera)",
    code == 103 ~ "Yersinia pestis",
    code == 104 ~ "Yersinia species (non-plague)",
    code == 105 ~ "Other",
    code == 0 ~ "None",
    .default = NA
  )
}

# create a function to map viral infection
label_viral_species <- function(code) {
  case_when(
    code == 1 ~ "Adenovirus",
    code == 2 ~ "B Virus",
    code == 3 ~ "BK Virus",
    code == 4 ~ "Chikungunya virus",
    code == 5 ~ "Coronavirus (all others except COVID-19)",
    code == 56 ~ "COVID-19 (SARS-CoV-2)",
    code == 6 ~ "Cytomegalovirus",
    code == 7 ~ "Dengue Virus",
    code == 8 ~ "Ebola virus",
    code == 9 ~ "Enterovirus",
    code == 10 ~ "Epstein-Barr Virus",
    code == 11 ~ "Hantavirus",
    code == 12 ~ "Hepatitis A",
    code == 13 ~ "Hepatitis B",
    code == 14 ~ "Hepatitis C",
    code == 15 ~ "Hepatitis D",
    code == 16 ~ "Hepatitis E",
    code == 17 ~ "Herpes Simplex Virus",
    code == 18 ~ "HHV 6",
    code == 19 ~ "HHV 7",
    code == 20 ~ "HHV 8",
    code == 21 ~ "HIV",
    code == 22 ~ "HTLV I/II",
    code == 23 ~ "Human Metapneumovirus (MPV)",
    code == 24 ~ "Human papillomavirus (HPV)",
    code == 25 ~ "Influenza",
    code == 26 ~ "Influenza A",
    code == 27 ~ "Influenza B",
    code == 28 ~ "Japanese Encephalitis Virus",
    code == 29 ~ "JC Virus",
    code == 30 ~ "Kyasanur virus",
    code == 31 ~ "Lassa Fever virus",
    code == 32 ~ "Marburg virus",
    code == 33 ~ "Measles",
    code == 34 ~ "Molluscum contagiosum",
    code == 35 ~ "Mumps",
    code == 36 ~ "Norovirus",
    code == 37 ~ "Omsk hemorrhagic fever virus",
    code == 38 ~ "Parainfluenza virus",
    code == 39 ~ "Parainfluenza virus-1",
    code == 40 ~ "Parainfluenza virus -2",
    code == 41 ~ "Parainfluenza virus -3",
    code == 42 ~ "Parainfluenza virus -4",
    code == 43 ~ "Parvovirus B19",
    code == 44 ~ "Poliovirus",
    code == 45 ~ "Rabies",
    code == 46 ~ "Respiratory syncytial virus",
    code == 47 ~ "Rhinovirus",
    code == 48 ~ "Rift Valley Virus",
    code == 49 ~ "Rubella",
    code == 50 ~ "Smallpox",
    code == 51 ~ "Vaccinia",
    code == 52 ~ "Varicella-zoster virus",
    code == 53 ~ "West Nile Virus",
    code == 54 ~ "Yellow fever virus",
    code == 55 ~ "Other",
    code == 0 ~ "None",
    .default = NA
  )
}

# create a function to map parasite infecion
label_parasitic_species <- function(code) {
  case_when(
    code == 1 ~ "Acanthamoeba",
    code == 2 ~ "Angiostrongylus",
    code == 3 ~ "Ascaris",
    code == 4 ~ "Babesia Species",
    code == 5 ~ "Baylisascariasis",
    code == 6 ~ "Cryptosporidia",
    code == 7 ~ "Cyclospora cayetanensis",
    code == 8 ~ "Cystoisospora belli",
    code == 9 ~ "Diphyllobothriasis",
    code == 10 ~ "Echinococcus spp.",
    code == 11 ~ "Entamoeba histolytica",
    code == 12 ~ "Enterobius (Pinworm)",
    code == 13 ~ "Fasciola hepatica",
    code == 14 ~ "Filariasis",
    code == 15 ~ "Giardia lamblia",
    code == 16 ~ "Hookworm",
    code == 17 ~ "Leishmania species",
    code == 18 ~ "Lice",
    code == 19 ~ "Loa Loa",
    code == 20 ~ "Microsporidium",
    code == 21 ~ "Onchocerciasis",
    code == 22 ~ "Plasmodium",
    code == 23 ~ "Plasmodium falciparum",
    code == 24 ~ "Plasmodium malariae",
    code == 25 ~ "Plasmodium ovale",
    code == 26 ~ "Plasmodium vivax",
    code == 27 ~ "Sarcoptes scabiei var hominis (Scabies)",
    code == 28 ~ "Schistosoma species",
    code == 29 ~ "Strongyloides stercoralis",
    code == 30 ~ "Taenia saginata",
    code == 31 ~ "Taenia solium",
    code == 32 ~ "Toxocariasis",
    code == 33 ~ "Toxoplasma gondii",
    code == 34 ~ "Trichinella species",
    code == 35 ~ "Trichomonas vaginalis",
    code == 36 ~ "Trypanosoma brucei gambiense/rhodesiense",
    code == 37 ~ "Trypanosoma cruzi",
    code == 38 ~ "Wucheria bancrofti",
    code == 39 ~ "Other",
    code == 0 ~ "None",
    .default = NA
  )
}

# Liver Donor Risk Index (LDRI / DRI) - Feng et al. 2006
# Inputs:
#  age   = donor age (years)
#  cod   = cause of death: "trauma", "anoxia", "cva", "other"
#  eth   = donor ethnicity: "white", "black", "other"
#  dcd   = 1 if DCD, 0 otherwise
#  split = 1 if split/partial graft, 0 otherwise
#  share = "local", "regional", or "national"
#  cit   = cold ischemia time (hours)
#  height= donor height (cm)

liver_dri_feng2006 <- function(age, cod, eth, dcd, split, share, cit, height) {
  # Basic checks
  n <- length(age)
  stopifnot(
    length(cod) == n, length(eth) == n, length(dcd) == n, length(split) == n,
    length(share) == n, length(cit) == n, length(height) == n
  )

  cod <- tolower(cod)
  eth <- tolower(eth)
  share <- tolower(share)

  # Age category effects
  agevar <- ifelse(age < 40, 0,
    ifelse(age < 50, 0.154,
      ifelse(age < 60, 0.274,
        ifelse(age < 70, 0.424, 0.501)
      )
    )
  )

  # Cause of death effects
  codvar <- ifelse(cod == "anoxia", 0.079,
    ifelse(cod == "cva", 0.145,
      ifelse(cod == "trauma", 0, 0.184)
    )
  )

  # Race/ethnicity effects
  racevar <- ifelse(eth == "black", 0.176,
    ifelse(eth == "white", 0, 0.126)
  )

  # Other components
  dcdvar <- 0.411 * dcd
  splitvar <- 0.422 * split
  sharevar <- ifelse(share == "local", 0,
    ifelse(share == "regional", 0.105, 0.244)
  )
  citvar <- 0.010 * (cit - 8)
  heightvar <- 0.066 * (170 - height) / 10

  # DRI (multiplicative) = exp(sum of components)
  exp(agevar + codvar + racevar + dcdvar + splitvar + sharevar + citvar + heightvar)
}


# create a funciton to generate infection data
# get all infections
infection_long <- function(data, events = NULL, origin = NULL) {
  df_date <- data |>
    dplyr::select(studyid, matches(events), infection1_aids, resolution) |>
    pivot_longer(
      cols = infection1_date, values_to = "infection_date",
      names_to = "infection_event"
    ) |>
    # create a key to ensure that we get a correct ordem in infection type
    mutate(key = 1:length(studyid)) |>
    dplyr::select(-infection_event)

  df_type <- data |>
    dplyr::select(
      studyid,
      matches("infection[0-9]{1,2}_type")
    ) |>
    pivot_longer(
      cols = infection1_type, values_to = "infection_type",
      names_to = "type"
    ) |>
    # create a key to ensure that we get a correct ordem in infection type
    mutate(key = 1:length(studyid)) |>
    dplyr::select(-type)


  # add site of infection
  df_site <- data |>
    dplyr::select(
      studyid,
      matches("infection[0-9]{1,2}_site")
    ) |>
    pivot_longer(
      cols = infection1_site, values_to = "infection_site",
      names_to = "site"
    ) |>
    # create a key to ensure that we get a correct ordem in infection type
    mutate(key = 1:length(studyid)) |>
    dplyr::select(-site)


  # Add if infection leads to hospitalization
  df_hosp <- data |>
    dplyr::select(
      studyid,
      matches("infection[0-9]{1,2}_hospital")
    ) |>
    pivot_longer(
      cols = infection1_hospital, values_to = "hospitalization",
      names_to = "foo"
    ) |>
    # create a key to ensure that we get a correct ordem in infection type
    mutate(key = 1:length(studyid)) |>
    dplyr::select(-foo)


  # combine all data
  df_date <- df_date |> left_join(df_type)
  df_date <- df_date |> left_join(df_site)
  df_date <- df_date |> left_join(df_hosp)


  # label infection site
  df_date <- df_date |>
    mutate(
      infection_site = label_condition(infection_site),
      infection_type = label_infection_type(infection_type)
    )


  # add dot to data
  df_date <- df_date |>
    left_join(data |>
      filter(!is.na(dot)) |>
      dplyr::select(studyid, dot, group))


  df_pathogen <- data |>
    select(
      studyid,
      matches("bacterial1_pathogen[0-9]"),
      matches("viral1_pathogen[0-9]"),
      matches("fungal1_pathogen[0-9]"),
      matches("parasitic1_pathogen[0-9]")
    ) |>
    # create a key to ensure that we get a correct order in infection type
    mutate(key = 1:length(studyid))


  # change collum names
  colnames(df_pathogen) <- c(
    "studyid", "first_bacteria", "second_bacteria",
    "third_bacteria", "first_virus", "second_virus",
    "third_virus", "first_fungal", "second_fungal",
    "third_fungal", "first_parasit",
    "second_parasit", "third_parasit", "key"
  )

  # add pathogen information
  df_date <- df_date |>
    left_join(df_pathogen, by = c("studyid", "key")) |>
    # label fungal species
    mutate(
      first_fungal = label_fungal_species(first_fungal),
      second_fungal = label_fungal_species(second_fungal),
      third_fungal = label_fungal_species(third_fungal),
      first_bacteria = label_bacterial_species(first_bacteria),
      second_bacteria = label_bacterial_species(second_bacteria),
      third_bacteria = label_bacterial_species(third_bacteria),
      first_virus = label_viral_species(first_virus),
      second_virus = label_viral_species(second_virus),
      third_virus = label_viral_species(third_virus),
      first_parasit = label_parasitic_species(first_parasit),
      second_parasit = label_parasitic_species(second_parasit),
      third_parasit = label_parasitic_species(third_parasit)
    )


  return(df_date)
}

# define a function to get year
get_year <- function(year_code) {
  return(year_value = case_when(
    year_code == 0 ~ NA, # Code 0 is converted to NA
    year_code == 1 ~ 2028,
    year_code == 2 ~ 2027,
    year_code == 3 ~ 2026,
    year_code == 4 ~ 2025,
    year_code == 5 ~ 2024,
    year_code == 6 ~ 2023,
    year_code == 7 ~ 2022,
    year_code == 8 ~ 2021,
    year_code == 9 ~ 2020,
    year_code == 10 ~ 2019,
    year_code == 11 ~ 2018,
    year_code == 12 ~ 2017,
    year_code == 13 ~ 2016,
    year_code == 14 ~ 2015,
    year_code == 15 ~ 2014,
    year_code == 16 ~ 2013,
    year_code == 17 ~ 2012,
    year_code == 18 ~ 2011,
    year_code == 19 ~ 2010,
    year_code == 20 ~ 2009,
    year_code == 21 ~ 2008,
    year_code == 22 ~ 2007,
    year_code == 23 ~ 2006,
    year_code == 24 ~ 2005,
    year_code == 25 ~ 2004,
    year_code == 26 ~ 2003,
    year_code == 27 ~ 2002,
    year_code == 28 ~ 2001,
    year_code == 29 ~ 2000,
    year_code == 30 ~ 1999,
    year_code == 31 ~ 1998,
    year_code == 32 ~ 1997,
    year_code == 33 ~ 1996,
    year_code == 34 ~ 1995,
    year_code == 35 ~ 1994,
    year_code == 36 ~ 1993,
    year_code == 37 ~ 1992,
    year_code == 38 ~ 1991,
    year_code == 39 ~ 1990,
    year_code == 40 ~ 1989,
    year_code == 41 ~ 1988,
    year_code == 42 ~ 1987,
    year_code == 43 ~ 1986,
    year_code == 44 ~ 1985,
    year_code == 45 ~ 1984,
    year_code == 46 ~ 1983,
    year_code == 47 ~ 1982,
    year_code == 48 ~ 1981,
    year_code == 49 ~ 1980,
    year_code == 50 ~ 1979,
    year_code == 51 ~ 1978,
    year_code == 52 ~ 1977,
    year_code == 53 ~ 1976,
    year_code == 54 ~ 1975,
    year_code == 55 ~ 1974,
    year_code == 56 ~ 1973,
    year_code == 57 ~ 1972,
    year_code == 58 ~ 1971,
    year_code == 59 ~ 1970,
    year_code == 60 ~ 1969,
    year_code == 61 ~ 1968,
    year_code == 62 ~ 1967,
    year_code == 63 ~ 1966,
    year_code == 64 ~ 1965,
    year_code == 65 ~ 1964,
    year_code == 66 ~ 1963,
    year_code == 67 ~ 1962,
    year_code == 68 ~ 1961,
    year_code == 69 ~ 1960,
    year_code == 70 ~ 1959,
    year_code == 71 ~ 1958,
    year_code == 72 ~ 1957,
    year_code == 73 ~ 1956,
    year_code == 74 ~ 1955,
    year_code == 75 ~ 1954,
    year_code == 76 ~ 1953,
    year_code == 77 ~ 1952,
    year_code == 78 ~ 1951,
    year_code == 79 ~ 1950,
    year_code == 80 ~ 1949,
    year_code == 81 ~ 1948,
    year_code == 82 ~ 1947,
    year_code == 83 ~ 1946,
    year_code == 84 ~ 1945,
    year_code == 85 ~ 1944,
    year_code == 86 ~ 1943,
    year_code == 87 ~ 1942,
    year_code == 88 ~ 1941,
    year_code == 89 ~ 1940,
    year_code == 90 ~ 1939,
    year_code == 91 ~ 1938,
    .default = NA
  ))
}


# CKD-EPI 2021 (race-free) eGFR from serum creatinine
# Output: eGFR in mL/min/1.73 m^2
#
# Ref (equation form):
# eGFR = 142 * min(Scr/k, 1)^a * max(Scr/k, 1)^(-1.200) * 0.9938^Age * (1.012 if female)
# where:
#   k = 0.7 (female), 0.9 (male)
#   a = -0.241 (female), -0.302 (male)

egfr_ckdepi2021 <- function(scr,
                            age,
                            sex = c("female", "male"),
                            scr_unit = c("mg/dL", "umol/L")) {
  sex <- match.arg(tolower(sex), c("female", "male"))
  scr_unit <- match.arg(scr_unit)

  # Convert creatinine to mg/dL if needed
  scr_mgdl <- if (scr_unit == "umol/L") scr / 88.4 else scr

  # Sex-specific constants
  k <- if (sex == "female") 0.7 else 0.9
  a <- if (sex == "female") -0.241 else -0.302
  sex_factor <- if (sex == "female") 1.012 else 1.0

  # Core terms
  ratio <- scr_mgdl / k
  egfr <- 142 *
    (pmin(ratio, 1)^a) *
    (pmax(ratio, 1)^(-1.200)) *
    (0.9938^age) *
    sex_factor

  return(egfr)
}

# create a function that organize data to run Win ratio analysis
generate_winratio_data <- function(data) {
  data |>
    mutate(

      # Death occurred if a death date is present.
      death_status = as.integer(!is.na(date_death)),

      # If death occurred, use the death date.
      # Otherwise use end_date as the censoring date.
      death_obs_date = if_else(
        death_status == 1,
        date_death,
        end_date
      ),

      # Time from transplant to death/censoring, in days.
      death_time = as.numeric(
        death_obs_date - dot
      ),

    
      # In this dataset all observed breakthrough
      # dates correspond to reported breakthrough.
      hiv_status = as.integer(!is.na(date_test1)),

      hiv_obs_date = if_else(
        hiv_status == 1,
        date_test1,
        end_date
      ),

      hiv_time = as.numeric(
        hiv_obs_date - dot
      ),

   
      # There is no separate graft-loss indicator,
      # so presence of the date defines the event.
      graft_status = as.integer(
        !is.na(graft_loss_date_l)
      ),

      graft_obs_date = if_else(
        graft_status == 1,
        graft_loss_date_l,
        end_date
      ),

      graft_time = as.numeric(
        graft_obs_date - dot
      ),

   
      # Presence of sae_date defines the SAE event.
      sae_status = as.integer(
        !is.na(sae_date)
      ),

      sae_obs_date = if_else(
        sae_status == 1,
        sae_date,
        end_date
      ),

      sae_time = as.numeric(
        sae_obs_date - dot
      ),


      # Only infection1_aids == 1 is considered
      # an opportunistic infection.
       oi_status = as.integer(
        coalesce(infection1_aids, 0) == 1),

      oi_obs_date = if_else(
        oi_status == 1,
        infection1_date,
        end_date
      ),

      oi_time = as.numeric(
        oi_obs_date - dot
      )
    ) |> 
    # change all NA values in status to 0
    mutate(
      across(
        c(death_status, hiv_status, graft_status, sae_status, oi_status),
        ~ replace_na(.x, 0)
      )
    )
}
