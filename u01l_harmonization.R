##Created by Moreno Rodrigues
##March 2026
"
This file is used to standardize the U01L dataset for the purpose of harmonization 
with other datasets in the project. It includes functions to recode 
recipient characteristics, and to clean specific variables such as rejection data. 
The recoding is based on the original variable names and values in the U01K dataset, 
and it aims to create standardized variables that can be easily 
compared across different datasets.
"


#define a function to get viral load
calculate_viral_load <- function(hiv_pcr) {
  return(case_when(
    str_detect(str_to_lower(hiv_pcr), "detected") ~ 0,
    str_detect(str_to_lower(hiv_pcr), "positive") ~ 0,
    str_detect(str_to_lower(hiv_pcr), "undet") ~ 0,
    str_detect(str_to_lower(hiv_pcr), "non") ~ 0,
    TRUE ~ as.numeric(str_extract(hiv_pcr, "\\d+\\.?\\d*"))
  ))
}


# define region where patient was born
define_us_region <- function(data){
  
  data %>% 
    mutate(
      birth_region = case_when(
        northeast == 1 ~ "Connecticut",
        northeast == 2 ~ "Maine",
        northeast == 3 ~ "Massachusetts",
        northeast == 4 ~ "New Hampshire",
        northeast == 5 ~ "Rhode Island Vermont",
        northeast == 6 ~ "New Jersey",
        northeast == 7 ~ "New York",
        northeast == 8 ~ "Pennsylvania",
        midwest == 1 ~ "Indiana",
        midwest == 2 ~ "Illinois",
        midwest == 3 ~ "Michigan",
        midwest == 4 ~ "Ohio",
        midwest == 5 ~ "Wisconsin",
        midwest == 6 ~ "Iowa",
        midwest == 7 ~ "Kansas",
        midwest == 8 ~ "Minnesota",
        midwest == 9 ~ "Missouri",
        midwest == 10 ~ "Nebraska",
        midwest == 11 ~ "North Dakota",
        midwest == 12 ~ "South Dakota",
        south == 1 ~ "Delaware",
        south == 2 ~ "District of Columbia",
        south == 3 ~ "Florida",
        south == 4 ~ "Georgia",
        south == 5 ~ "Maryland",
        south == 6 ~ "South Carolina",
        south == 7 ~ "North Carolina",
        south == 8 ~ "Virginia",
        south == 9 ~ "West Virginia",
        south == 10 ~ "Alabama",
        south == 11 ~ "Kentucky",
        south == 12 ~ "Mississippi",
        south == 13 ~ "Tennessee",
        south == 14 ~ "Arkansas",
        south == 15 ~ "Louisiana",
        south == 16 ~ "Oklahoma",
        south == 17 ~ "Texas",
        west == 1 ~ "Arizona",
        west == 1 ~ "Colorado",
        west == 1 ~ "Idaho",
        west == 1 ~ "New Mexico",
        west == 1 ~ "Montana",
        west == 1 ~ "Utah",
        west == 1 ~ "Nevada",
        west == 1 ~ "Wyoming",
        west == 1 ~ "Alaska",
        west == 1 ~ "California",
        west == 1 ~ "Hawaii",
        west == 1 ~ "Oregon",
        west == 1 ~ "Washington",
        islandareas == 1 ~ "Puerto Rico",
        islandareas == 1 ~ "U.S. Virgin Islands",
        islandareas == 1 ~ "Guam",
        islandareas == 1 ~ "American Samoa",
        islandareas == 1 ~ "Commonwealth of the Northern Mariana Islands"
      )
    )
  
}





harmonize_u01l <- function(u01l){
  
  #create a list of others cause of HIV
  other_cause_of_hiv <- c("hiv_risk_category___6","hiv_risk_category___3",
                          "hiv_risk_category___5", "hiv_risk_category___0")
  
  
  # read data with primary cause of liver transplant
  primary_cause <- read_csv(
    "//win.ad.jhu.edu/Users$/HOME/Data_Harmonization/auxiliary_files/primary_indication_for_LT.csv"
  ) 
  
  # add primary cause to data
  u01l <- u01l %>%  left_join(primary_cause %>% select(studyid, primary_cause_transp))
  u01l <- define_us_region(u01l)

  u01l <- u01l %>%
    #define the group variable based on the organ_hiv_status variable
    mutate(
      group = 
             case_when(
               organ_hiv_status == 1 ~ "HIV D+/R+",
               organ_hiv_status %in% c(2,3) ~ "HIV D-/R+"),
      race_combined =
        case_when(
          race___2 == 1 ~ "Black or African American",
          race___1 == 1 ~ "White",
          race___3 == 1 ~ "Asian",
          race___4 == 1 ~ "American Indian or Alaska Native",
          race___5 == 1 ~ "Native Hawaiian or Other Pacific Islander",
          race___0 == 1 ~ "Other",
          .default = NA,
        ),
      ethnicity = if_else(ethnicity == 1, "Hispanic or Latino", "Not Hispanic"),
      
      gender = if_else(gender == 0, "M", "F", NA),
      
      #Hepatitis B core ab positive
      hep_b_antibody = 
        case_when(
          hbv_core_ab == 1 ~ "Positive",
          hbv_core_ab == 0 ~ "Negative",
          .default =  NA),
      
      #hepatitis B surface antigen
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
      
      #Hepatitis B IGG or total core antibody
      hep_b_igg_core =
        case_when(
          hbv_core_ab == 1 ~ "Positive",
          hbv_core_ab == 0 ~ "Negative",
          .default = NA),
      
      #hepatitis C antibody
      hep_c_antibody = 
        case_when(
          hcv_ab == 1 ~ "Positive",
          hcv_ab == 0 ~ "Negative",
          .default = NA),
      
      #hepatitis C PCR
      hep_c_pcr = 
        case_when(
          hcv_pcr == 1 ~ "Positive",
          hcv_pcr == 0 ~ "Negative",
          .default = NA),
      
      #CMV IgG Ab
      cmv_igg_pos = 
        case_when(
          cmv == 1 ~ "Positive",
          cmv == 0 ~ "Negative",
          .default = NA),
      
      #donor CMV status
      donor_cmv_result =
        case_when(
          donorsero_cmvab == 1 ~ "Positive",
          donorsero_cmvab == 0 ~ "Negative",
          .default = NA),
      
      ebv_antibody = case_when(ebv_capsid == 1 ~ "Positive", 
                               ebv_capsid == 0 ~ "Negative", 
                               .default = NA),
      
      #blood_type
      blood_type = case_when(
        blood_type == 1 ~ "A",
        blood_type == 2 ~ "B",
        blood_type == 3 ~ "AB",
        blood_type == 4 ~ "O",
      ),
      #is patient a slk
      is_slk = case_when(
        donor_organs == 1 ~ "No",
        donor_organs == 2 ~ "Yes"
      ),
      #define diabetes
      diabetes = case_when(
        endocrine___1 == 1 ~ "Yes",
        endocrine___1 == 0 ~ "No",
        .default = "No"
        
      ),
      #define hypertension
      hypertension = case_when(
        renal___8 == 1 | cardiovascular___7 == 1 ~ "Yes",
        renal___8 == 0 & cardiovascular___7 == 0 ~ "No",
        .default = "No"
      ),
      hiv_rna_level = calculate_viral_load(hiv_pcr),
      viral_load_fu = 
        case_when(
          str_detect(str_to_lower(hiv_pcr), "detected") ~ 0,
          str_detect(str_to_lower(hiv_pcr), "undet") ~ 0,
          str_detect(str_to_lower(hiv_pcr), "non") ~ 0,
          .default = as.numeric(
            str_extract_all(
              hiv_pcr,
              "\\d+\\.?\\d*", 
              simplify = TRUE))),
      
      cd4_absolute = as.numeric(str_extract_all(cd4_absolute,"\\d+\\.?\\d*")),
      
      cd4_percent = as.numeric(str_extract_all(cd4_percent,"\\d+\\.?\\d*")),
      
      viral_load_transp_binary = if_else(hiv_rna_level < 200, "<= 200", ">200"),
      cancer = if_else(non_aids_cancer == 1, "Yes", "No", NA),
      msm = if_else(hiv_risk_category___1 == 1, "Yes", "No"),
      hetero_sex_risk = if_else(hiv_risk_category___2 == 1, "Yes", "No"),
      idu = if_else(hiv_risk_category___4 == 1, "Yes", "No"),
      other = case_when(
        any(across(other_cause_of_hiv) == 1) ~ "Yes",
        .default = "No"
      ),
      INSTI =
        case_when(
          art_class___1 == 1 ~ "Yes",
          .default =  "No"),
      NNRTI = 
        case_when(
          art_class___2 == 1 ~ "Yes",
          .default = "No"),
      pi_cobicistat = 
        case_when(
          art_class___3 == 1 ~ "Yes",
          .default = "No"),
      steroids = 
        case_when(
          induction_med___1 == 1 ~ "Yes",
          .default = "No"
        ),
      basiliximab = 
        case_when(
          induction_med___3 == 1 ~ "Yes",
          .default = "No"
        ),
      ks_hist = case_when(
        bl_kaposi %in% c(1,2) ~ "Yes",
        bl_kaposi == 0 ~ "No",
        .default = NA
      ),
      # add renal disorders
      fsgs = if_else(renal___1 == 1, "Yes", "No", NA),
      hivan = if_else(renal___2 == 1, "Yes", "No", NA),
      hfsgs = if_else(renal___3 == 1, "Yes", "No", NA),
      hickd = if_else(renal___4 == 1, "Yes", "No", NA),
      htm = if_else(renal___5 == 1, "Yes", "No", NA),
      glomerulonephritis = if_else(renal___7 == 1, "Yes", "No", NA),
      cystic_kidney = if_else(renal___9 == 1, "Yes", "No", NA),
      
      # Cardiovascular disorders
      arrythmia = if_else(cardiovascular___1 == 1, "Yes", "No", NA),
      aortic_aneurysm = if_else(cardiovascular___2 == 1, "Yes", "No", NA),
      cardiomyopathy = if_else(cardiovascular___3 == 1, "Yes", "No", NA),
      cerebrovascular_disease = if_else(cardiovascular___4 == 1, "Yes", "No", NA),
      coronary_artery_disease = if_else(cardiovascular___5 == 1, "Yes", "No", NA),
      congenital_heart_disease = if_else(cardiovascular___6 == 1, "Yes", "No", NA),
      hyperlipidemia = if_else(cardiovascular___8 == 1, "Yes", "No", NA),
      myocardial_infarction = if_else(cardiovascular___9 == 1, "Yes", "No", NA),
      peripheral_artery_disease = if_else(cardiovascular___10 == 1, "Yes", "No", NA),
      rheumatic_heart_disease = if_else(cardiovascular___11 == 1, "Yes", "No", NA),
      stroke = if_else(cardiovascular___12 == 1, "Yes", "No", NA),
      valvular_heart_disease = if_else(cardiovascular___13 == 1, "Yes", "No", NA),
      
      
      
      
      unos_donor_id = tolower(donor_unos_id_liver),
	    studyid = toupper(studyid),
      trial = "U01L"
      )
  u01l <- u01l %>%
    mutate(born_in_us = 
             case_when(
               place_birth == 1 ~ "Yes",
               !is.na(country_birth) ~ country_birth,
               .default = NA
             )
           ) %>% 
    rename(
      date_transplant = dot,
      date_death = date_death,
      date_graft_loss = graft_loss_date_l,
      date_last_fu = date_last_visit
    )

}


