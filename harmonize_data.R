##Created by Moreno Rodrigues
##June 2026
'
This file dispatches to the correct dataset-specific harmonization function.
Source this file once, then call harmonize_data(data, dataset = "donors").
'

harmonize_data <- function(data, dataset, 
                           harmonization_dir = 
                             "//win.ad.jhu.edu/Users$/HOME/Data_Harmonization/harmonization_files") {
  dataset_key <- tolower(trimws(dataset))

  harmonizers <- list(
    donors = list(file = "donor_harmonization.R", function_name = "harmonize_donors"),
    kltf = list(file = "kltf_harmonization.R", function_name = "harmonize_kltf"),
    klto = list(file = "klto_harmonization.R", function_name = "harmonize_klto"),
    lltf = list(file = "lltf_harmonization.R", function_name = "harmonize_lltf"),
    pilot = list(file = "pilot_harmonization.R", function_name = "harmonize_pilot"),
    u01k = list(file = "u01k_harmonization.R", function_name = "harmonize_u01k"),
    u01l = list(file = "u01l_harmonization.R", function_name = "harmonize_u01l"),
    strive = list(file = "strive_harmonization.R", function_name = "harmonize_strive"),
    epoc = list(file = "epoc_harmonization.R", function_name = "harmonize_epoc"),
    epoc_vaccine = list(file = "create_covid_vac_dates.R", 
                         function_name = "compile_covid_vax_dates")
  )

  harmonizer <- harmonizers[[dataset_key]]

  if (is.null(harmonizer)) {
    stop(
      "Unknown dataset: '", dataset, "'. Available datasets are: ",
      paste(sort(names(harmonizers)), collapse = ", "),
      call. = FALSE
    )
  }

  harmonization_file <- file.path(harmonization_dir, harmonizer$file)

  if (!file.exists(harmonization_file)) {
    stop("Cannot find harmonization file: ", harmonization_file, call. = FALSE)
  }

  harmonization_env <- new.env(parent = parent.frame())
  source(harmonization_file, local = harmonization_env)

  harmonization_function <- get0(
    harmonizer$function_name,
    envir = harmonization_env,
    mode = "function",
    inherits = FALSE
  )

  if (is.null(harmonization_function)) {
    stop(
      "Cannot find function '", harmonizer$function_name,
      "' in ", harmonization_file,
      call. = FALSE
    )
  }

  harmonization_function(data)
}
