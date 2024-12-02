
simulate_trial <- function(
  scenario,
  trial_configuration,
  fun.generate_cohort,
  population, #fun.sample_patients_for_scenario,
  model
)
{
  # Define cohort template
  cohort_template <- bind_cols(
    scenario[0, ],
    tibble(
      cohort_time = integer(), 
      num_patients = integer(),
      num_toxicities = integer()
    )
  )

  # Simulation main loop: if no trials are active anymore, the simulation is done
  cohort_iteration_index <- 1

  # Record initial state
  model <- model$set_initial_state(model, trial_configuration)
  
  while(trial_configuration$is_running)
  {
    # Create cohort from template (preserving factor levels for group_id, strata_id)
    current_cohorts <- cohort_template

    # Sample a cohort
    current_dosing <- dplyr::filter(scenario, dose == trial_configuration$current_dose)

    # Ensure dose exists
    assert_that(nrow(current_dosing) > 0)
    
    # Instantiate a cohort with the current dose of this trial
    cohort <- dplyr::bind_rows(cohort_template, current_dosing)

    stopifnot(is.finite(trial_configuration$current_date))
    
    # Sample a cohort (generate_cohort sets num_patients and num_toxicities)
    cohort <- fun.generate_cohort(cohort_iteration_index, cohort, trial_configuration$current_date)
    
    # Update model with the new data
    model <- model$update(model, trial_configuration, cohort, scenario)

    # Determine next trial configuration (next dose / whether to stop trial)
    trial_configuration <- model$get_next_configuration(model)

    cohort_iteration_index <- cohort_iteration_index + 1
  }

  # Return model & configuration
  return(list(
    trial_model = model,
    trial_configuration = trial_configuration
  ))
}

reduce_results_per_dose <- function(current_simulation, scenario_data)
{
  result_tibble <- NULL
  
  # current_data <- current_simulation$trial_model$data
  current_data <- filter(current_simulation$trial_model$data, group_id == "trial")
  final_configuration <- tail(current_simulation$trial_model$trial_configuration, n = 1)
  distinct_doses <- unique(current_data$dose)
  
  for (current_dose in distinct_doses) {
    current_dose_is_MTD <- !is.na(final_configuration$MTD_dose) && (final_configuration$MTD_dose == current_dose)
    
    if (current_dose_is_MTD) 
    {
      MTD_prob_true <- max(filter(scenario_data, dose == final_configuration$MTD_dose)$prob_true)
    } else {
      MTD_prob_true <- NA
    }
    
    current_dose_data <- filter(current_data, dose == current_dose)
    
    
    result_tibble <- bind_rows(
      result_tibble,
      tibble(
        dose = current_dose,
        is_MTD = current_dose_is_MTD,
        num_patients_start = sum(current_dose_data$num_patients_start),
        num_patients = sum(current_dose_data$num_patients),
        num_toxicity_events = sum(current_dose_data$num_toxicities)
      )
    )
  }


  

  return(result_tibble)
}

reduce_results_per_trial <- function(current_simulation, scenario_data)
{
  current_data <- filter(current_simulation$trial_model$data, group_id == "trial")
  final_configuration <- tail(current_simulation$trial_model$trial_configuration, n = 1)
  distinct_doses <- unique(current_data$dose)
  
  MTD_dose <- final_configuration$MTD_dose
  
  MTD_prob_true <- c(NA, NA, NA)
  MTD_prob_true_cumulative <- NA
  
  if (!is.na(MTD_dose)) 
  {
    MTD_prob_true <- inner_join(tibble(cycle_index = 1:3), filter(scenario_data, dose == MTD_dose))$prob_true
    MTD_prob_true_cumulative <- MTD_prob_true[1] + 
      (1-MTD_prob_true[1])*(MTD_prob_true[2] + 
      (1-MTD_prob_true[2])*(MTD_prob_true[3]))
  } 
  
  
  result_tibble <- tibble(
    # num_patients_underdose = sum(filter(current_data, prob_true <= 0.16)$num_patients),
    # num_patients_target = sum(filter(current_data, prob_true > 0.16 & prob_true <= 0.33)$num_patients),
    # num_patients_overdose = sum(filter(current_data, prob_true > 0.33)$num_patients),
    MTD_dose = MTD_dose,
    MTD_prob_true_cycle1 = MTD_prob_true[1],
    MTD_prob_true_cycle2 = MTD_prob_true[2],
    MTD_prob_true_cycle3 = MTD_prob_true[3],
    MTD_prob_true_cumulative = MTD_prob_true_cumulative,
    MTD_is_underdose_cumulative = !is.na(MTD_dose) && (MTD_prob_true_cumulative <= 0.16),
    MTD_is_target_cumulative    = !is.na(MTD_dose) && (MTD_prob_true_cumulative > 0.16 & MTD_prob_true_cumulative <= 0.33),
    MTD_is_overdose_cumulative  = !is.na(MTD_dose) && (MTD_prob_true_cumulative > 0.33),
    num_toxicities_cumulative   = sum(current_data$num_toxicities),
    num_patients_start          = sum(current_data$num_patients_start),
    num_patients_cumulative     = sum(current_data$num_patients),
    num_patients_dropout        = sum(current_data$num_patients_start) - sum(current_data$num_patients),
    num_patients_underdose_cumulative = sum(filter(current_data, prob_true_cumulative <= 0.16)$num_patients),
    num_patients_target_cumulative = sum(filter(current_data, prob_true_cumulative > 0.16 & prob_true_cumulative <= 0.33)$num_patients),
    num_patients_overdose_cumulative = sum(filter(current_data, prob_true_cumulative > 0.33)$num_patients),
    last_cohort_time = max(c(current_data$cohort_time, 0)),
    stopped_due_to_toxicity = final_configuration$stopped_due_to_toxicity,
    stopped_due_to_max_N = final_configuration$stopped_due_to_max_N,
    last_patient_cens_event_dt = max(pmin(current_data$cens_caltime_dt, current_data$event_caltime_dt)),
    start_to_last_patient_cens_event_dt = last_patient_cens_event_dt - min(current_data$cohort_start_dt),
    mean_dose = mean(current_data$dose)
  ) 
  
  for (i in 1:3) {
    result_tibble[[paste0("MTD_is_underdose_cycle", i)]] <- !is.na(MTD_dose) && (MTD_prob_true[i] <= 0.16)
    result_tibble[[paste0("MTD_is_target_cycle", i)]]    <- !is.na(MTD_dose) && (MTD_prob_true[i] > 0.16 & MTD_prob_true[i] <= 0.33)
    result_tibble[[paste0("MTD_is_overdose_cycle", i)]]  <- !is.na(MTD_dose) && (MTD_prob_true[i] > 0.33)
    
    
    result_tibble[[paste0("num_patients_cycle", i)]] <- sum(current_data[[paste0("num_patients_", i, "cycle")]])
    result_tibble[[paste0("num_toxicities_cycle", i)]] <- sum(current_data[[paste0("num_toxicities_", i, "cycle")]])
    
    var <- paste0("prob_true_cycle", i)
    result_tibble[[paste0("num_patients_underdose_cycle", i)]] <- 
      sum(filter(current_data, get(var) <= 0.16)[[paste0("num_patients_", i, "cycle")]])
    result_tibble[[paste0("num_patients_target_cycle", i)]] <- 
      sum(filter(current_data, get(var) > 0.16 & get(var) <= 0.33)[[paste0("num_patients_", i, "cycle")]])
    result_tibble[[paste0("num_patients_overdose_cycle", i)]] <- 
      sum(filter(current_data, get(var) > 0.33)[[paste0("num_patients_", i, "cycle")]])
    
  }
  
  result_tibble$num_overdose_patient_cycles = sum(
    result_tibble$num_patients_overdose_cycle1,
    result_tibble$num_patients_overdose_cycle2,
    result_tibble$num_patients_overdose_cycle3
  )
  

  return(result_tibble)
}

reduce_results_per_scenario <- function(per_trial_results)
{
  n_trials_by_exp <- per_trial_results %>% 
    group_by(experiment) %>%
    summarize(n_trials_tmp = n())

  
  result_tibble <- per_trial_results %>% 
    group_by(experiment) %>%
    inner_join(n_trials_by_exp) %>% 
    dplyr::summarize(
      n_trials                             = unique(n_trials_tmp),
      total_num_patients_cumulative        = sum(num_patients_cumulative),
      total_num_toxicities_cumulative      = sum(num_toxicities_cumulative),
      mean_num_cohorts                     = sum(last_cohort_time) / n_trials,
      mean_num_patients_cumulative         = sum(num_patients_cumulative) / n_trials,
      mean_num_patients_start              = sum(num_patients_start) / n_trials,
      mean_num_patients_dropout            = sum(num_patients_dropout) / n_trials,
      mean_num_toxicities_cumulative       = sum(num_toxicities_cumulative) / n_trials,
      prob_patients_underdose_cumulative   = sum(num_patients_underdose_cumulative) / sum(num_patients_cumulative),
      prob_patients_target_cumulative      = sum(num_patients_target_cumulative)    / sum(num_patients_cumulative),
      prob_patients_overdose_cumulative    = sum(num_patients_overdose_cumulative)  / sum(num_patients_cumulative),
      mean_MTD                             = mean(MTD_dose, na.rm = T),
      prob_MTD_found                       = sum(!is.na(MTD_dose))            / n_trials,
      prob_MTD_underdose_cumulative        = sum(MTD_is_underdose_cumulative) / n_trials,
      prob_MTD_target_cumulative           = sum(MTD_is_target_cumulative)    / n_trials,
      prob_MTD_overdose_cumulative         = sum(MTD_is_overdose_cumulative)  / n_trials,
      prob_stopped_due_to_toxicity         = sum(stopped_due_to_toxicity)     / n_trials,
      prob_stopped_due_to_max_N            = sum(stopped_due_to_max_N)        / n_trials,
      mean_num_overdose_cycles             = sum(num_overdose_patient_cycles) / sum(num_patients_start),
      mean_last_patient_cens_event_dt = mean(last_patient_cens_event_dt),
      mean_start_to_last_patient_event_cens_dt = mean(start_to_last_patient_cens_event_dt),
      mean_dose                            = mean(mean_dose),
      .groups = "drop"
    ) 
  
  
  for (i in 1:3) {
    result_tibble_temp <- per_trial_results %>% 
      group_by(experiment) %>%
      inner_join(n_trials_by_exp) %>% 
      dplyr::summarize(
        n_trials                                         = unique(n_trials_tmp),
        !!paste0("mean_num_patients_cycle", i)           := sum(!!rlang::sym(paste0("num_patients_cycle", i))) / n_trials,
        !!paste0("mean_num_toxicities_cycle", i)         := sum(!!rlang::sym(paste0("num_toxicities_cycle", i))) / n_trials,
        !!paste0("prob_patients_underdose_cycle", i)     := sum(!!rlang::sym(paste0("num_patients_underdose_cycle", i))) / sum(!!rlang::sym(paste0("num_patients_cycle", i))),
        !!paste0("prob_patients_target_cycle", i)        := sum(!!rlang::sym(paste0("num_patients_target_cycle", i))) / sum(!!rlang::sym(paste0("num_patients_cycle", i))),
        !!paste0("prob_patients_overdose_cycle", i)      := sum(!!rlang::sym(paste0("num_patients_overdose_cycle", i))) /sum(!!rlang::sym(paste0("num_patients_cycle", i))),
        !!paste0("prob_MTD_underdose_cycle", i)          := sum(!!rlang::sym(paste0("MTD_is_underdose_cycle", i))) / n_trials,
        !!paste0("prob_MTD_target_cycle", i)             := sum(!!rlang::sym(paste0("MTD_is_target_cycle", i))) / n_trials,
        !!paste0("prob_MTD_overdose_cycle", i)           := sum(!!rlang::sym(paste0("MTD_is_overdose_cycle", i))) / n_trials
      )
    result_tibble <- result_tibble %>% inner_join(result_tibble_temp)
  }

  return(result_tibble)
}


reduce_results_to_MTD_overview <- function(per_trial_results)
{
  n_trials_by_exp <- per_trial_results %>% 
    group_by(experiment) %>%
    summarize(n_trials = n())

  result_tibble <- per_trial_results %>%
    group_by(experiment, MTD_dose) %>%
    dplyr::summarize(
      MTD_is_underdose = unique(MTD_is_underdose),
      MTD_is_target    = unique(MTD_is_target),
      MTD_is_overdose  = unique(MTD_is_overdose),
      n             = n(),
      .groups = "drop"
    ) %>% 
    inner_join(n_trials_by_exp) %>% 
    mutate(prob = n/n_trials)

  return(result_tibble)
}
