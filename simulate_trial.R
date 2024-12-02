simulate_trial_clustermq <- function(experiment, experiment_args, job_id, replication_index)
{
  # simulate_error() # Simulate an error
  
  scenario <- experiment_args$scenario
  population <- experiment_args$fun.sample_patients_for_scenario(scenario) #, max_patients_per_dose = max_patients_per_dose_overall)
  
  fun.generate_cohort <- experiment_args$fun.generate_cohort(population)
  
  trial_configuration <- experiment_args$trial_configuration
  
  model <- experiment_args$model
  
  
  # The result is expected to be a list
  computed_output <- 
    simulate_trial(
      scenario,
      trial_configuration,
      fun.generate_cohort,
      population,
      model
    )
  

  result <- list()
  if (experiment_args$keep_full_results) {
    result$computed_output <- computed_output
  }
  
  # Pre-computing a tibble with stats for each simulation study 
  # on the workers allows for rapid result aggregation on the master node
  result$trial_tibble <- bind_cols(
    "job_id" = job_id,
    "computed_output" = reduce_results_per_trial(computed_output, scenario)
  )
  
  result$dose_tibble <- bind_cols(
    "job_id" = job_id,
    "computed_output" = reduce_results_per_dose(computed_output, scenario)
  )
  
  return(result)
  
}

simulate_error <- function()
{
  stop("This is a simulated error")
}