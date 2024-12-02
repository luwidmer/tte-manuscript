library("brms")

method_blrm_tte <- list(
  data = NULL,
  trial_configuration = NULL,
  tte_model = NULL,
  limit_toxicity = "per_cycle",
  last_analysis = ymd(NA),
  
  expand_scenario = function(
    scenario, 
    start_date = ymd("2010-01-10"),
    cycle_time = weeks(6)
  )
  {
    inf_date  <- as_date(Inf)
    inf_datetime  <- as_datetime(Inf)
    
    expanded_scenario <- mutate(
      scenario,
      cycle1_caltime_dt = start_date,
      cycle2_caltime_dt = start_date + cycle_time,
      cycle3_caltime_dt = start_date + 2*cycle_time,
      cens_caltime_dt   = start_date + cycle_index*cycle_time,
      num_cycles        = cycle_index,
      event_caltime_dt  = inf_date
    )
    
    return(expanded_scenario)
  },
  
  update = function(model, trial_configuration, cohort, scenario) 
  {
    # Add new cohort to data
    model$data <- dplyr::bind_rows(model$data, cohort)
    
    if (sum(filter(model$data, group_id == "trial")$num_patients_1cycle) >= trial_configuration$max_patients) {
      # We have to stop the trial here
      
      # Update model after all data is observed
      analysis_time <- max(filter(model$data, group_id == "trial")$cens_caltime_dt)
      
      model$last_analysis <- analysis_time
    } else {
      analysis_time <- max(pmin(filter(model$data, group_id == "trial")$cycle2_caltime_dt, filter(model$data, group_id == "trial")$cens_caltime_dt))
      
      model$last_analysis <- analysis_time
    }
    
    # current_data <- prepare_data(model$data, analysis_time)
    
    current_schedule_data <- model$expand_scenario(scenario)
    
    prepared_data <-  prepare_data(model$data, analysis_time)
    
    if (nrow(prepared_data) > 0)
    {
      # Not all patients dropped out before being eligible for cycle 1
    
      cycle_data_available <- length(unique(prepared_data$ocycle))
      
      
      if(cycle_data_available < 3) {
        # Patch missing cycles (brms model fit crashes with wrong simplex dimensions)
        prepared_data_2 <- prepared_data[c(seq_len(cycle_data_available), rep(1, 3 - cycle_data_available)), ] 
        
        prepared_data_2[(cycle_data_available+1):3, c("num_patients", "num_patients_2", "num_toxicities")] <- 0
        prepared_data_2[(cycle_data_available+1):3, c("follow_up")] <- 1e-4
        prepared_data_2[(cycle_data_available+1):3, ]$ocycle <- (cycle_data_available+1):3
        prepared_data_2[(cycle_data_available+1):3, ]$cycle <- (cycle_data_available+1):3
      } else {
        prepared_data_2 <- prepared_data
      }
    
      # prepared_data_end <- prepare_data(model$data)
      #model$combo2_tte_model_now <- update(model$combo2_tte_model, newdata=prepared_data, control=list(adapt_delta=0.95), cores=1, seed=356456, refresh=0, drop_unused_levels = FALSE)
      model$combo2_tte_model <- update(model$combo2_tte_model, newdata=prepared_data_2, control=list(adapt_delta=0.95), cores=1, seed=356456, refresh=0, drop_unused_levels = FALSE)
      #model$combo2_tte_model_end <- update(model$combo2_tte_model, newdata=prepared_data_end, control=list(adapt_delta=0.95), cores=1, seed=356456, refresh=0, drop_unused_levels = FALSE)
      
      
      model$post_crisk <- predict_crisk(model$combo2_tte_model, current_schedule_data)
      model$post_risk  <- predict_risk (model$combo2_tte_model, current_schedule_data)
      model$sum_crisk  <- summarise_risk(model$post_crisk) %>% mutate(ewoc_ok=q75 < 0.33)
      model$sum_risk   <- summarise_risk(model$post_risk)  %>% mutate(ewoc_ok=q75 < 0.33)
      
      model$dlt_crisk <- bind_cols(select(current_schedule_data, sid, drugA, drugB, dose, num_cycles), model$sum_crisk)
      model$dlt_risk <-  bind_cols(select(current_schedule_data, sid, drugA, drugB, dose, num_cycles), model$sum_risk)
      
      
      current_dose <- unique(cohort$dose)
      
  
      
      
      if (model$limit_toxicity == "per_cycle") {
        # We only choose the max. dose for cycle 1:
        # Make sure EWOC is okay over all cycles when choosing a dose to enroll at
        possible_doses_with_ewoc_in_all_cycles <- model$dlt_crisk %>% 
          group_by(dose) %>% 
          summarize(ewoc_ok_sum = sum(ewoc_ok), n_cycles = n()) %>%
          filter(ewoc_ok_sum == n_cycles)
        
        possible_doses <- dplyr::filter(model$dlt_crisk, dose %in% possible_doses_with_ewoc_in_all_cycles$dose, num_cycles == 1)
        current_target_probability <- max(filter(model$dlt_crisk, dose == current_dose)$prob_target)
      } else {
        # Limit cumulative toxicity over 3 cycles
        possible_doses <- model$dlt_risk %>% 
          filter(num_cycles == 3, ewoc_ok)
        
        current_target_probability <- filter(model$dlt_risk, dose == current_dose, num_cycles == 3)$prob_target
      }
      
      
      possible_doses <- dplyr::filter(possible_doses, dose <= model$max_escalation_factor*current_dose)
      
      if (nrow(possible_doses) > 0)
      {
        next_dose <- slice_max(possible_doses, drugB, n = 1)$drugB[1]
        
        assert_that(nrow(trial_configuration) == 1)
        assert_that(length(sum(model$data$num_patients_1cycle)) == 1)
        assert_that(length(trial_configuration$MTD_min_patients) == 1)
        assert_that(length(trial_configuration$MTD_min_prob_target) == 1)
        assert_that(length(current_target_probability) == 1) 
        
        
        if (
          (sum(dplyr::filter(model$data, dose == trial_configuration$current_dose)$num_patients_1cycle) >= trial_configuration$MTD_min_patients_on_dose &&
           next_dose == trial_configuration$current_dose && 
           (
             sum(model$data$num_patients_1cycle) >= trial_configuration$MTD_min_patients ||
             current_target_probability >= trial_configuration$MTD_min_prob_target
           )
          )
        )
        {
          # We can declare MTD 
          trial_configuration$MTD_dose <- trial_configuration$current_dose
          trial_configuration$current_dose <- NA
          trial_configuration$is_running <- FALSE
          model$combo2_tte_model <- NULL
          model$post_crisk <- NULL
        } else {
          if (sum(model$data$num_patients_1cycle) >= trial_configuration$max_patients) {
            # Can't declare MTD due to max N
            trial_configuration$current_dose <- NA
            trial_configuration$is_running <- FALSE
            trial_configuration$stopped_due_to_max_N <- TRUE
            
            
            # Check if at least an RDE can be declared (MTD criteria fulfilled, except for being the "maximum")
            data_by_dose <- model$data %>% 
              dplyr::group_by(dose) %>% 
              dplyr::summarize(
                num_patients_total = sum(num_patients_1cycle),
                .groups = "drop_last"
              )
            
            possible_RDE_doses <- filter(data_by_dose, num_patients_total >= trial_configuration$MTD_min_patients_on_dose)
            
            RDE_dose <- slice_max(
              filter(possible_doses, dose %in% possible_RDE_doses$dose),
              prob_target,
              n = 1
            )
            
            if (nrow(RDE_dose) == 1) {
              # If a lower dose with 6 patients and <= 1 toxicities exists, call that the RD at max N
              trial_configuration$RDE_dose <- RDE_dose$dose
            } 
            model$combo2_tte_model <- NULL
            model$post_crisk <- NULL
          } else {
            trial_configuration$current_dose <- next_dose
          }
        }
      } else {
        trial_configuration$current_dose <- NA
        trial_configuration$is_running <- FALSE
        trial_configuration$stopped_due_to_toxicity <- TRUE
        model$combo2_tte_model <- NULL
        model$post_crisk <- NULL
      }
    }
    
    trial_configuration$current_date <- model$last_analysis
    
    # Add decision to trial configurations
    model$trial_configuration <- dplyr::bind_rows(model$trial_configuration, trial_configuration)
    return(model)
  },
  
  get_next_configuration = function(model)
  {
    return(model$trial_configuration[nrow(model$trial_configuration), ])
  },
  
  set_initial_state = function(model, trial_configuration)
  {
    model$trial_configuration <- trial_configuration
    return(model)
  },
  
  setup = function(model, hist_data, drug_info, dose_info, combo2_tte_model, limit_toxicity = "per_cycle", max_escalation_factor = 2) {
    model$data <- hist_data
    
    model$combo2_tte_model <- combo2_tte_model
    model$max_escalation_factor <- max_escalation_factor
    model$limit_toxicity <- limit_toxicity
    
    return(model)
  }
)
