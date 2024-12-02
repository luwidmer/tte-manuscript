method_3cycle_blrm_trial_cie <- list(
  data = NULL,
  trial_configuration = NULL,
  trial = NULL,

  update = function(model, trial_configuration, cohort, scenario)
  {
    cohort_3cycle <- dplyr::select(cohort, -num_patients, -num_toxicities) %>%
      dplyr::mutate(
        num_patients = as.numeric(num_patients_3cycle),
        num_toxicities = as.numeric(num_toxicities_3cycle)
      )
    
    
    cohort_new <- cohort_3cycle %>% mutate(stratum_id="all")
    
    # Add new cohort to data
    model$data <- dplyr::bind_rows(model$data, cohort_new)

    cat("Added cohort:\n")
    print(cohort_new)
    print(summary(model$trial))
    
    
    model$trial <- update(model$trial, add_data = cohort_new)

    ##cat("Model data:\n")
    ##print(model$data)

    dose_prediction <- summary(model$trial, summarize="dose_prediction")
    ##cat("Dose prediction posterior:\n")
    ##print(dose_prediction)


    possible_doses <- dplyr::filter(dose_prediction, ewoc_ok)
    current_dose <- unique(cohort_new$drugB)
    possible_doses <- dplyr::filter(possible_doses, drugB <= model$max_escalation_factor*current_dose)
    
    
    current_dose_prediction <- filter(dose_prediction, drugB == current_dose)
    
    ##cat("Possible doses:\n")
    ##print(possible_doses)

    if (nrow(possible_doses) > 0)
    {
      next_dose <- slice_max(possible_doses, drugB, n = 1)$drugB[1]
      
      if (
        (sum(dplyr::filter(model$data, dose == trial_configuration$current_dose)$num_patients) >= trial_configuration$MTD_min_patients_on_dose &&
         next_dose == trial_configuration$current_dose && 
         (
           sum(model$data$num_patients) >= trial_configuration$MTD_min_patients ||
           current_dose_prediction$prob_target >= trial_configuration$MTD_min_prob_target
         )
        )
      )
      {
        # We can declare MTD 
        trial_configuration$MTD_dose <- trial_configuration$current_dose
        trial_configuration$current_dose <- NA
        trial_configuration$is_running <- FALSE
        model$trial <- NULL
      } else {
        if (sum(model$data$num_patients) >= trial_configuration$max_patients) {
          # Can't declare MTD due to max N
          trial_configuration$current_dose <- NA
          trial_configuration$is_running <- FALSE
          trial_configuration$stopped_due_to_max_N <- TRUE
          
          
          # Check if at least an RDE can be declared (MTD criteria fulfilled, except for being the "maximum")
          data_by_dose <- model$data %>% 
            dplyr::group_by(dose) %>% 
            dplyr::summarize(
              num_patients_total = sum(num_patients), 
              num_toxicities_total = sum(num_toxicities),
              .groups = "drop_last"
            )
          
          possible_RDE_doses <- filter(data_by_dose, num_patients_total >= trial_configuration$MTD_min_patients_on_dose)
          
          RDE_dose <- slice_max(
            filter(possible_doses, drugB %in% possible_RDE_doses$dose),
            prob_target,
            n = 1
          )
          
          if (nrow(RDE_dose) == 1) {
            # If a lower dose with 6 patients and <= 1 toxicities exists, call that the RD at max N
            trial_configuration$RDE_dose <- RDE_dose$drugB
          } 
          model$trial <- NULL
        } else {
          trial_configuration$current_dose <- next_dose
        }
      }
    } else {
      trial_configuration$current_dose <- NA
      trial_configuration$is_running <- FALSE
      trial_configuration$stopped_due_to_toxicity <- TRUE
      model$trial <- NULL
    }
    
    trial_configuration$current_date <- max(model$data$cens_caltime_dt)

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

  setup = function(
    model, 
    hist_data, 
    drug_info, 
    dose_info,
    ## trial,
    max_escalation_factor = 2
  ) 
  {
    model$max_escalation_factor <- max_escalation_factor
    
    prob_3_cycle <- \(p) {p + (1-p)*p + (1-p)^2*p}
    
    if(TRUE) {
      
      
      model$trial <- blrm_trial(
        data = hist_data,
        drug_info = drug_info,
        dose_info = dose_info,
        formula_generator=function(...) blrm_formula_saturating(..., max_interaction_level=1)
        # ,interval_prob = prob_3_cycle(c(0, 0.16, 0.33, 1))
      )
      
      dims <- summary(model$trial, "dimensionality")
    
      model$trial <- update(
        model$trial,
    
        prior_EX_mu_mean_comp = matrix(
          c(
            logit(0.09), log(1.0), # hyper-mean of (intercept, log-slope) for drug A
            logit(0.12), log(1.0)  # hyper-mean of (intercept, log-slope) for SOC
          ),
          nrow = dims$num_components,
          ncol = 2,
          byrow = TRUE
        ),
        prior_EX_mu_sd_comp = matrix(
          c(
            1.3, 1.2, # sd of (intercept, log-slope) for drug A
            0.7, log(4)/1.96  # sd of (intercept, log-slope) for SOC
          ), 
          # c(
          #   exp(0.22*1.25), 1.2,#log(4)/2, # sd of (intercept, log-slope) for drug A
          #   exp(0.22*1.25), log(4)/2  # sd of (intercept, log-slope) for SOC
          # ), 
          nrow = dims$num_components,
          ncol = 2,
          byrow = TRUE
        ),
        prior_EX_mu_mean_inter = rep(0, dims$num_interaction_terms),
        prior_EX_mu_sd_inter = rep(1.121, dims$num_interaction_terms),
        ## Here we take tau as known and as zero.
        ## This disables the hierarchical prior which is
        ## not required in this example as we analyze a
        ## single trial.
        prior_EX_tau_mean_comp = matrix(
          0,
          nrow = dims$num_components,
          ncol = 2,
          byrow = TRUE
        ),
        prior_EX_tau_sd_comp = matrix(
          1,
          nrow = dims$num_components,
          ncol = 2,
          byrow = TRUE
        ),
        prior_EX_tau_mean_inter = matrix(
          0, 
          nrow = dims$num_strata, 
          ncol = dims$num_interaction_terms
        ),
        prior_EX_tau_sd_inter =  matrix(
          1, 
          nrow = dims$num_strata, 
          ncol = dims$num_interaction_terms
        ),
        prior_is_EXNEX_comp  = rep(FALSE, dims$num_components),
        prior_EX_prob_comp   = matrix(1, nrow = dims$num_groups, ncol = dims$num_components),
        prior_is_EXNEX_inter = rep(FALSE, dims$num_interaction_terms),
        prior_EX_prob_inter  = matrix(1, nrow = dims$num_groups, ncol = dims$num_interaction_terms),
        prior_tau_dist = 0,
        prior_PD = FALSE,
        iter = 2000,
        warmup = 1000,
        chains = 4
      )
    }
    
    # model$trial <- trial
    
    dose_prediction <- summary(model$trial, summarize="dose_prediction")
    
    # assert_that(sum(dose_prediction$ewoc_ok) != 0, msg="No dose is safe at trial start!")
    
    return(model)
  }
)
