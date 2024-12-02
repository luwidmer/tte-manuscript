set.seed(90234198)



source("cluster_engine.R")
source("load_packages.R")
source("load_functions.R")
source("plot_functions.R")

source("combo2-tte-model.R")

source("combo2-tte-oc-methods-1cycle-blrm-trial.R")
source("combo2-tte-oc-methods-3cycle-blrm-trial-cie.R")
source("combo2-tte-oc-methods-blrm-tte.R")

time_stamp <- gsub(Sys.time(), pattern = "[ :]", replacement = "-")
path.export <- file.path(normalizePath("results"), time_stamp)
dir.create(path.export)
theme_set(theme_bw())


lambda_scenario <- list()

doses <- c(10, 20, 40, 80, 160, 320, 640, 1280)
n_cycles <- 3
n_doses <- length(doses)
doses_by_cycle <- unlist(lapply(doses, \(dose) {rep(dose, n_cycles)}))


lambda_scenario[["constant"]] <- tibble( 
  drugA = 1,
  drugB = doses_by_cycle, 
  sid = unlist(lapply(1:n_doses, function(x){rep(x, n_cycles)})),
  cycle_index = rep(seq_len(n_cycles), n_doses ), 
  prob_true = c(
    rep(0.05, n_cycles),
    rep(0.06, n_cycles),
    rep(0.07, n_cycles),
    rep(0.09, n_cycles),    
    rep(0.11, n_cycles), 
    rep(0.21, n_cycles),
    rep(0.35, n_cycles),
    rep(0.47, n_cycles)
  ),
  dropout_rate = 0,
  dropout_total = 0,
  dropout_type = "none",
  group_id = as_group_factor_nonhierarchical("trial")
) %>% mutate(
  dose = drugB,
  lambda = exp(cloglog(prob_true) - log(cycle_length_d))
)

increasing_risk_logit_offset <- c(-1.3, 0, 0.6)
lambda_scenario[["increasing"]] <- lambda_scenario[["constant"]] %>% 
  mutate(
    prob_true = inv_logit(logit(prob_true) + increasing_risk_logit_offset),
  ) %>% mutate(
    lambda = exp(cloglog(prob_true) - log(cycle_length_d))
  )


decreasing_risk_logit_offset <- rev(increasing_risk_logit_offset)

lambda_scenario[["decreasing"]] <- lambda_scenario[["constant"]] %>% 
  mutate(
    prob_true = inv_logit(logit(prob_true) + decreasing_risk_logit_offset),
  ) %>% mutate(
    lambda = exp(cloglog(prob_true) - log(cycle_length_d))
  )


lambda_scenario_noninfo <- lambda_scenario
print(plot_scenarios(lambda_scenario_noninfo))

print(plot_scenarios_cumulative(lambda_scenario_noninfo))

ggsave(file.path(path.export, "scenarios_conditional.pdf"), plot = plot_scenarios(lambda_scenario_noninfo)           , width=20*1.2, height=10, units="cm")
ggsave(file.path(path.export, "scenarios_cumulative.pdf"),  plot = plot_scenarios_cumulative(lambda_scenario_noninfo), width=20*1.2, height=10, units="cm")
ggsave(file.path(path.export, "scenarios_conditional.png"), plot = plot_scenarios(lambda_scenario_noninfo)           , width=20*1.2, height=10, units="cm")
ggsave(file.path(path.export, "scenarios_cumulative.png"),  plot = plot_scenarios_cumulative(lambda_scenario_noninfo), width=20*1.2, height=10, units="cm")

dropout_rate <- \(dropout_fraction) {exp(cloglog(dropout_fraction) - log(3*cycle_length_d))}
dropout_fraction <- \(dropout_rate) {inv_cloglog(log(dropout_rate) + log(3*cycle_length_d))}


scenario_names <- names(lambda_scenario)
for (scenario_name in scenario_names)
{

  # Scenario with 2*dropout rate of 0.33
  lambda_scenario[[paste0(scenario_name, "_5511")]] <- mutate(
    lambda_scenario[[scenario_name]], 
    dropout_total = 0.5511, 
    dropout_rate = dropout_rate(dropout_total),  
    dropout_type = "constant"
  )
  lambda_scenario[[paste0(scenario_name, "_33")]] <- mutate(
    lambda_scenario[[scenario_name]], 
    dropout_total = 0.33, 
    dropout_rate = dropout_rate(dropout_total),
    dropout_type = "constant"
  )
}

scenario_name <- "constant_33"  
lambda_scenario[[paste0(scenario_name, "_info_incr_33")]] <- mutate(
  lambda_scenario[[scenario_name]], 
  # dropout_rate = (drugB - 10) * 1.682e-05,
  dropout_rate = if_else(sid > 4, dropout_rate * 2, 0),
  dropout_type = "increasing"
)

lambda_scenario[[paste0(scenario_name, "_info_decr_33")]] <- mutate(
  lambda_scenario[[scenario_name]], 
  # dropout_rate = (1280 - drugB) * 3.092e-06, 
  dropout_rate = if_else(sid <= 4, dropout_rate * 2, 0),
  dropout_type = "decreasing"
)


fun.create_patient_generator_time_poisson_piecewise_reldate <- function(
    scenario_data, 
    cycle_length = weeks(6),
    accrual_rate = 1/10 # 1 patient in 10 days
) 
{
  
  patient_id <- 0
  cycle_length 
  cycle_length_d <- as.numeric(cycle_length, "days")
  epochs <- cumsum(rep(cycle_length_d, 2)) # Cycle lengths in day
  scenario_data
  accrual_rate
  
  enrolment_doses <- unique(scenario_data$dose)
  
  function(
    n_patients_per_dose,
    patient_doses
  ) {
    assert_that(n_patients_per_dose == length(patient_doses))
    
    patients <- NULL
    
    patient_id_local <- 0
    patient_start_offset <- days(0)
    for (current_patient in seq_len(n_patients_per_dose)) {
      patient_id_local <- patient_id_local + 1
      patient_id <- patient_id + 1
      
      patient_starting_dose <- patient_doses[patient_id_local]
      
      
      patient_cycles <- dplyr::filter(scenario_data, dose==patient_starting_dose)
      
      patient_start_offset <- patient_start_offset + days(ceiling(rexp(1, rate = accrual_rate)))
      
      event_time <- ceiling(rpexp(lambda = patient_cycles$lambda, epochs = epochs))
      
      
      cycle2_offset     <- patient_start_offset + cycle_length
      cycle3_offset     <- patient_start_offset + 2*cycle_length
      end_of_last_cycle_offset <- patient_start_offset + 3*cycle_length
      
      eligible_1cycle_offset <- cycle2_offset
      eligible_2cycle_offset <- cycle3_offset
      eligible_3cycle_offset <- end_of_last_cycle_offset
      
      event_offset <- patient_start_offset + days(event_time)
      
      if (all(patient_cycles$dropout_rate == 0))
      {
        dropout_time <- Inf
        dropout_offset <- period(Inf, "seconds")
        cens_offset  <- end_of_last_cycle_offset
      } else {
        dropout_time   <- ceiling(rpexp(lambda = patient_cycles$dropout_rate, epochs = epochs))
        dropout_offset <- patient_start_offset + days(dropout_time)
        cens_offset  <- pmin(end_of_last_cycle_offset, dropout_offset)
      }
      
      if (cens_offset == 0)
      {
        browser()
        stop("Censoring should not happen before day 1")
      }
      
      current_patient <- tibble(
        dose              = patient_starting_dose, 
        patient_id        = patient_id, 
        patient_seen_in_cohort_index = 0,
        toxicity_time_d   = event_time,
        cycle1_offset     = patient_start_offset,
        cycle2_offset     = cycle2_offset,
        cycle3_offset     = cycle3_offset,
        end_of_last_cycle_offset = end_of_last_cycle_offset,
        eligible_1cycle_offset = eligible_1cycle_offset,
        eligible_2cycle_offset = eligible_2cycle_offset,
        eligible_3cycle_offset = eligible_3cycle_offset,
        dropout_rate_cycle1 = patient_cycles$dropout_rate[1],
        dropout_rate_cycle2 = patient_cycles$dropout_rate[2],
        dropout_rate_cycle3 = patient_cycles$dropout_rate[3],
        dropout_offset = dropout_offset,
        cens_offset    = cens_offset,
        event_offset  = event_offset,
        prob_true_cycle1  = patient_cycles$prob_true[1],
        prob_true_cycle2  = patient_cycles$prob_true[2],
        prob_true_cycle3  = patient_cycles$prob_true[3],
        prob_true_cumulative = patient_cycles$prob_true[1] + 
          (1 - patient_cycles$prob_true[1]) * (patient_cycles$prob_true[2] + 
          (1 - patient_cycles$prob_true[2]) *  patient_cycles$prob_true[3]), 
        num_patients       = 1,
        num_patients_start = 1,
        num_toxicities        = event_offset <= cens_offset,
        num_toxicities_1cycle = event_offset <= pmin(eligible_1cycle_offset, cens_offset),
        num_toxicities_2cycle = event_offset <= pmin(eligible_2cycle_offset, cens_offset),
        num_toxicities_3cycle = event_offset <= pmin(eligible_3cycle_offset, cens_offset),
        num_patients_1cycle   = (eligible_1cycle_offset <= cens_offset) | num_toxicities_1cycle,
        num_patients_2cycle   = (eligible_2cycle_offset <= cens_offset) | num_toxicities_2cycle,
        num_patients_3cycle   = (eligible_3cycle_offset <= cens_offset) | num_toxicities_3cycle, 
        group_id = as_group_factor_nonhierarchical("trial")
      ) %>% mutate(drugA = 1, drugB = dose)
      
      patients <- dplyr::bind_rows(patients, current_patient)
    }
    
    patient_id <<- patient_id
    
    return(patients)
  }
}


fun.generate_cohort_TTE_reldate <- function(
    population_generator, 
    fun.cohort_num_patients = function(cohort_iteration_index, cohort) {3}
)
{
  population_generator
  fun.cohort_num_patients
  
  function(cohort_iteration_index,
           cohort,
           cohort_start_date
  )
  {
    dose_to_recruit <- unique(cohort$dose)
    assert_that(length(dose_to_recruit) == 1)
    
    num_patients_in_cohort <- fun.cohort_num_patients(cohort_iteration_index, cohort)
    
    patients_sampled_for_cohort <- population_generator(num_patients_in_cohort, rep(dose_to_recruit, num_patients_in_cohort))
    
    patients_sampled_for_cohort <- patients_sampled_for_cohort %>% mutate(
      cohort_start_dt   = cohort_start_date,
      patient_start_dt  = cohort_start_date + cycle1_offset,
      cycle1_caltime_dt = cohort_start_date + cycle1_offset,
      cycle2_caltime_dt = cohort_start_date + cycle2_offset,
      cycle3_caltime_dt = cohort_start_date + cycle3_offset,
      end_of_last_cycle_caltime_dt = cohort_start_date + end_of_last_cycle_offset,
      eligible_caltime_1cycle_dt = cohort_start_date + eligible_1cycle_offset,
      eligible_caltime_2cycle_dt = cohort_start_date + eligible_2cycle_offset,
      eligible_caltime_3cycle_dt = cohort_start_date + eligible_3cycle_offset,
      cens_caltime_dt   = cohort_start_date + cens_offset,
      event_caltime_dt  = cohort_start_date + event_offset
    )
    
    infinite_indices <- is.infinite(patients_sampled_for_cohort$dropout_offset)
    
    patients_sampled_for_cohort$dropout_caltime_dt <- as.Date(NA)
    patients_sampled_for_cohort$dropout_caltime_dt[infinite_indices]  <- as.Date(Inf)
    patients_sampled_for_cohort$dropout_caltime_dt[!infinite_indices] <- 
      patients_sampled_for_cohort$cohort_start_dt[!infinite_indices] + 
      patients_sampled_for_cohort$dropout_offset[!infinite_indices]
    
    assert_that(class(patients_sampled_for_cohort$dropout_caltime_dt) == "Date")

    if (nrow(patients_sampled_for_cohort) != num_patients_in_cohort) {
      stop("Insufficient patients in population!")
    }
    
    # cohort$num_patients   <- num_patients_in_cohort
    # # cohort$num_toxicities <- sum(patients_sampled_for_cohort$toxicity_event)
    patients_sampled_for_cohort$cohort_time    <- cohort_iteration_index
    
    patients_sampled_for_cohort
  }
}


patient_generator_reldate <- fun.create_patient_generator_time_poisson_piecewise_reldate(lambda_scenario$increasing)

cohort_generator_TTE_reldate <- fun.generate_cohort_TTE_reldate(patient_generator_reldate)

#
cohort_reldate <- NULL
for (i in 1:1000)
{
  cohort_reldate <- bind_rows(
    cohort_reldate,
    cohort_generator_TTE_reldate(0, tibble(num_patients = 0, dose = 160), cohort_start_date = ymd("2022-02-01"))
  )
}
mean(cohort_reldate$num_toxicities_1cycle)
mean(cohort_reldate$num_toxicities_2cycle)
mean(cohort_reldate$num_toxicities_3cycle)

prob_over_3_cycles <- \(p) {p[1] + (1-p[1])*(p[2] + (1-p[2])*p[3])}

prob_over_3_cycles(filter(lambda_scenario$decreasing, drugB == 160)$prob_true)



cohort_reldate <- cohort_generator_TTE_reldate(0, tibble(dose = 160), cohort_start_date = ymd("2022-02-01"))

full_patient <- cohort_reldate[1, ]
full_patient_prepped <- prepare_data(full_patient, full_patient$cens_caltime_dt)

assert_that(nrow(full_patient_prepped) == 3)
assert_that(all(full_patient_prepped$follow_up == 42))

cycle1_only_patient <- prepare_data(
  full_patient %>% mutate(cens_caltime_dt = cycle2_caltime_dt), 
  full_patient$cens_caltime_dt
)

assert_that(nrow(cycle1_only_patient) == 1)
assert_that(all(cycle1_only_patient$follow_up == 42))

cycle1_half_patient_result <- prepare_data(
  full_patient %>% mutate(cens_caltime_dt = cycle1_caltime_dt + weeks(5) + days(6)) %>% select(-contains("offset")), 
  full_patient$cens_caltime_dt
)

assert_that(nrow(cycle1_half_patient_result) == 0)

cycle2_half_patient_result <- prepare_data(
  full_patient %>% 
    mutate(cens_caltime_dt = cycle2_caltime_dt + weeks(5) + days(6)) %>% 
    select(-contains("offset")), 
  full_patient$cens_caltime_dt
)

assert_that(nrow(cycle2_half_patient_result) == 1)
assert_that(all(cycle2_half_patient_result$follow_up == 42))

cycle1_only_patient_partial_result <- prepare_data(
  full_patient, 
  analysis_caltime_dt = full_patient$cycle2_caltime_dt - days(1)
)

assert_that(nrow(cycle1_only_patient_partial_result) == 0)

dlt_day1_full_patient <- full_patient %>% mutate(
  event_caltime_dt = cycle1_caltime_dt + days(1)
)

dlt_day1_full_patient_result <- prepare_data(
  dlt_day1_full_patient,
  dlt_day1_full_patient$cens_caltime_dt
)

assert_that(nrow(dlt_day1_full_patient_result) == 1)
assert_that(dlt_day1_full_patient_result$num_toxicities == 1)
assert_that(all(dlt_day1_full_patient_result$follow_up == 42))

dropout_in_cycle1_patient <- full_patient %>% mutate(
  dropout_caltime_dt = cycle1_caltime_dt + days(1),
  cens_caltime_dt = dropout_caltime_dt
)

dropout_in_cycle1_patient_result <- prepare_data(
  dropout_in_cycle1_patient,
  dropout_in_cycle1_patient$cens_caltime_dt
)

assert_that(nrow(dropout_in_cycle1_patient_result) == 0)

# cohort
# cohort_reldate

analysis_time <- ymd("2025-03-01")
#prepare_data(cohort,  max(cohort$cycle2_caltime_dt))




# Run ---------------------------------------------------------------------



trial_configuration_TTE <- tibble(
  group_id                 = "trial",
  is_running               = TRUE,
  current_date             = ymd("2022-02-01"),
  current_dose             = 20, # starting dose
  MTD_dose                 = NA,
  max_patients             = 60, #30,
  MTD_min_patients         = 21, #12,
  MTD_min_prob_target      = 0.5,
  MTD_min_patients_on_dose = 6, #6,
  stopped_due_to_toxicity  = FALSE,
  stopped_due_to_max_N     = FALSE
)

## Uncomment the section below to generate fake patients and write these data
## into RDS files
for (i in 1:10) {
  fake_patient_generator <- fun.create_patient_generator_time_poisson_piecewise_reldate(lambda_scenario$constant_33)
  fake_patient_generator_reldate <- fun.generate_cohort_TTE_reldate(fake_patient_generator)

  fake_population <- NULL
  fake_population_cur_date <- ymd("2022-01-01")
  for (j in 1:3) {
    fake_population <- bind_rows(
      fake_population,
      fake_patient_generator_reldate(j, tibble(dose = doses[j + 3]), fake_population_cur_date)
    )
    fake_population_cur_date <- max(pmin(fake_population$cycle2_caltime_dt, fake_population$cens_caltime_dt))
  }
  
  fake_population <- fake_population |> inner_join(tibble(dose = doses[4:6], dose_shifted = doses[2:4]))
  saveRDS(fake_population, paste0("fake/fake_patients-", i,".RDS"))
}

## doses which we will be testing in our trial
dose_info <- tibble(
  stratum_id="all",
  group_id = as_group_factor_nonhierarchical("trial"),
  drugB = doses,
  drugA = 1
)

drug_info <- tibble(
  drug_name = c("drugB", "drugA"),
  dose_ref = c(drugB_dose_ref, drugA_dose_ref),
  dose_unit = c("mg", "a.u.")
)

# We will compare 4 models:
# - TTE that limits cumulative toxicity over 3 cycles
# - TTE that limits conditional toxicity in each of the 3 cycles given patients make it that far
# - 3-cycle BLRM with prior aligned to cumulative TTE model (over 3 cycles)
# - 1-cycle BLRM with prior aligned to cumulative/conditional TTE model (over 1st cycle)


TTE_model_per_cycle <- method_blrm_tte$setup(
  method_blrm_tte, 
  hist_data = hist_data_nonhierarchical_fake[0, ], 
  drug_info, 
  dose_info, 
  combo2_tte_model_nonhierarchical,
  limit_toxicity = "per_cycle"
) 

TTE_model_cumulative <- method_blrm_tte$setup(
  method_blrm_tte, 
  hist_data = hist_data_nonhierarchical_fake[0, ], 
  drug_info, 
  dose_info, 
  combo2_tte_model_nonhierarchical,
  limit_toxicity = "cumulative"
) 

BLRM_1cycle_model <- method_1cycle_blrm_trial$setup(
  method_1cycle_blrm_trial, 
  hist_data = NULL, #blrm_empty_hist_data,
  drug_info, 
  dose_info
) 

BLRM_3cycle_model <- method_3cycle_blrm_trial_cie$setup(
  method_3cycle_blrm_trial_cie, 
  hist_data = NULL, #blrm_empty_hist_data,
  drug_info, 
  dose_info
) 



tte_schedules <- TTE_model_cumulative$expand_scenario(lambda_scenario$constant)
tte_schedules <- tte_schedules %>% 
  bind_rows(tte_schedules %>% mutate(drugA = 0, sid = sid + max(sid))) %>%
  bind_rows(
    unique(tte_schedules %>% mutate(drugB = 0, dose = 0, sid = 0, prob_true = NA, lambda = NA))
  )



tte_prior_fit <- blrm_combo2_tte_fit_nonhierarchical(combo2_tte_model_nonhierarchical, tibble(), ymd("2025-01-10"))

ylims <- c(0, 1)

p_tte_prior <- plot_risk(tte_prior_fit, tte_schedules) +
  coord_cartesian(ylim = ylims) + 
  ggtitle("Overall risk for a DLT, TTE")


p_blrm_1cycle_3cycle <- plot_risk_ob2(
  BLRM_1cycle_model$trial, 
  BLRM_3cycle_model$trial, 
  mutate(tte_schedules, stratum_id="all")
) + coord_cartesian(ylim=ylims)

prior_plot <- arrangeGrob(p_blrm_1cycle_3cycle, p_tte_prior, ncol = 1)
grid.arrange(prior_plot)

ggsave(file.path(path.export, "prior_overview.png"), plot = prior_plot, width=20*1.5, height=10*1.5, units="cm")
ggsave(file.path(path.export, "prior_overview.pdf"), plot = prior_plot, width=20*1.5, height=10*1.5, units="cm")



p_max_blrm_1cycle <- plot_max_dose_ob2(BLRM_1cycle_model$trial, dose_col = "drugB")
p_max_blrm_3cycle <- plot_max_dose_ob2(BLRM_3cycle_model$trial, dose_col = "drugB", cycle_num = 3)




tite_plot  <-  plot_max_dose_tte(tte_prior_fit, "drugB", tte_schedules, drug_info)
tite_cplot <- plot_max_cdose_tte(tte_prior_fit, "drugB", tte_schedules, drug_info)


select_grobs <- function(lay) {
  id <- unique(c(t(lay))) 
  id[!is.na(id)]
} 


hlay <- rbind(c(1,NA,2),
              c(3,4,5))

gs <- c(
  list(p_max_blrm_1cycle, p_max_blrm_3cycle),
  tite_plot
)


prior_plot_continuous <- arrangeGrob(grobs=gs[select_grobs(hlay)], layout_matrix=hlay)
grid.arrange(prior_plot_continuous)

ggsave(
  file.path(path.export, "prior_overview_continuous.png"), 
  plot = prior_plot_continuous, 
  width = 20*2, 
  height = 10*2, 
  units = "cm"
)

ggsave(
  file.path(path.export, "prior_overview_continuous.pdf"), 
  plot = prior_plot_continuous, 
  width = 20*2, 
  height = 10*2, 
  units = "cm"
)
# browser()

hlay_c <- rbind(c(1,NA,NA),
              c(2,3,4))


gs_c <- c(
  list(p_max_blrm_1cycle),
  tite_cplot
)

prior_cplot_continuous <- arrangeGrob(grobs=gs_c[select_grobs(hlay_c)], layout_matrix=hlay_c)
grid.arrange(prior_cplot_continuous)

ggsave(file.path(path.export, "prior_cond_overview_continuous.png"), plot = prior_cplot_continuous, width=20*2, height=10*2, units="cm")
ggsave(file.path(path.export, "prior_cond_overview_continuous.pdf"), plot = prior_cplot_continuous, width=20*2, height=10*2, units="cm")


dropout_plot <- plot_dropout(lambda_scenario)
ggsave(file.path(path.export, "dropout_scenarios.png"), plot = dropout_plot, width=20*1.5, height=10, units="cm")
ggsave(file.path(path.export, "dropout_scenarios.pdf"), plot = dropout_plot, width=20*1.5, height=10, units="cm")

# print_prior_dose_comparison(summary_BLRM_1cycle, 1)
# print_prior_dose_comparison(summary_BLRM_3cycle, 3)

experiments <- list()
experiment_mapping <- NULL

  
models <- list(
  "blrm_1cycle" = BLRM_1cycle_model,
  "blrm_3cycle_cie" = BLRM_3cycle_model, 
  "tte_per_cycle" = TTE_model_per_cycle,
  "tte_cumulative" = TTE_model_cumulative
)

for (model_name in names(models))
{
  for (scenario_name in names(lambda_scenario)) #  "constant_33") #
  {
    max_patients <- trial_configuration_TTE$max_patients
    exp_name <- paste0(model_name, "_", scenario_name, "_", max_patients, "_patients")
    
    experiments[[exp_name]] <- list(
      model = models[[model_name]],
      seed = 123,
      keep_full_results = FALSE,
      trial_configuration = trial_configuration_TTE,
      fun.generate_cohort = fun.generate_cohort_TTE_reldate,
      fun.sample_patients_for_scenario = fun.create_patient_generator_time_poisson_piecewise_reldate,
      scenario = lambda_scenario[[scenario_name]]
    )

  }
}

if(T)
{
  ## Macbook M3 Max compatible settings
  # replications <- 50; n_cores <- 14 
  # options(clustermq.scheduler="multiprocess")
  
  ## High performance compute cluster settings
  replications <- 1000; n_cores <- 200 
  
  qsave(
    list(
      experiments = experiments, 
      lambda_scenario = lambda_scenario,
      cycle_length = cycle_length,
      dose_info = dose_info,
      replications = replications
    ), 
    file = file.path(path.export, paste0("inputs.qs"))
  )
  
  
  gc()
  # options(error = recover)
  result <- run_batch(
    experiments = experiments,
    batch_replications = replications,
    job_function = simulate_trial_clustermq,
    worker_setup_function = default_worker_setup(),
    n_jobs = 14
    # , test_single_job_through_clustermq = FALSE
    # , test_single_job_per_experiment = TRUE
    # , test_single_job_index = c(1004 ) # Use this option to debug single jobs that fail
  )
  
  
  gc()
  
  
  qsave(
    result$outcome_overall, 
    file = file.path(path.export, "result-overall.qs")
  )
  
  qsave(
    result, 
    file = file.path(path.export, "result.qs")
  )
  
  
  library(lobstr)
  print(obj_size(result))
  
}

print(sessionInfo())


