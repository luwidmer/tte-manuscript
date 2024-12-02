cloglog <- function(p) {
  binomial(cloglog)$linkfun(p)
}

inv_cloglog <- function(cl) {
  binomial(cloglog)$linkinv(cl)
}

## avoids underflow to 0
inv_log <- function(l) {
    poisson()$linkinv(l)
}


## define Stan functions once more as R variants
log_combo2_inv_link <- function(mu, phi, fup, finiteA, finiteB) {
     ## log(lambda)
     log_lambda <- log(fup)
     N <- length(mu)
     log_mu_phi <- apply(matrix(cbind(mu, phi), N, 2), 1, matrixStats::logSumExp)
     log_lambda <- log_lambda +
         case_when(finiteA == 1 & finiteB == 1 ~ log_mu_phi,
                   finiteA == 1 & finiteB != 1 ~ mu,
                   finiteA != 1 & finiteB == 1 ~ phi,
                   finiteA != 1 & finiteB != 1 ~ rep(-Inf, times=N)
                   )
     return(log_lambda)
 }


blrm_combo2_tte_lpmf <- function(y, mu, phi, fup, finiteA, finiteB) {
    theta <- log_combo2_inv_link(mu, phi, fup, finiteA, finiteB)
    #// poisson density without normalization y!
    ##return y * theta - exp(theta);
    return(sum(y * theta - exp(theta)))
}

blrm_combo2_tte_rng <- function(mu, phi, fup, finiteA, finiteB) {
    N <- length(mu)
    return(rpois(N, exp(log_combo2_inv_link(mu, phi, fup, finiteA, finiteB))))
}

log_lik_blrm_combo2_tte <- function(i, prep) {
  mu  <- brms::get_dpar(prep, "mu", i= i)
  phi  <- brms::get_dpar(prep, "phi", i= i)
  fup <- prep$data$vreal1[i]
  finiteA <- prep$data$vint1[i]
  finiteB <- prep$data$vint2[i]
  ##trials <- draws$data$trials[i]
  y <- prep$data$Y[i]
  blrm_combo2_tte_lpmf(y, mu, phi, fup, finiteA, finiteB)
}

posterior_predict_blrm_combo2_tte <- function(i, prep, ...) {
  mu  <- brms::get_dpar(prep, "mu", i= i)
  phi  <- brms::get_dpar(prep, "phi", i= i)
  fup <- prep$data$vreal1[i]
  finiteA <- prep$data$vint1[i]
  finiteB <- prep$data$vint2[i]
  ##trials <- draws$data$trials[i]
  blrm_combo2_tte_rng(mu, phi, fup, finiteA, finiteB)
}


posterior_epred_blrm_combo2_tte <- function(prep) {
  mu  <- brms::get_dpar(prep, "mu")
  phi  <- brms::get_dpar(prep, "phi")
  nc <- ncol(mu)
  nr <- nrow(mu)
  fup <- matrix(prep$data$vreal1, nrow=nr, ncol=nc, byrow=TRUE)
  finiteA <- matrix(prep$data$vint1, nrow=nr, ncol=nc, byrow=TRUE)
  finiteB <- matrix(prep$data$vint2, nrow=nr, ncol=nc, byrow=TRUE)
  ##trials <- matrix(draws$data$trials, nrow=nr, ncol=nc, byrow=TRUE)
  res <- matrix(NA, nrow=nr, ncol=nc)
  for(i in seq_len(nc)) {
    res[,i] <- inv_log(log_combo2_inv_link(mu[,i], phi[,i], fup[,i], finiteA[,i], finiteB[,i]))
  }
  res
}



## Note: Whenever a drug is absent, then we have to handle this in a
## special way. That is, the log(0) would be negative infinity which
## throws over Stan numerics. This is why we set in these cases the
## standardized drug predictor (log of scaled dose) to be 0 AND we do
## mark this with the finite columns that the drug is not present by
## setting it to 0 and 1 otherwise.
mutate_dose_cov <- function(data, dose_col, dref) {
    std_dose_col <- paste0("std_", dose_col)
    finite_dose_col <- paste0("finite_", dose_col)
    mutate(data,
           !!std_dose_col := if_else(!!sym(dose_col) == 0, 0, log(!!sym(dose_col) / dref)),
           !!finite_dose_col := 1 * (!!sym(dose_col) != 0)
           )
}

as_group_factor_nonhierarchical <- function(group_id) {
  levs <- c("trial")
  assert_that(all(group_id %in% levs))
  factor(group_id, levels=levs)
}

get_rel_time <- function(dt, ref_dt) {
    ## return time counted since the reference time-point in full days
    ceiling(as.numeric(as_datetime(dt) - as_datetime(ref_dt), "days"))
}

## filters NA's from a vector and defines it as accrual date + 3
## cycles
impute_default_cens_dt <- function(cens_dt, accrual_dt) {
    is_na <- is.na(cens_dt)
    cens_dt[is_na] <- accrual_dt[is_na] + 3 * cycle_length
    cens_dt
}

## given data with 1 row per patient in date-time format, the function
## converts it to a suitably long format for fitting purpose. That
## means to convert the date-time times into correct times measured in
## units of days.
prepare_data <- function(data_dt,  analysis_caltime_dt=today(),  summarise=TRUE) {

    assert_that(all(
      data_dt[["cens_caltime_dt"]] <= data_dt[["dropout_caltime_dt"]]
      )
    )
  
    ddw <- mutate(data_dt,
                  cens_time =get_rel_time(cens_caltime_dt, cycle1_caltime_dt),
                  event_time=get_rel_time(event_caltime_dt, cycle1_caltime_dt),
                  analysis_time=get_rel_time(analysis_caltime_dt, cycle1_caltime_dt),
                  cycle2_time=get_rel_time(cycle2_caltime_dt, cycle1_caltime_dt),
                  cycle3_time=get_rel_time(cycle3_caltime_dt, cycle1_caltime_dt),
                  id=seq_len(nrow(data_dt))) %>%
        select(id, group_id, drugA, drugB, cens_time, event_time, cens_time, analysis_time, cycle2_time, cycle3_time)
    ## build up data set using tmerge from survival, see
    ## https://cran.r-project.org/web/packages/survival/vignettes/timedep.pdf
    ## we take the cycle start as an event here (which ends the previous cycle)
    ##return(ddw)
    ddl <- tmerge(ddw, ddw, id=id, tstop=pmin(event_time, cens_time, analysis_time))
    if(any(is.finite(ddw$cycle2_time)))
        ddl <- tmerge(ddl, ddw, id=id, cycle_end=event(cycle2_time))
    if(any(is.finite(ddw$cycle3_time)))
        ddl <- tmerge(ddl, ddw, id=id, cycle_end=event(cycle3_time))
    ## enumerate the cycles
    
    if(length(ddl$tstart) > 0)
    {
      ddl <- tmerge(ddl, ddl, id=id, cycle=cumtdc(tstart))
    } else {
      ddl$cycle <- numeric(0)
    }
    
    ## finally add any event should there be one
    
    if(length(ddl$event_time) > 0)
    {
      ddl <- tmerge(ddl, ddw, id=id, dlt=event(event_time))
    } else {
      ddl$dlt <- numeric(0)
    }
    
    
    ddl <- ddl %>% mutate(
      tstop_orig = tstop,
      tstop = if_else(group_id=="trial", pmax(tstop, tstart + cycle_length_d), tstop),
      cycle_incomplete = if_else(
        group_id=="trial", 
        tstop != tstop_orig,
        FALSE
      ),
      drop = if_else(
        group_id=="trial", 
        cycle_incomplete & (dlt == 0), 
        FALSE
      )
   ) 
    
    ddl <- ddl %>%
      filter(!drop) %>%
      mutate(drop=NULL)
   
    if (summarise) {
      ddl <- ddl %>%
        group_by(group_id, drugA, drugB, tstart, tstop, cycle) %>%
        summarise(
          num_toxicities=sum(dlt), 
          num_patients=n(), 
          follow_up=sum(tstop - tstart), 
          num_patients_2=length(unique(id)), 
          .groups="drop"
        )
    } else { 
        ddl <- ddl %>%
            mutate(
              cycle_end=NULL, 
              follow_up = tstop - tstart, 
              num_patients=1, 
              num_toxicities=dlt
            )
    }
    ddl <- ddl %>%
        ungroup() %>%
        mutate(ocycle=factor(cycle, 1:3, ordered=TRUE)) %>% #, group_id=as_group_factor(group_id)) %>%
        mutate_dose_cov("drugA", drugA_dose_ref) %>%
        mutate_dose_cov("drugB", drugB_dose_ref)
    
    
    ddl %>% as.data.frame()
}



blrm_combo2_tte_fit_nonhierarchical  <- function(combo2_tte_model_nonhierarchical, data_dt, analysis_date=today(), seed=356456, ...) {
  if(nrow(data_dt) == 0) {
    return(
      update(
        combo2_tte_model_nonhierarchical, 
        sample_prior = "only", 
        control=list(adapt_delta=0.95), 
        cores=1, seed=seed, refresh=0, drop_unused_levels = FALSE, 
        ...
      )
    )
  } 
  update(
    combo2_tte_model_nonhierarchical, 
    sample_prior = "no",
    newdata=prepare_data(data_dt, analysis_date), 
    control=list(adapt_delta=0.95),
    cores=1, seed=seed, refresh=0, drop_unused_levels = FALSE, 
    ...
  )
  
}


## predicts the risk of an event up to a time-point
predict_risk <- function(model, newdata) {
    prep_newdata <- prepare_data(newdata, summarise=FALSE)## %>%
    ## obtain log(lambda)
    post_logH_pieces <- log(fitted(model, newdata=prep_newdata, allow_new_levels=TRUE, sample_new_levels="gaussian", summary=FALSE, scale="response"))
    post_logH <- sapply(seq_len(max(prep_newdata$id)), function(piece) {
        apply(post_logH_pieces[, prep_newdata$id == piece, drop=FALSE], 1, matrixStats::logSumExp)
    })
    ## P(T =< t) = 1-P(T > t) = inv_cloglog(log(H(t)))
    inv_cloglog(post_logH)
}


## predicts the risk of an event up to a time-point, conditional on
## having survived the previous cycles
predict_crisk <- function(model, newdata) {
    prep_newdata <- prepare_data(newdata, summarise=FALSE)
    ##print(prep_newdata)
    ##
    
    post_logH_pieces <- log(fitted(model, newdata=prep_newdata, allow_new_levels=TRUE, sample_new_levels="gaussian", summary=FALSE, scale="response"))
    ##print(post_logH_pieces)
    post_logH_delta <- sapply(seq_len(max(prep_newdata$id)), function(piece) {
        post_logH_case <- post_logH_pieces[, prep_newdata$id == piece, drop=FALSE]
        ## logH up to the end of a period
        post_logH_end <- apply(post_logH_case, 1, matrixStats::logSumExp)
        ## logH up to the start of a period
        num_periods <- ncol(post_logH_case)
        ##print(num_periods)
        if(num_periods == 1) {
            post_logH_start  <- -Inf
        } else {
            post_logH_start <- apply(cbind(post_logH_case[, -ncol(post_logH_case), drop=FALSE]), 1, matrixStats::logSumExp)
        }
        ##print(post_logH_start)
        ## (Fstart - Fend) / Fstart
        ## = 1- exp(- Delta H)
        ## logH_delta <- log( exp(post_logH_end) - exp(post_logH_start) )
        logH_delta <- post_logH_end + log1p( - exp(post_logH_start-post_logH_end) )
    })
    inv_cloglog(post_logH_delta)
}

summarise_risk <- function(post_risk) {
    df_1 <- t(apply(post_risk, 2, function(r) {
        quants  <- unname(quantile(r, c(0.25, 0.5, 0.75)))
        s <- c(mean=mean(r), sd=sd(r), q25=quants[1], q50=quants[2], q75=quants[3])
    })) %>% as.data.frame()
    
    interval_probs <- posterior::summarize_draws(
      post_risk,
      list(
        function(x){
          prop.table(table(cut(x, breaks = c(0, 0.16, 0.33, 1))))
        }
      )
    ) %>% rename(
      "prob_underdose" = `(0,0.16]`,
      "prob_target" = `(0.16,0.33]`,
      "prob_overdose" = `(0.33,1]`
    ) %>% select(-variable)
    
    bind_cols(df_1, interval_probs) %>% mutate(ewoc_ok = q75 < 0.33)
}



## piece-wise simulation of piece-wise constant hazard
rpexp <- function(lambda, epochs) {
  num_periods <- length(epochs) + 1
  
  assert_that(!anyDuplicated(epochs))
  assert_that(all(epochs == sort(epochs)))
  
  assert_that(length(lambda) == num_periods)
  
  ## sample for each period an event time
  event_times_nonzero <- rexp(sum(lambda != 0), lambda[lambda != 0])
  
  event_times <- numeric(num_periods)
  event_times[lambda != 0] <- event_times_nonzero
  event_times[lambda == 0] <- Inf
  ## check if the sampled time per period falls inside the epoch
  epochs_duration <- diff(c(0, epochs, Inf))
  
  ## ... and find first epoch with an event
  event_happened_in_epoch <- event_times < epochs_duration
  
  ## sampled time is the start of the epoch plus the within-epoch time
  if (any(event_happened_in_epoch)) {
    first_event_epoch <- which(event_happened_in_epoch)[1]
    sampled_time <- event_times[first_event_epoch] + c(0,epochs)[first_event_epoch]
  } else {
    sampled_time <- Inf
  }
  
  
  if (any(is.na(sampled_time)) || length(sampled_time) != 1) {
    browser()
    stop("rpexp error")
  }
  
  return(sampled_time)
}

plot_crisk <- function(model, schedules) {
    post_crisk <- predict_crisk(model, schedules)
    ##
    bind_cols(schedules, summarise_risk(post_crisk)) %>%
        mutate(ewoc_ok=factor(q75 < 0.33, c(TRUE, FALSE), c("OK", "Not OK")), drugA=factor(drugA)) %>%
        ggplot(aes(factor(drugB), q50, colour=ewoc_ok, shape=drugA, linetype=drugA)) +
        geom_pointrange(aes(ymin=q25, ymax=q75), position=position_dodge(width=0.3)) +
        facet_wrap(~num_cycles, labeller=label_both) +
        scale_y_continuous(breaks=seq(0,1,by=0.1)) +
        coord_cartesian(ylim=c(0, 0.5)) +
        hline_at(0.33, linetype=I(2)) +
        ggtitle("Conditional risk for one DLT per cycle", "Risk is conditional on survival up to cycle start.\nShown is the median (dot) and central 50% CrI (line).") +
        xlab("Dose [mg]") +
        ylab("P(DLT in cycle| survival to cycle start)")
}

plot_risk_ob2 <- function(blrm_1cycle, blrm_3cycle, schedules) {
  post_risk_1cycle <- summary(blrm_1cycle, "dose_prediction", newdata = schedules, prob=c(0.25, 0.5, 0.75))
  post_risk_3cycle <- summary(blrm_3cycle, "dose_prediction", newdata = schedules, prob=c(0.25, 0.5, 0.75))
  ##
  post_risk <- bind_rows(
    mutate(post_risk_1cycle, num_cycles = 1),
    mutate(post_risk_1cycle[1, ], num_cycles = 2, ewoc_ok = NA, `50%` = NA, `25%` = NA, `75%` = NA),
    mutate(post_risk_3cycle, num_cycles = 3)
  )
  
  post_risk$cycle_label <- paste0("Cycle ", post_risk$num_cycles)
  
  
  post_risk %>%
    mutate(ewoc_ok=factor(ewoc_ok, c(TRUE, FALSE), c("OK", "Not OK")), drugA=factor(drugA)) %>%
    ggplot(aes(factor(drugB), `50%`, colour=ewoc_ok, shape=drugA, linetype=drugA)) +
    geom_pointrange(aes(ymin=`25%`, ymax=`75%`), position=position_dodge(width=0.3)) +
    facet_wrap(~cycle_label) +
    #facet_wrap(~num_cycles, labeller=label_both) +
    scale_y_continuous(breaks=seq(0,1,by=0.1)) +
    coord_cartesian(ylim=c(0, 0.5)) +
    hline_at(0.33, linetype=I(2)) +
    ggtitle("Overall risk for a DLT, BLRM", "Shown is the median (dot) and central 50% CrI (line).") +
    xlab("Dose [mg]") +
    ylab("P(DLT)")
}

plot_risk <- function(model, schedules) {
  post_risk <- predict_risk(model, schedules)
  
  risk_summary <- bind_cols(schedules, summarise_risk(post_risk))
  
  risk_summary$cycle_label <- paste0("Cycle ", risk_summary$num_cycles)
  ##
  risk_summary %>%
    mutate(ewoc_ok=factor(q75 < 0.33, c(TRUE, FALSE), c("OK", "Not OK")), drugA=factor(drugA)) %>%
    ggplot(aes(factor(drugB), q50, colour=ewoc_ok, shape=drugA, linetype=drugA)) +
    geom_pointrange(aes(ymin=q25, ymax=q75), position=position_dodge(width=0.3)) +
    facet_wrap(~cycle_label) +
    scale_y_continuous(breaks=seq(0,1,by=0.1)) +
    coord_cartesian(ylim=c(0, 0.5)) +
    hline_at(0.33, linetype=I(2)) +
    ggtitle("Overall risk for a DLT, TTE", "Shown is the median (dot) and central 50% CrI (line).") +
    xlab("Dose [mg]") +
    ylab("P(DLT)")
}

plot_crisk_time <- function(model, schedules) {
    schedule_data_time <- schedules %>%
        expand_grid(day=round(seq(1, cycle_length_d, length.out=3*cycle_length_w))) %>%
        mutate(cens_caltime_dt = cycle1_caltime_dt + (num_cycles-1) * cycle_length + days(day))

    post_crisk_time <- predict_crisk(model, schedule_data_time)

    bind_cols(schedule_data_time, summarise_risk(post_crisk_time)) %>%
        mutate(ewoc_ok=q75 < 0.33) %>%
        ggplot(aes(day/7, q50)) +
        geom_ribbon(aes(ymin=q25, ymax=q75), fill="grey70") +
        geom_line() +
        facet_grid(vars(drugB), vars(drugA, num_cycles), labeller=label_both) +
        scale_y_continuous(breaks=seq(0,1,by=0.1)) +
        scale_x_continuous(breaks=seq(0,cycle_length_w,by=1)) +
        coord_cartesian(ylim=c(0, 0.5)) +
        hline_at(0.33, linetype=I(2)) +
        ggtitle("P(DLT during t| survival to cycle start)", "Risk is conditional on survival up to cycle start.\nShown is the median (line) and central 50% CrI (grey area).") +
        xlab("Time since cycle start [week]") +
        ylab("P(DLT during t| survival to cycle start)")
}

plot_risk_time <- function(model, schedules) {
    max_cycles <- max(schedules$num_cycles)
    schedule_data_time <- schedules %>%
        group_by(group_id, sid, drugA, drugB) %>%
        filter(num_cycles==max(num_cycles)) %>%
        ungroup() %>%
        expand_grid(day=round(seq(1, max_cycles*cycle_length_d, length.out=max_cycles*cycle_length_w))) %>%
        mutate(cens_caltime_dt = cycle1_caltime_dt + days(day))

    post_risk_time <- predict_risk(model, schedule_data_time)

    bind_cols(schedule_data_time, summarise_risk(post_risk_time)) %>%
        ggplot(aes(day/7, q50)) +
        geom_ribbon(aes(ymin=q25, ymax=q75), fill="grey70") +
        geom_line() +
        facet_grid(vars(drugB), vars(drugA, num_cycles), labeller=label_both) +
        scale_y_continuous(breaks=seq(0,1,by=0.1)) +
        scale_x_continuous(breaks=unique(c(round(seq(0,max_cycles*cycle_length_w,length.out=11)), 0, cycle_length_w, 2*cycle_length_w, 3*cycle_length_w))) +
        coord_cartesian(ylim=c(0, 1.0)) +
        hline_at(0.33, linetype=I(2)) +
        vline_at(c(cycle_length_w, 2*cycle_length_w, 3*cycle_length_w), linetype=I(3)) +
        ggtitle("P(DLT during t)", "Risk for a DLT up to time-point t.\nShown is the median (line) and central 50% CrI (grey area).") +
        xlab("Time since treatment start [week]") +
        ylab("P(DLT during t)")
}

maximal_dose <- function(model, predict_fun, newdata, dose_col) {
  samp <- predict_fun(model,  newdata=newdata)
  doses <- pull(newdata, {{ dose_col }})
  overdose_metric <- c(-0.5, summarize_draws(samp, ~quantile2(.x, 0.75) - 0.33)$q75, 1)
  floor(uniroot(approxfun(c(-1, doses, Inf), overdose_metric), range(doses))$root)
}

plot_max_dose_ob2 <- function(model, dose_col, cycle_num = 1)
{
  dose_set <- summary(model, "dose_info")[[dose_col]]
  
  dose_min <- min(dose_set)
  dose_max <- max(dose_set)
  
  fake_dens <- tibble(stratum_id="all", group_id="trial", drugA = 1, drugB = seq(from = dose_min, to = dose_max, length.out = 101))
  # fake_dens$drug_A = fake_dens[[dose_col]]
  
  prior_blrm_samp <- posterior_linpred(model, transform = TRUE, newdata=fake_dens)
  max_dose_ok <- maximal_dose(model, predict_fun = \(x, newdata){posterior_linpred(x, transform = TRUE, newdata=newdata)}, fake_dens, dose_col)
  
  dref <- filter(summary(model, "drug_info"), drug_name == dose_col)[["dose_ref"]]
  
  max_dose_plot <- ppd_ribbon(prior_blrm_samp, fake_dens[[dose_col]]) + scale_x_log10(breaks=dose_set) +
    scale_y_continuous(breaks=c(0,0.16, 0.33, seq(0.5,1,by=0.25))) +
    hline_at(c(0.16,0.33), linetype=2) + vline_at(max_dose_ok, linetype=1) +
    annotate("label", max_dose_ok, 0.9, label=paste0("Maximal dose:\n", max_dose_ok, "mg")) +
    annotate("label", 15, 0.9, label=paste0("B", cycle_num)) +
    xlab("Dose [mg]") + ylab("P(DLT)") +
    coord_cartesian(ylim=c(0,1)) +
    ggtitle("Logistic model prior toxicity probability", paste0("Reference dose ", dref, ", cycle ", cycle_num))
  
  return(max_dose_plot)
}

plot_max_dose_tte <- function(model, dose_col, tte_schedules, drug_info)
{
  dose_set <- tte_schedules[[dose_col]]
  
  dose_min <- min(dose_set[dose_set > 0])
  dose_max <- max(dose_set)
  
  tite_fake_dens <- NULL
  tite_fake_dens_template <- filter(tte_schedules, drugB == dose_min, drugA == 1)
  
  i <- 1
  for (drugB_dose in seq(from = dose_min, to = dose_max, length.out = 101)) {
    tite_fake_dens <- tite_fake_dens %>% bind_rows(
      tite_fake_dens_template %>% mutate(
        drugB = drugB_dose,
        dose = drugB_dose,
        sid = i
      )
    )
    i <- i + 1
  }
  
  tite_fake_dens_prepared <- prepare_data(tite_fake_dens)
  prior_tite_samp <-  predict_risk(model, newdata=tite_fake_dens)
  
  tite_max_dose <- tite_fake_dens %>% group_by(cycle_index) %>%
    summarize(max_ok = maximal_dose(model, predict_risk, pick(everything()), drugB)) %>%
    mutate(group=cycle_index)
  
  tite_plot <- list()
  for (cycle_index in 1:3) {
    cur_cycle_index <- tite_fake_dens$cycle_index == cycle_index
    
    tite_plot[[cycle_index]] <- ppd_ribbon(prior_tite_samp[,cur_cycle_index], tite_fake_dens$drugB[cur_cycle_index]) + 
      scale_x_log10(breaks=dose_set) +
      scale_y_continuous(breaks=c(0,0.16, 0.33, seq(0.5,1,by=0.25))) +
      hline_at(c(0.16,0.33), linetype=2) +
      geom_vline(aes(xintercept=max_ok), tite_max_dose[cycle_index,]) +
      geom_label(aes(x=max_ok, label=paste0("Maximal dose:\n", max_ok, "mg")), tite_max_dose[cycle_index,], y=0.9, inherit.aes=FALSE) +
      annotate("label", 15, 0.9, label=ifelse(cycle_index == 1, "TCO\nTCU", "TCU")) +
      xlab("Dose [mg]") + 
      ylab("P(DLT)") +
      coord_cartesian(ylim=c(0,1), xlim=c(dose_min, dose_max)) +
      ggtitle("TTE model prior toxicity probability", paste0("Reference dose ", filter(drug_info, drug_name == dose_col)$dose_ref, ", up to end of cycle ", cycle_index))
  }
  
  return(tite_plot)
}

plot_max_cdose_tte <- function(model, dose_col, tte_schedules, drug_info)
{
  dose_set <- tte_schedules[[dose_col]]
  
  dose_min <- min(dose_set[dose_set > 0])
  dose_max <- max(dose_set)
  
  tite_fake_dens <- NULL
  tite_fake_dens_template <- filter(tte_schedules, drugB == dose_min, drugA == 1)
  
  i <- 1
  for (drugB_dose in seq(from = dose_min, to = dose_max, length.out = 101)) {
    tite_fake_dens <- tite_fake_dens %>% bind_rows(
      tite_fake_dens_template %>% mutate(
        drugB = drugB_dose,
        dose = drugB_dose,
        sid = i
      )
    )
    i <- i + 1
  }
  
  
  prior_tite_csamp <-  predict_crisk(model, newdata=tite_fake_dens)
  
  tite_max_cdose <- tite_fake_dens %>% group_by(cycle_index) %>%
    summarize(max_ok = maximal_dose(model, predict_crisk, pick(everything()), drugB)) %>%
    mutate(group=cycle_index)
  
  tite_cplot <- list()
  
  for (cycle_index in 1:3) {
    cur_cycle_index <- tite_fake_dens$cycle_index == cycle_index
    
    tite_cplot[[cycle_index]] <- ppd_ribbon(prior_tite_csamp[,cur_cycle_index], tite_fake_dens$drugB[cur_cycle_index]) + 
      scale_x_log10(breaks=dose_set) +
      scale_y_continuous(breaks=c(0,0.16, 0.33, seq(0.5,1,by=0.25))) +
      hline_at(c(0.16,0.33), linetype=2) +
      geom_vline(aes(xintercept=max_ok), tite_max_cdose[cycle_index,]) +
      geom_label(aes(x=max_ok, label=paste0("Maximal dose:\n", max_ok, "mg")), tite_max_cdose[cycle_index,], y=0.9, inherit.aes=FALSE) +
      annotate("label", 15, 0.9, label="TCO") +
      xlab("Dose [mg]") + 
      ylab("P(DLT)") +
      coord_cartesian(ylim=c(0,1), xlim=c(dose_min, dose_max)) +
      ggtitle("TTE model prior conditional toxicity probability", paste0("Reference dose ", filter(drug_info, drug_name == dose_col)$dose_ref, ", cycle ", cycle_index))
  }
  
  return(tite_cplot)
}

