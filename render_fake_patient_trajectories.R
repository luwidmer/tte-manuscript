library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(purrr)
library(survival)
library(assertthat)
library(bayesplot)
theme_set(theme_bw(12))


## cycle length in weeks
cycle_length <- weeks(6)
cycle_length_d <- as.numeric(cycle_length, "days")
cycle_length_w <- as.numeric(cycle_length, "weeks")
cycle_length_h <- as.numeric(cycle_length, "hours")

## pseudo date used to mark distant / unkown future time-point
inf_date  <- as_date(Inf)
inf_datetime  <- as_datetime(Inf)


for (i in 1:10) {
  fake_data_file <- paste0("fake/fake_patients-", i, ".RDS")
  export_basename <- str_remove(fake_data_file, ".RDS$")
  
  fake <- readRDS(fake_data_file) |> 
    mutate(label=factor(paste(patient_id, dose_shifted, sep="/")))
  
  get_rel_time <- function(dt, ref_dt) {
    ## return time counted since the reference time-point in full days
    ceiling(as.numeric(as_datetime(dt) - as_datetime(ref_dt), "days"))
  }
  
  ## in each cycle these things can happen:
  ## - event
  ## - no event
  ##
  ## requirements:
  ## - cycle needs to be eligible (100% of dose taken)
  ## - no dropout must have happened
  ## - analysis cut-off must be late enough
  ## - patient still being followed up
  ##
  ## each patient is either at risk or not at risk in every time window.
  ## whenever an event happens in a time-window, then all following
  ## windows are not at risk.
  
  
  
  ## takes in subject history data with one data row per patient and
  ## outputs counting format style data. That is we split for each
  ## patient row his data over cycles marking start and stop. For each
  ## patient the output data contains always the maximally observable
  ## number of cycles, but marks which ones were actually observed
  ## (status column) and during which time-periods the patient was
  ## at-risk. What is considered as an event is configurable via the
  ## incidence argument which defines columns to search for events.
  counting_format <- function(subject_history, id="patient_id", analysis_caltime_dt=as_date(Inf), incidences) {
    num_cycles <- max(as.numeric(str_remove(str_extract(names(fake), "[0-9]+cycle"), "cycle")), na.rm=TRUE)
    counting_long <- expand_grid(id=pull(subject_history, {{id}}), cycle=1:num_cycles)
    names(counting_long)[1] <- id
    
    subject_history$analysis_caltime_dt <- analysis_caltime_dt
    all_incidences <- c(incidences, "analysis_caltime_dt")
    
    ## extract from wide data cycle specific dates and times
    cycle_timing_long <- select(subject_history, all_of(id), starts_with("cycle") & ends_with("caltime_dt"), end_of_last_cycle_caltime_dt) |>
      pivot_longer(-all_of(id), names_to="cycle", values_to="start_caltime_dt") %>%
      mutate(cycle=as.integer(str_extract(cycle, "[0-9]"))) |>
      group_by(!!sym(id)) |>
      mutate(end_caltime_dt=c(start_caltime_dt[-1], Inf)) |>
      filter(!is.na(cycle))
    
    ## take over cycle specific dates like eligible
    cycle_data_long_date <- select(subject_history, all_of(id), !starts_with("cycle") & matches("[0-9]+cycle") & where(is.date)) |>
      pivot_longer(-all_of(id), names_to="cycle_var", values_to="value") |>
      mutate(cycle=as.integer(str_remove(str_extract(cycle_var, "[0-9]cycle"), "cycle")),
             cycle_var=str_replace(str_remove(cycle_var, "[0-9]+cycle"), "__", "_")) |>
      pivot_wider(id_cols = all_of(c(id, "cycle")), names_from="cycle_var")
    
    ## find out what was actually observed
    obs_incidence <- select(subject_history, all_of(id), all_of(all_incidences)) |>
      pivot_longer(-all_of(id), names_to="obs_incidence", values_to="obs_caltime_dt") |>
      group_by(!!sym(id)) |>
      slice_min(obs_caltime_dt) |>
      mutate(obs_incidence=str_remove(obs_incidence, "_caltime_dt"))
    
    baseline_data <- select(subject_history, all_of(id), !(contains("cycle") | any_of(all_incidences)))
    
    ## now merge this data to full data set
    counting <- reduce(list(counting_long, cycle_timing_long, cycle_data_long_date), left_join, by=c(id, "cycle")) |>
      left_join(obs_incidence, by=id)
    
    counting <- counting |>
      mutate(
        at_risk = (start_caltime_dt <= obs_caltime_dt), 
        eligible = end_caltime_dt <= obs_caltime_dt | obs_incidence == "event", 
        at_risk_eligible = (at_risk & eligible),
        status = case_when(
          !at_risk ~ "unobserved",
          obs_caltime_dt == end_caltime_dt & obs_incidence == "end_of_last_cycle" ~ "censored",
          ## we add 1s to the event to make sure it is in the cycle
          !( (obs_caltime_dt+seconds(1)) %within% interval(start_caltime_dt, end_caltime_dt)) ~ "at-risk",
          obs_incidence == "event"  & (obs_caltime_dt <= analysis_caltime_dt) ~ "event",
          obs_incidence == "event"  & (obs_caltime_dt > analysis_caltime_dt) ~ "censored",
          !eligible ~ "non-eligible",
          obs_incidence == "dropout" ~ "censored",
          obs_incidence == "analysis" ~ "censored",
          TRUE ~ "unknown" ## should never happen!!!
        )
      )
    
    counting |> 
      group_by(!!sym(id)) |>
      mutate(
        status_next = c(status[-1], "unknown"), 
        status = ifelse(status_next == "non-eligible", "censored", status)
      ) |>
      select(-status_next) |>
      ungroup() |>
      left_join(baseline_data, by=id)
  }
  
  study_history_counting <- counting_format(
    mutate(fake, origin_caltime_dt = min(cycle1_caltime_dt)), 
    incidences = c("event_caltime_dt", "dropout_caltime_dt", "end_of_last_cycle_caltime_dt")
  ) |>
    filter(at_risk) |> 
    select(-at_risk)
  
  study_history_counting |> as.data.frame()
  
  analysis_counting_format <- function(subject_history, analysis_caltime_dt, origin="cycle1_caltime_dt") {
    counting_format(mutate(subject_history, origin_caltime_dt=!!sym(origin)),
                    analysis_caltime_dt=analysis_caltime_dt,
                    incidences=c("event_caltime_dt", "dropout_caltime_dt", "end_of_last_cycle_caltime_dt")) |>
      filter(at_risk) |> select(-at_risk) |>
      mutate(across(ends_with("caltime_dt"), ~get_rel_time(.x, origin_caltime_dt)/7)) |>
      rename_with(~str_replace(.x, "caltime_dt$", "w"), ends_with("caltime_dt")) |> 
      select(-origin_w)
  }
  
  study_history_counting_w <- analysis_counting_format(mutate(fake, fpfv_caltime_dt=min(cycle1_caltime_dt)), as_date(Inf), origin="fpfv_caltime_dt")
  study_history_w_terminal <- filter(study_history_counting_w, status != "at-risk") |>
    mutate(obs_incidence=factor(obs_incidence, c("event", "dropout", "end_of_last_cycle"), c("DLT", "dropout", "last cycle\nend")))
  
  fpfv_caltime_dt <- min(fake$cycle1_caltime_dt)
  dem1_caltime_dt <- max(pmin(fake$cycle2_caltime_dt[1:3], fake$dropout_caltime_dt[1:3])) 
  dem2_caltime_dt <- max(pmin(fake$cycle2_caltime_dt[4:6], fake$dropout_caltime_dt[4:6])) 
  fake
  fpfv_caltime_dt + weeks(12)
  
  dem_dates <- data.frame(
    dem_dt = c(dem1_caltime_dt, dem2_caltime_dt)
  ) |> 
    mutate(Analysis=1:n())
  
  
  study_time <- ggplot(study_history_counting_w , aes(x=label)) +
    geom_linerange(linewidth=1, aes(ymin=start_w, ymax=end_w, colour=factor(cycle))) + 
    geom_hline(data=dem_dates, aes(yintercept=get_rel_time(dem_dt, fpfv_caltime_dt)/7), linetype=I(2)) +
    geom_point(data=study_history_w_terminal, aes(y=obs_w, shape=obs_incidence), size=3) +
    geom_label(
      data=dem_dates, 
      aes(y=get_rel_time(dem_dt, fpfv_caltime_dt)/7, label=paste("Analysis", Analysis)), 
      x=c(8.5, 8.5), 
      size=3
    ) +
    scale_y_continuous(breaks=seq(0, 200, by=cycle_length_w)) +
    coord_flip() +
    scale_colour_brewer("Cycle", type="qual", palette="Dark2") +
    scale_shape_manual("Incidence", values=c(0, 1, 2), drop=FALSE) +
    ylab("Study time since first patient cycle 1 start [weeks]") +
    xlab("Patient id / Dose [mg]") +
    theme(legend.position="right") +
    ggtitle("Study history")
  
  
  study_time
  
  long_fake_milstone <- dem_dates |>
    rowwise() |>
    mutate(
      analysis_sub = list(analysis_counting_format(fake, dem_dt) |> filter(at_risk_eligible)), 
      analysis_terminal_sub=list(filter(analysis_sub, status != "at-risk"))
    )
  
  
  long_fake_milstone_history <- long_fake_milstone %>%
    select(-analysis_terminal_sub) %>%
    unnest(analysis_sub)
  
  long_fake_milstone_terminal <- long_fake_milstone %>%
    select(-analysis_sub) %>%
    unnest(analysis_terminal_sub) %>%
    mutate(status=factor(status, c("event", "censored")))
  
  
  long_fake_milstone_terminal
  
  subset(long_fake_milstone_history, label=="4/150")
  
  analysis_counting_format(fake, ymd("2022-05-24"))
  
  
  patient_time_snap <- long_fake_milstone_history %>%
    ggplot(aes(x=label)) +
    facet_grid(. ~ Analysis, labeller=label_both) +
    geom_linerange(aes(ymin=start_w, ymax=end_w, colour=factor(cycle))) +
    geom_point(data=long_fake_milstone_terminal, aes(y=end_w, shape=status), size=3) +
    scale_y_continuous(breaks=seq(0, 200, by=cycle_length_w)) +
    scale_x_discrete(drop=FALSE) +
    ##geom_hline(data=dem_dates, aes(yintercept=get_rel_time(dem_dt, fpfv_caltime_dt)/7), linetype=I(2)) +
    coord_flip() +
    ##scale_shape_discrete("Class") +
    scale_shape_manual("Class", values=c(15, 16, 17), drop=FALSE) +
    scale_colour_brewer("Cycle", type="qual", palette="Dark2") +
    ylab("Patient time [weeks]") +
    ##xlab("Patient id / Dose [mg]") +
    xlab(NULL) +
    ggtitle("Patient history per analysis") +
    guides(colour="none") +
    theme(
      legend.position = "inside", 
      legend.position.inside = c(0.18, 0.75), 
      legend.box.background = element_rect(colour = "black")
    )
  
  patient_time_snap
  
  plot_combined <- bayesplot_grid(study_time, patient_time_snap, grid_args=list(nrow=1, widths=c(1.4,1)))
  
  plot_combined
  
  ggsave(paste0(export_basename, "-combined_history.pdf"), plot_combined, width=1.63*3*2, height=3.1)
  
}
