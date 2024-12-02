
plot_dropout <- function(lambda_scenario) {
  plot_scenario_data <- NULL
  
  for (scenario_name in names(lambda_scenario)) {
    plot_scenario_data <- bind_rows(
      plot_scenario_data, 
      mutate(lambda_scenario[[scenario_name]], scenario_name = scenario_name, ewoc_ok = prob_true < 0.33)
    )
  }
  plot_scenario_data$scenario_name <- factor(plot_scenario_data$scenario_name, levels = names(lambda_scenario))
  plot_scenario_data$cycle_label <- paste0("Cycle ", plot_scenario_data$cycle_index)
  plot_scenario_data <- plot_scenario_data %>% mutate(
    prob_dropout = inv_cloglog(log(dropout_rate) + log(cycle_length_d))
  )
  doses <- unique(lambda_scenario[[scenario_name]]$dose)
  
  ggplot(
    plot_scenario_data, 
    aes(drugB, prob_dropout, linetype=factor(dropout_total) , color=dropout_type )
  ) +
    geom_line() +
    facet_wrap(~cycle_label) +
    scale_x_continuous(trans='log10', breaks=doses) +
    ggtitle("Dropout risk for a DLT, per-cycle / conditional", "True probability") +
    xlab("Dose [mg]") +
    ylab("Dropout probability") + 
    guides(linetype=guide_legend(title="Total Dropout"), color=guide_legend(title="Dropout Scenario")) +
    scale_colour_brewer(palette = "Dark2")
}

plot_scenarios <- function(lambda_scenario) {
  plot_scenario_data <- NULL
  
  for (scenario_name in names(lambda_scenario)) {
    
    plot_scenario_data <- bind_rows(
      plot_scenario_data, 
      mutate(lambda_scenario[[scenario_name]], scenario_name = scenario_name, ewoc_ok = prob_true < 0.33)
    )
  }
  plot_scenario_data$scenario_name <- factor(plot_scenario_data$scenario_name, levels = names(lambda_scenario))
  plot_scenario_data$cycle_label <- paste0("Cycle ", plot_scenario_data$cycle_index)
  
  doses <- unique(lambda_scenario[[scenario_name]]$dose)
  
  ggplot(plot_scenario_data, aes(drugB, prob_true, linetype=scenario_name, color=scenario_name)) +
    geom_line() +
    facet_wrap(~cycle_label) +
    scale_y_continuous(breaks=c(0,0.16, 0.33, seq(0.5,1,by=0.25))) +
    scale_x_continuous(trans='log10', breaks=doses) +
    coord_cartesian(ylim=c(0, 1)) +
    hline_at(c(0.16, 0.33), linetype=I(2)) +
    ggtitle("Scenario risk for a DLT, per-cycle / conditional", "True probability") +
    xlab("Dose [mg]") +
    ylab("P(DLT | no DLT so far)") + 
    guides(linetype=guide_legend(title="Scenario"), color=guide_legend(title="Scenario")) +
    scale_colour_brewer(palette = "Dark2")
}

plot_scenarios_cumulative <- function(lambda_scenario) {
  plot_scenario_data <- NULL
  for (scenario_name in names(lambda_scenario)) {
    plot_scenario_data <- bind_rows(
      plot_scenario_data, 
      mutate(lambda_scenario[[scenario_name]], scenario_name = scenario_name, ewoc_ok = prob_true < 0.33)
    )
  }
  plot_scenario_data$scenario_name <- factor(plot_scenario_data$scenario_name, levels = names(lambda_scenario))
  
  plot_scenario_data_cumulative <- filter(plot_scenario_data, cycle_index == 1 ) 
  plot_scenario_data_cumulative <- plot_scenario_data_cumulative %>% bind_rows(
    filter(plot_scenario_data, cycle_index == 2 ) %>%
      inner_join(
        filter(plot_scenario_data_cumulative, cycle_index == 1) %>% 
          select("drugB", "drugA", "scenario_name", "prob_true") %>%
          rename(prob_true_1 = prob_true), 
        by = c("drugB", "drugA", "scenario_name")
      ) %>%
      mutate(prob_true = prob_true_1 + (1-prob_true_1)*prob_true)
  )
  
  plot_scenario_data_cumulative <- plot_scenario_data_cumulative %>% bind_rows(
    filter(plot_scenario_data, cycle_index == 3 ) %>%
      inner_join(
        filter(plot_scenario_data_cumulative, cycle_index == 2) %>% 
          select("drugB", "drugA", "scenario_name", "prob_true") %>%
          rename(prob_true_2 = prob_true), 
        by = c("drugB", "drugA", "scenario_name")
      ) %>%
      mutate(prob_true = prob_true_2 + (1-prob_true_2)*prob_true)
  )
  
  plot_scenario_data_cumulative$cycle_label <- paste0("Cycle ", plot_scenario_data_cumulative$cycle_index)
  
  doses <- unique(lambda_scenario[[scenario_name]]$dose)
  
  ggplot(plot_scenario_data_cumulative, aes(drugB, prob_true, linetype=scenario_name, color=scenario_name)) +
    geom_line() +
    facet_wrap(~cycle_label) +
    scale_y_continuous(breaks=c(0,0.16, 0.33, seq(0.5,1,by=0.25))) +
    scale_x_continuous(trans='log10', breaks=doses) +
    coord_cartesian(ylim=c(0, 1)) +
    hline_at(c(0.16, 0.33), linetype=I(2)) +
    ggtitle("Scenario risk for a DLT, cumulative", "True probability") +
    xlab("Dose [mg]") +
    ylab("P(DLT)") + 
    guides(linetype=guide_legend(title="Scenario"), color=guide_legend(title="Scenario")) +
    scale_colour_brewer(palette = "Dark2")
}