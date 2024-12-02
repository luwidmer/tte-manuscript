# README

This repository hosts the TeX and R code for the manuscript

*TITE-CLRM: Towards efficient time-to-event dose-escalation guidance of multi-cycle cancer therapies* <https://arxiv.org/abs/TODO>

Reproduce by firing up `run_model_and_compile_manuscript.R` in R in this folder (or in RStudio after loading the project).

## Required software

-   R 4.3.1
-   CmdStan version 2.32.2 (cmdstanr 0.8.1)
-   renv 1.0.7
-   R packages will be restored by {renv}

## How to run the simulation study and generate the results plots and tables

Run `run_simulation_and_render_results.R` with R 4.3.1. This will generate a results folder with the time and date of the simulation in the `results/` folder. Note that this will run for approximately one hour on 200 CPU cores (a high-performance compute cluster is highly recommended). Changing this can be done in `combo2-tte-oc-main.R`. After simulations complete, a report is rendered and patient trajectory plots are generated in the `fake` folder. Note that if running `combo2-tte-oc-main.R` is skipped, the results are re-generated from the provided ones.