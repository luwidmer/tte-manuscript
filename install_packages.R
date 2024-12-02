renv::restore()

required_packages <- c(
  "here",
  "clustermq",
  "tidyverse",
  "assertthat",
  "lubridate",
  "bayesplot",
  "egg",
  "qs",
  "RColorBrewer",
  "OncoBayes2",
  "posterior",
  "brms",
  "survival",
  "knitr",
  "tidyr",
  "quarto",
  "remotes"
)

if (!all(required_packages %in% installed.packages()))
{
  install.packages(required_packages[!(required_packages %in% installed.packages())])
}

remotes::install_github("stan-dev/cmdstanr@v0.8.1", dependencies="never")
if (cmdstanr::cmdstan_version() != "2.32.2") {
  cmdstanr::install_cmdstan(version = "2.32.2")
}

cmdstan_path <- cmdstanr::cmdstan_path()
if (basename(cmdstan_path) != "cmdstan-2.32.2") {
  cmdstanr::set_cmdstan_path(path = file.path(dirname(cmdstan_path), "cmdstan-2.32.2"))
  cmdstanr::check_cmdstan_toolchain()
}
