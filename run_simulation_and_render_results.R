source("renv/activate.R")
source("install_packages.R")


# Comment to below to use the provided results instead of re-running the 
# operating characteristics simulation (approx. 2 hours on 300 cores)
source("combo2-tte-oc-main.R") 

# Render quarto results report
source("render_results.R")

# Render the example patient trajectory figure
source("render_fake_patient_trajectories.R")