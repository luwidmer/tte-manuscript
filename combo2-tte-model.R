# automatic caching of models when using "file" argument of brm

################################
## STAN MODEL DEFINTION
#########################################################

## define a custom family to get control over things
## additional input terms to formula are
## 1. number of trials (integer)
## 2. total follow-up time (not as log!)
## 3. if drug_A is present indicator 1 = present, 0 = absent
## 4. if drug_B is present indicator 1 = present, 0 = absent
blrm_combo2_tte <- custom_family(
  "blrm_combo2_tte", dpars = c("mu", "phi"),
  type = "int", vars = c("vreal1", "vint1", "vint2"),
  loop=FALSE
)


## Stan functions to make custom family come to live
blrm_combo2_tte_stan_code <- "
  // given log-rate of single agent 1 & 2, follow up time and event status
  // log(lambda) = log(exp(mu) + exp(phi)) + eta + offset(log(fup))
  // no interaction here, so eta=0 constant
  real log_combo2_inv_link(real mu, real phi, real fup, int finiteA, int finiteB) {
     // log(lambda)
     real log_lambda = log(fup);
     if(finiteA == 1 && finiteB == 1) {
        log_lambda += log_sum_exp(mu, phi);
     } else if(finiteA == 1 && finiteB != 1) {
        log_lambda += mu;
     } else if(finiteA != 1 && finiteB == 1) {
        log_lambda += phi;
     } else if(finiteA != 1 && finiteB != 1) {
        log_lambda = negative_infinity();
     }
     return log_lambda;
  }
/*
  real binomial_cloglog_lpmf(int y, real theta) {
     // cloglog(p) = theta = log(-log(1-p)) => need log(p) & log(1-p)
     real log_np = -exp(theta);
     real log_p = log1m_exp(log_np);
     int nr = T - y;
     // drop normalization here, we would need log (y over T)
     return y * log_p + nr * log_np;
  }
*/
  real blrm_combo2_tte_lpmf(int[] y, vector mu, vector phi, real[] fup, int[] finiteA, int[] finiteB) {
    int num_obs = size(y);
    vector[num_obs] theta;
    for(i in 1:num_obs)
      theta[i] = log_combo2_inv_link(mu[i], phi[i], fup[i], finiteA[i], finiteB[i]);
    // poisson density without normalization y!
    return dot_product(to_vector(y), theta) - sum(exp(theta));
  }
  int blrm_combo2_tte_rng(real mu, real phi, real fup, int finiteA, int finiteB) {
    return poisson_rng(exp(log_combo2_inv_link(mu, phi, fup, finiteA, finiteB)));
  }
"

blrm_combo2_tte_stanvar <- stanvar(scode = blrm_combo2_tte_stan_code, block = "functions")

print("Stan model defined")



## get unit information prior for log(lambda)...use the fact that a
## negative binomial with kappa=0 is equal to a Poisson
neg_bin_fisher_info <- function(log_lambda, kappa, exp_time) {
    lambda_eff <- exp_time * exp(log_lambda)
    lambda_eff / (1 + kappa * lambda_eff)
}


## model formula
blrm_tte_model_nonhierarchical  <-
  bf(num_toxicities | vreal(follow_up) + vint(finite_drugB, finite_drugA)  ~ interIdrugB + exp(slopeIdrugB) * std_drugB, nl=TRUE, loop=TRUE) +
  lf(interIdrugB ~ 1 ) +
  lf(slopeIdrugB ~ 1 ) +
  lf(phi ~ 1 + mo(ocycle) ) +
  blrm_combo2_tte


## Prior:
## -7.2 approx cloglog(0.09) - log(3*6*7) for drugB (assuming no further risk increase)
## -7   approx cloglog(0.11) - log(3*6*7) for drugA


blrm_tte_prior_nonhierarchical <-
  prior(normal(-7.2, 1), nlpar="interIdrugB") +
  prior(normal(0, log(4)/1.96), nlpar="slopeIdrugB") +
  prior(normal(-7, 0.5), class="Intercept", dpar="phi") +
  prior(normal(0, 0.5), class="b", dpar="phi")


hist_data_nonhierarchical_fake <- tibble(
  group_id = factor("trial", levels = c("trial")),
  drugA = 1,
  drugB = 0,
  cycle1_caltime_dt = ymd("2019-09-14"),
  cycle2_caltime_dt = ymd("2019-09-14") + cycle_length,
  cycle3_caltime_dt = ymd("2019-09-14") + 2*cycle_length,
  cens_caltime_dt   = ymd("2019-09-14") + 3*cycle_length,
  dropout_caltime_dt = as_date(Inf),
  event_caltime_dt  = as_date(Inf)              
)

hist_data_nonhierarchical_prepared <- prepare_data(hist_data_nonhierarchical_fake)
get_prior(blrm_tte_model_nonhierarchical, data=hist_data_nonhierarchical_prepared)
print("Prior defined")

#########################################################
## MODEL SETUP
#########################################################

if(FALSE) {
    ## can be used to clean the cached model
    file.remove("brms_combo2_tte_model.rds")
}


combo2_tte_stanmodel_nonhierarchical <- make_stancode(
  blrm_tte_model_nonhierarchical,
  data = hist_data_nonhierarchical_prepared,
  stanvars = blrm_combo2_tte_stanvar,
  prior=blrm_tte_prior_nonhierarchical
)

combo2_tte_standata_nonhierarchical <- make_standata(
  blrm_tte_model_nonhierarchical,
  data = hist_data_nonhierarchical_prepared,
  stanvars = blrm_combo2_tte_stanvar,
  prior = blrm_tte_prior_nonhierarchical
)

cat(combo2_tte_stanmodel_nonhierarchical, file="blrm_combo2_tte_nonhierarchical.stan")



names(combo2_tte_stanmodel_nonhierarchical)

print("Make_stancode / make_standata done")

## setup model such that we cache the compiled stan model, no actual fit yet.
combo2_tte_model_nonhierarchical <- brm(
  blrm_tte_model_nonhierarchical,
  data = hist_data_nonhierarchical_prepared %>% mutate(follow_up = 1e-8),
  stanvars = blrm_combo2_tte_stanvar,
  prior=blrm_tte_prior_nonhierarchical,
  chains=0
)

stancode(combo2_tte_model_nonhierarchical )

print("brm call done")

print("expose function call done")

hist_combo2_tte_fit_nonhierarchical <- blrm_combo2_tte_fit_nonhierarchical(
  combo2_tte_model_nonhierarchical, hist_data_nonhierarchical_fake[0, ]
)

print("hist data fit done")



print("Model file loaded")
