setwd("~/BRUCE-CHAIN-main/Network MS")
source("Network Mech Model Pomp_SETUP_prior.R")
registerDoParallel()
registerDoRNG(2488420)

#Make One dataset
day <- seq(1, 200, 1)
Data <- as.data.frame(day)
obs_names <- paste0("reports", 1:5)

repeat {
  # Matrix A generated then hard coded into step
  p=0.3
  A_mat <- matrix(rbinom(25,1,p), nrow=5)
  diag(A_mat) <- 1
  A_mat[lower.tri(A_mat)] <- t(A_mat)[lower.tri(A_mat)]
  #A_mat == t(A_mat)   # should be all TRUE
  A_vec <- as.vector(A_mat)
  
  
  # True Parameter Values
  w=0.5
  beta_par=0.000025*10000
  gamma=1/10
  rho=0.5
  k=3
  pars <- c(w, beta_par, gamma, rho, k, A_vec)
  names(pars) <- paramnames
  
  Data %>% pomp(
    times="day",t0=0,
    rprocess=euler(step,delta.t=tau),
    rinit=rinit,
    accumvars=accumvars,
    statenames=statenames,
    paramnames=paramnames,
    rmeasure=rmeas,
    dmeasure=dmeas,
    rprior = prior_sampler,
    dprior = prior_density,
    obsnames=obs_names
  ) -> gen_data
  
  gen_data %>%
    simulate(
      params=pars,
      nsim=1,
      format="data.frame",include.data=FALSE) -> sims
  
  # Check if all reports1–reports5 have at least one non-zero value
  nonzero_check <- colSums(sims[ , paste0("reports", 1:5)]) > 0
  
  if (all(nonzero_check)) {
    break
  }
}

###################Check Slice of the Likelihood
sims %>%select(day,all_of(obs_names))%>%
  pomp(
    times      = "day",
    t0         = 0,
    rprocess   = euler(step, delta.t = tau),
    rinit      = rinit,
    statenames = statenames,
    accumvars  = accumvars,
    paramnames = paramnames,
    rmeasure   = rmeas,
    dmeasure   = dmeas,
    rprior = prior_sampler,
    dprior = prior_density,
    obsnames=obs_names
  )-> sim_dat

sim_dat %>%
  pomp(params=pars) -> make_prof


#################Diagnostics
# --- RUN THIS AFTER CREATING 'make_prof' ---

# 1. Define configuration
test_np <- 1000    # The particle count you want to test
n_reps  <- 15    # 15 reps gives a highly stable variance estimate

cat("Running particle filter diagnostics... please wait.\n")

# 2. Replicate the filters at the baseline/center parameters
pf_reps <- replicate(n_reps, {
  make_prof %>% pfilter(Np = test_np)
}, simplify = FALSE)

# 3. Calculate metrics
loglik_values   <- sapply(pf_reps, logLik)
loglik_variance <- var(loglik_values)
mean_loglik     <- mean(loglik_values)

# Pull ESS from the first run across the time series
ess_time_series <- eff_sample_size(pf_reps[[1]])
min_ess         <- min(ess_time_series)
avg_ess         <- mean(ess_time_series)

# 4. Print clean report to console
cat("\n=== POMP PARTICLE FILTER DIAGNOSTIC REPORT ===\n")
cat(sprintf("Tested Particle Count (Np) : %d\n", test_np))
cat(sprintf("Replications Evaluated     : %d\n", n_reps))
cat("----------------------------------------------\n")
cat(sprintf("Mean Log-Likelihood        : %.2f\n", mean_loglik))
cat(sprintf("Log-Likelihood Variance    : %.4f\n", loglik_variance))
cat(sprintf("Minimum ESS Found          : %.1f\n", min_ess))
cat(sprintf("Average ESS Across Time    : %.1f\n", avg_ess))
cat("==============================================\n\n")

# 5. Diagnostic Guidance Interpretations
if (loglik_variance > 3.0) {
  cat("⚠️ WARNING: Log-Likelihood Variance is > 3.0. Your PMCMC chain will likely lock up. Increase Np.\n")
} else if (loglik_variance < 0.5) {
  cat("💡 NOTE: Log-Likelihood Variance is very low (< 0.5). You might be able to lower Np to save time.\n")
} else {
  cat("✅ SUCCESS: Log-Likelihood Variance is in the optimal 0.5 - 3.0 target window for PMCMC.\n")
}

if (min_ess < 10) {
  cat("⚠️ WARNING: Particle depletion detected! Minimum ESS dropped below 10. Consider increasing Np.\n")
}

