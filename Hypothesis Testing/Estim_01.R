setwd("~/BRUCE-CHAIN-main/Network MS/July/Hypothesis Testing")
source("SETUP Mech Model_FINAL.R")
registerDoParallel()
registerDoRNG(2488220)

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
  w=0.1
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


##############Set Up Estimation

# Set Up Chains
nchains=5
w = runif(nchains, 0.4,0.6)
beta_par = runif(nchains, 0.00003*10000, 0.00004*10000)
gamma = rep(pars["gamma"], nchains)
rho = rep(pars["rho"], nchains)
k = runif(nchains, 3,5)
A_frame =  as.data.frame(matrix(rep(A_vec, 5), nrow = 5, byrow = TRUE))
colnames(A_frame) <- A_names
theta.start <- data.frame(cbind(w,beta_par,gamma,rho,k, A_frame))

#############################################################

M=10000 # the number of mcmc iterations to run
start_time <- Sys.time()
foreach (theta.start=iter(theta.start,"row"), .inorder=FALSE) %dopar% {
  library(pomp)
  library(magrittr)
  sim_dat %>% pmcmc(Nmcmc=M,
                    proposal=proposal,
                    Np=1000,
                    params=theta.start
  ) -> pmcmc
  results <- as.data.frame(traces(pmcmc))
} -> results_pmcmc
end_time <- Sys.time()
end_time - start_time 
################################################################
list_results_pmcmc <- results_pmcmc[c(seq(1:nrow(theta.start)))]
chain.no <- seq(1:nchains)
iter.no <- seq(1, M+1, by=1)

for(i in seq_along(list_results_pmcmc)){
  list_results_pmcmc[[i]]$chain <- rep(chain.no[i],nrow(list_results_pmcmc[[i]]))
  list_results_pmcmc[[i]]$iter <- iter.no
}

posterior <- do.call(rbind, list_results_pmcmc)

#write.csv(posterior, file="estim_01.csv")

