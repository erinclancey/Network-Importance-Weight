setwd("~/BRUCE-CHAIN-main/Network MS/Sim_Study_FINAL_Sep2026")
source("SETUP Mech Model_July.R")
registerDoParallel()
registerDoRNG(2488620)

setwd("~/BRUCE-CHAIN-main/Network MS/Sim_Study_FINAL_Sep2026/50_datasets")
data_list <- lapply(
  1:50,
  function(i) read.csv(paste0("file_", i, ".csv"))
)
names(data_list) <- paste0("file_", 1:50)
setwd("~/BRUCE-CHAIN-main/Network MS/Sim_Study_FINAL_Sep2026")
pars <- read.csv(file="metadata_50_datasets.csv", header=TRUE)
oldsim <- read.csv(file="simstudy_original_50.csv", header=TRUE)
bad <- subset(oldsim, oldsim$Gel_w>1.05 | oldsim$Gel_beta>1.05 | oldsim$Gel_k>1.05)
bad$logLik <- rep(NA, nrow(bad))
obs_names <- paste0("reports", 1:5)
max_lik <- read.csv(file="iter_sim_study_Fullcomparison.csv", header = TRUE)

rw.var <- matrix(
  c(0.01, 0, 0,
    0, 0.001, 0,
    0, 0, 0.5),
  nrow = 3,
  dimnames = list(
    c("w", "beta_par", "k"),
    c("w", "beta_par", "k")
  )
)

proposal <- mvn_rw_adaptive(
  rw.var=rw.var,
  scale.start = 1,
  scale.cooling = 0.999,
  shape.start = 500,
  target = 0.32,
  max.scaling = 100
)

#############################################################
 for(i in 1:nrow(bad)){
   badindex=bad$X[i]
sim_dat <- data_list[[badindex]] %>% select(day, paste0("reports", 1:5))

sim_dat %>%select(day,all_of(obs_names))%>%
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
  )-> pomp_dat

# Set Up Chains
nchains=5
w = rep(max_lik$w_hat[badindex], nchains)
beta_par = rep(max_lik$beta_par_hat[badindex], nchains)
gamma = rep(pars["gamma"][badindex,], nchains)
rho = rep(pars["rho"][badindex,], nchains)
k = rep(max_lik$k_hat[badindex], nchains)
A_frame =  pars[badindex, 7:31]
rownames(A_frame) <- NULL
theta.start <- data.frame(cbind(w,beta_par,gamma,rho,k, A_frame))

M=5000
start_time <- Sys.time()
foreach (theta.start=iter(theta.start,"row"), .inorder=FALSE, .options.RNG = 2488620) %dopar% {
  library(pomp)
  library(magrittr)
  pomp_dat %>% pmcmc(Nmcmc=M,
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

for(j in seq_along(list_results_pmcmc)){
  list_results_pmcmc[[j]]$chain <- rep(chain.no[j],nrow(list_results_pmcmc[[j]]))
  list_results_pmcmc[[j]]$iter <- iter.no
}

posterior <- do.call(rbind, list_results_pmcmc)

post_name <- paste0("post_new_", i, ".csv")
write.csv(posterior, file = post_name, row.names = FALSE)

##################
post <- posterior %>% 
  dplyr::select(w, beta_par, k, chain, iter)
post.list <- split(post, f = post$chain)  # converts the dataframes back into a list for post-processing
mcmc.list <- mcmc.list(list())
### Fills the list with MCMC objects
for(k in seq_along(post.list)){
  mcmc.list[[k]] <- mcmc(post.list[[k]])
}

gel <- gelman.diag(mcmc.list, confidence = 0.95, transform = FALSE, autoburnin = TRUE, multivariate = FALSE)
bad$Gel_w[i] <- gel$psrf[1,1]
bad$Gel_beta[i] <- gel$psrf[2,1]
bad$Gel_k[i] <- gel$psrf[3,1]

#Post-process the chains
processed <- window(mcmc.list, start=2000, end=M+1, thin=1)
processed <- data.frame(do.call(rbind, processed))
processed.long <- processed %>%
  pivot_longer(cols = c(w, beta_par,k),
               names_to = "variable",
               values_to = "value")

# Posterior summaries for w
bad$w_mode[i] <- as.vector(posterior.mode(mcmc(processed$w), adjust=1))
bad$w_low_H[i]  <- ci(processed$w, ci=0.95, method="HDI")$CI_low
bad$w_hi_H[i]   <- ci(processed$w, ci=0.95, method="HDI")$CI_high

bad$w_mean[i] <- mean(processed$w)
bad$w_low_E[i]  <- ci(processed$w, ci=0.95, method="ETI")$CI_low
bad$w_hi_E[i]   <- ci(processed$w, ci=0.95, method="ETI")$CI_high

# Posterior summaries for beta_par
bad$beta_par_mode[i] <- as.vector(posterior.mode(mcmc(processed$beta_par), adjust=1))
bad$beta_par_low_H[i]  <- ci(processed$beta_par, ci=0.95, method="HDI")$CI_low
bad$beta_par_hi_H[i]   <- ci(processed$beta_par, ci=0.95, method="HDI")$CI_high

bad$beta_par_mean[i] <- mean(processed$beta_par)
bad$beta_par_low_E[i]  <- ci(processed$beta_par, ci=0.95, method="ETI")$CI_low
bad$beta_par_hi_E[i]   <- ci(processed$beta_par, ci=0.95, method="ETI")$CI_high

# Posterior summaries for w
bad$k_mode[i] <- as.vector(posterior.mode(mcmc(processed$k), adjust=1))
bad$k_low_H[i]  <- ci(processed$k, ci=0.95, method="HDI")$CI_low
bad$k_hi_H[i]   <- ci(processed$k, ci=0.95, method="HDI")$CI_high

bad$k_mean[i] <- mean(processed$k)
bad$k_low_E[i]  <- ci(processed$k, ci=0.95, method="ETI")$CI_low
bad$k_hi_E[i]   <- ci(processed$k, ci=0.95, method="ETI")$CI_high

mode.loglik <- data.frame(cbind(w=bad$w_mode[i],
                                beta_par=bad$beta_par_mode[i],
                                gamma=pars["gamma"][badindex,],
                                rho=pars["rho"][badindex,],
                                k=bad$k_mode[i],
                                pars[badindex, 7:31]))
pomp_dat %>%pfilter(params=mode.loglik,Np=500) -> pf
bad$logLik[i] <-  logLik(pf)

print(i)

 }

write.csv(bad, file="simstudy_new_22rows.csv")

