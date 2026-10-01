setwd("~/BRUCE-CHAIN-main/Network MS/Sim_Study_FINAL_Sep2026")

library(bayestestR)
library(logspline)

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



##############Data Frame for Sim Study####################
sim_study <- data.frame(matrix(ncol = 25, nrow = nrow(pars)))
x <- c("w_true","w_mode","w_low_H","w_hi_H","w_mean","w_low_E","w_hi_E",
       "beta_par_true","beta_par_mode","beta_par_low_H","beta_par_hi_H",
       "beta_par_mean","beta_par_low_E","beta_par_hi_E",
       "k_true","k_mode","k_low_H","k_hi_H", "k_mean","k_low_E","k_hi_E",
       "BF", "Gel_w", "Gel_beta","Gel_k")
colnames(sim_study) <- x
#############################################################
 for(i in 1:n){
   
sim_dat <- data_list[[i]] %>% select(day, paste0("reports", 1:5))

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


#set true parameter values
sim_study$w_true[i] <- pars$w[i]
sim_study$beta_par_true[i] <- pars$beta_par[i]
sim_study$k_true[i] <- pars$k[i]

# Set Up Chains
nchains=5
w = runif(nchains, 0.4,0.6)
beta_par = runif(nchains, 0.00003*10000, 0.00004*10000)
gamma = rep(pars["gamma"][i,], nchains)
rho = rep(pars["rho"][i,], nchains)
k = runif(nchains, 3,5)
A_frame =  pars[i, 6:30]
rownames(A_frame) <- NULL
theta.start <- data.frame(cbind(w,beta_par,gamma,rho,k, A_frame))

M=5000
start_time <- Sys.time()
foreach (theta.start=iter(theta.start,"row"), .inorder=FALSE) %dopar% {
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
post_name <- paste0("post_", i, ".csv")
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
sim_study$Gel_w[i] <- gel$psrf[1,1]
sim_study$Gel_beta[i] <- gel$psrf[2,1]
sim_study$Gel_k[i] <- gel$psrf[3,1]



#Post-process the chains
processed <- window(mcmc.list, start=2000, end=M+1, thin=1)
processed <- data.frame(do.call(rbind, processed))

processed.long <- processed %>%
  pivot_longer(cols = c(w, beta_par,k),
               names_to = "variable",
               values_to = "value")


# Posterior summaries for w
sim_study$w_mode[i] <- as.vector(posterior.mode(mcmc(processed$w), adjust=1))
sim_study$w_low_H[i]  <- ci(processed$w, ci=0.95, method="HDI")$CI_low
sim_study$w_hi_H[i]   <- ci(processed$w, ci=0.95, method="HDI")$CI_high

sim_study$w_mean[i] <- mean(processed$w)
sim_study$w_low_E[i]  <- ci(processed$w, ci=0.95, method="ETI")$CI_low
sim_study$w_hi_E[i]   <- ci(processed$w, ci=0.95, method="ETI")$CI_high

# Posterior summaries for beta_par
sim_study$beta_par_mode[i] <- as.vector(posterior.mode(mcmc(processed$beta_par), adjust=1))
sim_study$beta_par_low_H[i]  <- ci(processed$beta_par, ci=0.95, method="HDI")$CI_low
sim_study$beta_par_hi_H[i]   <- ci(processed$beta_par, ci=0.95, method="HDI")$CI_high

sim_study$beta_par_mean[i] <- mean(processed$beta_par)
sim_study$beta_par_low_E[i]  <- ci(processed$beta_par, ci=0.95, method="ETI")$CI_low
sim_study$beta_par_hi_E[i]   <- ci(processed$beta_par, ci=0.95, method="ETI")$CI_high

# Posterior summaries for w
sim_study$k_mode[i] <- as.vector(posterior.mode(mcmc(processed$k), adjust=1))
sim_study$k_low_H[i]  <- ci(processed$k, ci=0.95, method="HDI")$CI_low
sim_study$k_hi_H[i]   <- ci(processed$k, ci=0.95, method="HDI")$CI_high

sim_study$k_mean[i] <- mean(processed$k)
sim_study$k_low_E[i]  <- ci(processed$k, ci=0.95, method="ETI")$CI_low
sim_study$k_hi_E[i]   <- ci(processed$k, ci=0.95, method="ETI")$CI_high


print(i)
}


write.csv(sim_study, file="Sim_Study_July_original_50.csv")




