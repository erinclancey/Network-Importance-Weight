

############# pMCMC Setup ############# 
rw.var <- matrix(
  c(0.0012,  0,        0,     # Dropped from 0.0020
    0,       0.00007,  0,     # Dropped from 0.00011
    0,       0,        0.050),# Dropped from 0.085
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
  target = 0.234,
  max.scaling = 50
)


##############Data Frame for Sim Study####################
sim_study <- data.frame(matrix(ncol = 15, nrow = nrow(pars)))
x <- c("w_true","w_mode","w_low","w_hi",
       "beta_par_true","beta_par_mode","beta_par_low","beta_par_hi",
       "k_true","k_mode","k_low","k_hi",
       "Gel_w", "Gel_beta","Gel_k")
colnames(sim_study) <- x

##########################Run Sim Study
for (i in 1:nrow(pars)){
  #prepare simulated data
  sim_dat <- data_list[[i]] %>% select(day, paste0("reports", 1:5))
  sim_dat %>%
    pomp(
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
    ) -> pomp_dat
  
  #set true parameter values
  sim_study$w_true[i] <- pars$w[i]
  sim_study$beta_par_true[i] <- pars$beta_par[i]
  sim_study$k_true[i] <- pars$k[i]

  # Set Up Chains
  nchains=5
  w = runif(nchains, 0.4,0.6)
  beta_par = runif(nchains, 0.00002*10000, 0.00003*10000)
  gamma = rep(pars["gamma"][i,], nchains)
  rho = rep(pars["rho"][i,], nchains)
  k = runif(nchains, 1,2)
  A_frame =  as.data.frame(matrix(rep(A_vec, 5), nrow = 5, byrow = TRUE))
  colnames(A_frame) <- A_names
  theta.start <- data.frame(cbind(w,beta_par,gamma,rho,k, A_frame))
  
  
  M=5000
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
  
  list_results_pmcmc <- results_pmcmc[c(seq(1:nrow(theta.start)))]
  chain.no <- seq(1:nchains)
  iter.no <- seq(1, M+1, by=1)
  
  for(j in seq_along(list_results_pmcmc)){
    list_results_pmcmc[[j]]$chain <- rep(chain.no[j],nrow(list_results_pmcmc[[j]]))
    list_results_pmcmc[[j]]$iter <- iter.no
  }
  
  posterior <- do.call(rbind, list_results_pmcmc)
  
  accepts <- posterior %>%
    dplyr::group_by(chain) %>%
    dplyr::summarise(
      # A step is accepted if ANY of the estimated parameters changed value
      accept_rate = mean(diff(w) != 0 | diff(beta_par) != 0 | diff(k) != 0)
    )
  
  # 2. Identify the chains that meet your 0.09 threshold using the new metric
  good_chains <- accepts %>%
    filter(accept_rate >= 0.09) %>%
    pull(chain) 
  
  if (length(good_chains) <2) {
    next # Skips to next iteration if good_chains is empty
  }
  
  # 3. Subset your posterior data to keep only those good chains
  posterior_filtered <- posterior %>%
    filter(chain %in% good_chains)
  
  # Select both parameters plus chain and iter
  post <- posterior_filtered  %>% 
    dplyr::select(w, beta_par, k, chain, iter)
  # Split by chain
  post.list <- split(post, f = post$chain)  # converts the dataframes back into a list for post-processing
  mcmc.list <- mcmc.list(list())
### Fills the list with MCMC objects
  ##########need to chagne to j
  for(j in seq_along(post.list)){
  mcmc.list[[j]] <- mcmc(post.list[[j]])
    }

  gel <- gelman.diag(mcmc.list, confidence = 0.95, transform = FALSE, autoburnin = TRUE, multivariate = FALSE)
  sim_study$Gel_w[i] <- gel$psrf[1,1]
  sim_study$Gel_beta[i] <- gel$psrf[2,1]
  sim_study$Gel_k[i] <- gel$psrf[3,1]
  
  # Process window
  processed <- window(mcmc.list, start=1000, end=M+1, thin=1) 
  processed <- data.frame(do.call(rbind, processed))
  
  # Posterior summaries for w
  sim_study$w_mode[i] <- as.vector(posterior.mode(mcmc(processed$w), adjust=1))
  sim_study$w_low[i]  <- ci(processed$w, ci=0.95, method="HDI")$CI_low
  sim_study$w_hi[i]   <- ci(processed$w, ci=0.95, method="HDI")$CI_high
  
  # Posterior summaries for beta_par
  sim_study$beta_par_mode[i] <- as.vector(posterior.mode(mcmc(processed$beta_par), adjust=1))
  sim_study$beta_par_low[i]  <- ci(processed$beta_par, ci=0.95, method="HDI")$CI_low
  sim_study$beta_par_hi[i]   <- ci(processed$beta_par, ci=0.95, method="HDI")$CI_high
  
  # Posterior summaries for w
  sim_study$k_mode[i] <- as.vector(posterior.mode(mcmc(processed$k), adjust=1))
  sim_study$k_low[i]  <- ci(processed$k, ci=0.95, method="HDI")$CI_low
  sim_study$k_hi[i]   <- ci(processed$k, ci=0.95, method="HDI")$CI_high
  
  print(i)
}  

write.csv(sim_study, file="network_simstudy_June28.csv")



# sim_df <- sim_study
# sim_df <- rbind(sim_study, sim_df_apr2)
# sim_df <- subset(sim_df, Gel_w<1.05 & Gel_beta<1.05 & Gel_k<1.05)


sim_df <- sim_df %>%
  mutate(
    w_inHDI = ifelse(w_true >= w_low & w_true <= w_hi, 1, 0),
    beta_inHDI = ifelse(beta_par_true >= beta_par_low & beta_par_true <= beta_par_hi, 1, 0),
    k_inHDI = ifelse(k_true >= k_low & k_true <= k_hi, 1, 0)
  )

View(sim_df)
nrow(sim_df)
mean(sim_df$w_true >= sim_df$w_low & sim_df$w_true <= sim_df$w_hi)
mean(sim_df$beta_par_true >= sim_df$beta_par_low & sim_df$beta_par_true <= sim_df$beta_par_hi)
mean(sim_df$k_true >= sim_df$k_low & sim_df$k_true <= sim_df$k_hi)

#####SET 1
# Make linear prediction for plots
wmod <- summary(lm(w_mode ~ w_true , data = sim_df))
betamod <- summary(lm(beta_par_mode ~ beta_par_true , data = sim_df))
kmod <- summary(lm(k_mode ~ k_true , data = sim_df))




library(ggplot2)
library(latex2exp)
library(patchwork)
library(scales)

## ---- Plot 1: w recovery ----
p1 <- ggplot(sim_df, aes(x = w_true, y = w_mode)) +
  geom_point(color = "darkorange2", shape = 20, size = 3, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, 
              linetype = 2, color = "black", size = 1) +
  geom_abline(intercept = wmod$coefficients[1,1], 
              slope = wmod$coefficients[2,1],
              linetype = 1, color = "darkorange4", size = 1) +
  xlim(0, 1) + ylim(0, 1) +
  labs(x = TeX("$w$"), y = TeX("$\\hat{w}$")) +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 16))



## ---- Plot 2: beta recovery ----

sci_hybrid <- function(x) {
  labs <- label_scientific()(x)   # scientific for nonzero
  labs[!is.na(x) & x == 0] <- "0.00"   # override zero
  labs
}

p2 <- ggplot(sim_df, aes(x = beta_par_true / 10000,
                         y = beta_par_mode / 10000)) +
  geom_point(color = "steelblue", shape = 20, size = 3, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1,
              linetype = 2, color = "black", size = 1) +
  geom_abline(intercept = betamod$coefficients[1,1] / 10000,
              slope = betamod$coefficients[2,1],
              linetype = 1, color = "steelblue4", size = 1) +
  scale_x_continuous(
    limits = c(0, 0.8/10000),
    labels = sci_hybrid
  ) +
  scale_y_continuous(
    limits = c(0, 0.8/10000),
    labels = sci_hybrid
  ) +
  labs(x = TeX("$\\beta$"),
       y = TeX("$\\hat{\\beta}$")) +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 16))

## ---- Plot 1: w recovery ----
p3 <- ggplot(sim_df, aes(x = k_true, y = k_mode)) +
  geom_point(color = "darkolivegreen4", shape = 20, size = 3, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, 
              linetype = 2, color = "black", size = 1) +
  geom_abline(intercept = kmod$coefficients[1,1], 
              slope = kmod$coefficients[2,1],
              linetype = 1, color = "darkolivegreen", size = 1) +
  xlim(0, 10) + ylim(0, 10) +
  labs(x = TeX("$k$"), y = TeX("$\\hat{k}$")) +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 16))


## ---- Combine side-by-side ----
p1_tagged <- p1 + labs(tag = "A") +
  theme(plot.tag = element_text(size = 16, face = "bold"))

p2_tagged <- p2 + labs(tag = "B") +
  theme(plot.tag = element_text(size = 16, face = "bold"))

p3_tagged <- p3 + labs(tag = "C") +
  theme(plot.tag = element_text(size = 16, face = "bold"))

p1_tagged + p2_tagged + p3_tagged







