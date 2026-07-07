
# standard.devs <-c(0.01, 0.001, 0.05)
# params.est <- c("w", "beta_par", "k")
# rw.sd <- setNames(standard.devs, params.est)
# # proposal for the first iteration
# proposal <- mvn_diag_rw(rw.sd)
# 

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
sim_study <- data.frame(matrix(ncol = 22, nrow = nrow(pars)))
x <- c("w_true","w_mode","w_low_H","w_hi_H","w_mean","w_low_E","w_hi_E",
       "beta_par_true","beta_par_mode","beta_par_low_H","beta_par_hi_H",
       "beta_par_mean","beta_par_low_E","beta_par_hi_E",
       "k_true","k_mode","k_low_H","k_hi_H", "k_mean","k_low_E","k_hi_E",
       "BF")
colnames(sim_study) <- x
#############################################################
i=2

# 2. Build your pomp object
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

pomp_dat %>%pfilter(params=theta.start[5,],Np=50) -> pf
logLik(pf)

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
# post_name <- paste0("post_", i, ".csv")
# write.csv(posterior, file = post_name, row.names = FALSE)

param_cov <- cov(posterior[, c("w", "beta_par", "k")])

###########################################################  
accepts <- posterior %>%
  dplyr::group_by(chain) %>%
  dplyr::summarise(
    # A step is accepted if ANY of the estimated parameters changed value
    accept_rate = mean(diff(w) != 0 | diff(beta_par) != 0 | diff(k) != 0)
  )
accepts
####################################################

### Post-process the chains

# Select the parameters of interest
post <- posterior %>% 
  dplyr::select(w, beta_par,k, chain, iter)

# Ensure chain is treated as character
chains <- post %>% mutate(chain = as.character(chain))

# Pivot both w and beta_par into long format
chains.long <- chains %>%
  pivot_longer(cols = c(w, beta_par,k), 
               names_to = "variable", 
               values_to = "value")


plot_names <- as_labeller(c(
  'w' = "paste(hat(w))",
  'beta_par' = "paste(hat(beta))",
  'k' = "paste(hat(k))"
), label_parsed)


# Plot both parameters in facets
ggplot(chains.long, aes(x = iter, y = value, group = chain)) + 
  geom_line(aes(color = chain)) + 
  theme_minimal() +
  #scale_y_continuous(limits = c(c(0, 1),c(0,1),c(0,10)))+
  scale_color_manual(values = c("#F0E442","#000000", "#009E73","#999999","#0072B2")) + 
  facet_wrap(vars(variable), labeller = plot_names, scales = 'free', ncol = 1) +
  labs(title = "",
       x = "Iteration", y = "Parameter Value") + 
  theme(
    title = element_text(size = 18),
    strip.text = element_text(size = 20),
    axis.title = element_text(size = 18),
    panel.spacing = unit(0, "lines"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 'cm')
  )
##################
post.list <- split(post, f = post$chain)  # converts the dataframes back into a list for post-processing
mcmc.list <- mcmc.list(list())
### Fills the list with MCMC objects
for(k in seq_along(post.list)){
  mcmc.list[[k]] <- mcmc(post.list[[k]])
}

effectiveSize(mcmc.list)
#needs to be over 200

### Perform Diagnistic Tests for convergence and burn-in
raftery.diag(mcmc.list, q=0.025, r=0.005, s=0.95, converge.eps=0.001)
geweke.diag(mcmc.list, frac1=0.5, frac2=0.25)
gelman.diag(mcmc.list, confidence = 0.95, transform = FALSE, autoburnin = TRUE, multivariate = FALSE)

#Post-process the chains
processed <- window(mcmc.list, start=2000, end=M+1, thin=1)
processed <- data.frame(do.call(rbind, processed))

processed.long <- processed %>%
  pivot_longer(cols = c(w, beta_par,k),
               names_to = "variable",
               values_to = "value")

nrow(processed)

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

library(bayestestR)
library(logspline)
prior_w <- runif(nrow(processed), min = 0, max = 1)
bf <- bayesfactor_parameters(processed$w, prior = prior_w, null = 0)
sim_study$BF[i] <- exp(bf$log_BF)



# Combine into a summary data frame
summaries <- data.frame(
  variable = c("w", "beta_par","k"),
  true     = c(w=pars$w[i],
               beta_par=pars$beta_par[i],
               k=pars$k[i]),
  mode     = c(sim_study$w_mode[i], sim_study$beta_par_mode[i],sim_study$k_mode[i]),
  low      = c(sim_study$w_low_H[i], sim_study$beta_par_low_H[i], sim_study$k_low_H[i] ),
  high     = c(sim_study$w_hi_H[i] , sim_study$beta_par_hi_H[i], sim_study$k_hi_H[i] )
)

summaries <- summaries %>% 
  mutate(across(where(is.numeric), ~ signif(.x, digits = 3)))




# 1. Dynamically build levels: w first, then beta_par, then k, then any leftovers
existing_vars <- unique(as.character(processed.long$variable))
lvl <- c("w", "beta_par", "k")
lvl <- c(lvl, setdiff(existing_vars, lvl)) # Adds any other variables at the end

processed.long$variable <- factor(processed.long$variable, levels = lvl)
summaries$variable <- factor(summaries$variable, levels = lvl)

# 1. Create a dummy data frame to force the x-axis limits for 'w'
blank_data <- data.frame(variable = "w", value = c(0, 1))
blank_data$variable <- factor(blank_data$variable, levels = levels(processed.long$variable))

# 2. Update the plot
P <- ggplot(processed.long, aes(x = value, fill = variable, color = variable)) + 
  theme_minimal() + 
  # This hidden layer forces the x-axis for 'w' to 0 and 1
  geom_blank(data = blank_data) + 
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.2, bins = 50) + 
  geom_density(alpha = .2, adjust = 2) + 
  geom_vline(data = summaries, aes(xintercept = mode), color = "black", linewidth = 0.75, linetype = 2) + 
  geom_vline(data = summaries, aes(xintercept = true), color = "red", linewidth = 0.75, linetype = 2) + 
  geom_rect(data = summaries, aes(xmin = low, xmax = high, ymin = 0, ymax = Inf), 
            inherit.aes = FALSE, alpha = 0.1, fill = "grey20") + 
  facet_wrap(vars(variable), labeller = plot_names, scales = 'free', ncol = 1) + 
  scale_fill_manual(values = c("w" = "#0072B2", "beta_par" = "grey40", "k" = "grey40")) + 
  scale_color_manual(values = c("w" = "#0072B2", "beta_par" = "grey40", "k" = "grey40")) + 
  labs(title = "B", x = "Parameter Value", y = "Density") + 
  # Note: removed 'expand' here or set it to 0 if you want the limits to be EXACTLY 0 and 1
  scale_x_continuous(n.breaks = 6) + 
  scale_y_continuous(expand = expansion(mult = 0.2), n.breaks = 6) + 
  theme(title = element_text(size = 20, face="bold"), 
        strip.text = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20), 
        axis.text = element_text(size = 15, color = "black"), 
        panel.spacing = unit(0, "lines"),
        legend.position = "none")

plot(P)

clean_summary_df <- function(summary_obj) {
  df <- as.data.frame(summary_obj)
  
  # 1. Ensure parameter names are explicitly captured as the first column
  if (!"Parameter" %in% colnames(df)) {
    df <- cbind(Parameter = rownames(df), df)
  }
  rownames(df) <- NULL
  
  # 2. Identify and format columns that contain numbers (even if stored as text)
  # We skip the text-based columns like "Parameter" and "variable"
  cols_to_convert <- c("true", "mode", "low", "high")
  
  for (col in cols_to_convert) {
    if (col %in% colnames(df)) {
      # Force text like "2.5e-05" into actual numbers
      numeric_values <- as.numeric(as.character(df[[col]]))
      
      # Format strictly as decimals with up to 6 decimal places, no scientific notation
      df[[col]] <- format(round(numeric_values, 10), scientific = FALSE, drop0trailing = TRUE)
    }
  }
  
  return(df)
}


cleaned_summary <- clean_summary_df(summaries)

# View your newly formatted data
print(cleaned_summary)

#write.csv(sim_study, file="manual_simstudy_July6.csv")



sim_df <- sim_study
sim_df <- sim_df %>%
  mutate(
    w_inHDI = ifelse(w_true >= w_low_H & w_true <= w_hi_H, 1, 0),
    beta_inHDI = ifelse(beta_par_true >= beta_par_low_H & beta_par_true <= beta_par_hi_H, 1, 0),
    k_inHDI = ifelse(k_true >= k_low_H & k_true <= k_hi_H, 1, 0)
  )

mean(sim_df$w_true >= sim_df$w_low_H & sim_df$w_true <= sim_df$w_hi_H)
mean(sim_df$beta_par_true >= sim_df$beta_par_low_H & sim_df$beta_par_true <= sim_df$beta_par_hi_H)
mean(sim_df$k_true >= sim_df$k_low_H & sim_df$k_true <= sim_df$k_hi_H)

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











