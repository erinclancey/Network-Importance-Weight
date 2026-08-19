setwd("~/BRUCE-CHAIN-main/Network MS/July")
source("SETUP Mech Model_July.R")
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
  w=0
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

slice_design(
  center=coef(make_prof),
  w=rep(seq(from=0,to=1,length=100),each=3),
  beta_par=rep(seq(from=0.00001*10000,to=0.0001*10000,length=100),each=3),
  k=rep(seq(from=0.5,to=10,length=100),each=3)
) -> p

foreach (theta=iter(p,"row"), .combine=rbind,
         .inorder=FALSE) %dopar%
  {
    library(pomp)
    library(dplyr)
    make_prof %>% pfilter(params=theta,Np=50) -> pf
    theta$loglik <- logLik(pf)
    theta
  } -> p


p <- subset(p, loglik!=-Inf)


p %>% 
  gather(variable, value, w, beta_par, k) %>%
  filter(variable == slice) %>%
  mutate(
    variable = fct_relevel(variable, "w", "beta_par", "k"),
    variable = recode(variable,
                      "beta_par" = "β",
                      "w" = "w",
                      "k" = "k")
  ) %>%
  ggplot(aes(x = value, y = loglik, color = variable)) +
  geom_point() +
  facet_wrap(~ variable, scales = "free") +
  guides(color = "none") +
  labs(
    x = "parameter value",
    y = "log-likelihood",
    color = ""
  ) + geom_vline(data = data.frame(variable = c("w","β","k"),
                                   x0 = c(w, beta_par, k)),
                 aes(xintercept = x0),
                 color = "black",
                 inherit.aes = FALSE)+
  ggtitle("Profile Likelihood- Network 30% Connected")


# 2, 3, 4, 5 should not work

##############Set Up Estimation

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
#Time difference of 3.313588 hours
################################################################
list_results_pmcmc <- results_pmcmc[c(seq(1:nrow(theta.start)))]
chain.no <- seq(1:nchains)
iter.no <- seq(1, M+1, by=1)

for(i in seq_along(list_results_pmcmc)){
  list_results_pmcmc[[i]]$chain <- rep(chain.no[i],nrow(list_results_pmcmc[[i]]))
  list_results_pmcmc[[i]]$iter <- iter.no
}

posterior <- do.call(rbind, list_results_pmcmc)

#write.csv(posterior, file="estim_0_Aug11.csv")

posterior <- read.csv(file="estim_0_Aug11.csv", header=TRUE)
posterior$beta_par <- posterior$beta_par/10000

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
for(i in seq_along(post.list)){
  mcmc.list[[i]] <- mcmc(post.list[[i]])
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
w_mode <- as.vector(posterior.mode(mcmc(processed$w), adjust=2))
w_low  <- ci(processed$w, ci=0.95, method="HDI")$CI_low
w_hi   <- ci(processed$w, ci=0.95, method="HDI")$CI_high

# Posterior summaries for beta_par
beta_par_mode <- as.vector(posterior.mode(mcmc(processed$beta_par), adjust=2))
beta_par_low  <- ci(processed$beta_par, ci=0.95, method="HDI")$CI_low
beta_par_hi   <- ci(processed$beta_par, ci=0.95, method="HDI")$CI_high

k_mode <- as.vector(posterior.mode(mcmc(processed$k), adjust=2))
k_low  <- ci(processed$k, ci=0.95, method="HDI")$CI_low
k_hi   <- ci(processed$k, ci=0.95, method="HDI")$CI_high

library(bayestestR)
library(logspline)
prior_w <- runif(nrow(processed), min = 0, max = 1)
bf <- bayesfactor_parameters(processed$w, prior = prior_w, null = 0)
print(bf)



# Combine into a summary data frame
summaries <- data.frame(
  variable = c("w", "beta_par","k"),
  true     = c(w=0,
               beta_par=0.000025,
               k=3),
  mode     = c(w_mode, beta_par_mode, k_mode),
  low      = c(w_low, beta_par_low, k_low),
  high     = c(w_hi, beta_par_hi, k_hi)
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
  labs(title = "A", x = "Parameter Value", y = "Density") + 
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


cleaned_summary <- clean_summary_df(summaries)[,-2]

# View your newly formatted data
print(cleaned_summary)

# ggsave(
#   filename = "post_w0_aug.pdf",
#   width    = 9,
#   height   = 10,
#   units    = "in",
#   device   = "pdf"
# )

#######################

HDI_post <- posterior %>%
  filter(
    w >= as.numeric(cleaned_summary$low[1]) & w <= as.numeric(cleaned_summary$high[1]),
    beta_par >= as.numeric(cleaned_summary$low[2]) & beta_par <= as.numeric(cleaned_summary$high[2]),
    k >= as.numeric(cleaned_summary$low[3]) & k <= as.numeric(cleaned_summary$high[3]) )

# Sample 1000 rows without replacement
HDI_sample <- HDI_post %>% 
  slice_sample(n = 1000)

HDI_sample$Incidence_alt <- rep(NA,nrow(HDI_sample))
HDI_sample$Cases_alt <- rep(NA, nrow(HDI_sample))
HDI_sample$Incidence_null <- rep(NA,nrow(HDI_sample))
HDI_sample$Cases_null <- rep(NA, nrow(HDI_sample))

max_attempts <- 100  # Set maximum retries per scenario

for (i in 1:nrow(HDI_sample)) {
  w_alt <- HDI_sample$w[i]
  w_null <- 0
  beta_post <- HDI_sample$beta_par[i] * 10000
  k_post <- HDI_sample$k[i]
  
  # --- ALTERNATIVE SCENARIO ---
  pars_alt <- c(w_alt, beta_post, gamma, rho, k_post, A_vec)
  names(pars_alt) <- paramnames
  
  attempts <- 0
  sims_alt <- NULL
  
  repeat {
    attempts <- attempts + 1
    
    sim <- gen_data %>%
      simulate(
        params = pars_alt,
        nsim = 1,
        format = "data.frame",
        include.data = FALSE
      )
    
    nonzero_check <- colSums(sim[, paste0("reports", 1:5)]) > 0
    
    if (all(nonzero_check)) {
      sims_alt <- sim
      break
    }
    
    if (attempts >= max_attempts) {
      warning(paste("Alternative scenario reached max attempts at row", i))
      break
    }
  }
  
  if (!is.null(sims_alt)) {
    HDI_sample$Incidence_alt[i] <- sum(sims_alt[, c("H1", "H2", "H3", "H4", "H5")])
    HDI_sample$Cases_alt[i] <- sum(sims_alt[, paste0("reports", 1:5)])
  } else {
    HDI_sample$Incidence_alt[i] <- NA
    HDI_sample$Cases_alt[i] <- NA
  }
  
  # --- NULL SCENARIO ---
  pars_null <- c(w_null, beta_post, gamma, rho, k_post, A_vec)
  names(pars_null) <- paramnames
  
  attempts <- 0
  sims_null <- NULL
  
  repeat {
    attempts <- attempts + 1
    
    sim <- gen_data %>%
      simulate(
        params = pars_null,
        nsim = 1,
        format = "data.frame",
        include.data = FALSE
      )
    
    nonzero_check <- colSums(sim[, paste0("reports", 1:5)]) > 0
    
    if (all(nonzero_check)) {
      sims_null <- sim
      break
    }
    
    if (attempts >= max_attempts) {
      warning(paste("Null scenario reached max attempts at row", i))
      break
    }
  }
  
  if (!is.null(sims_null)) {
    HDI_sample$Incidence_null[i] <- sum(sims_null[, c("H1", "H2", "H3", "H4", "H5")])
    HDI_sample$Cases_null[i] <- sum(sims_null[, paste0("reports", 1:5)])
  } else {
    HDI_sample$Incidence_null[i] <- NA
    HDI_sample$Cases_null[i] <- NA
  }
}

# Calculate per-capita incidence (ignoring failed runs)
HDI_sample$attack_rate_alt <- HDI_sample$Incidence_alt / (5000 * 5)
HDI_sample$attack_rate_null <- HDI_sample$Incidence_null / (5000 * 5)
HDI_sample$diff_AR <- HDI_sample$attack_rate_null - HDI_sample$attack_rate_alt
hist(HDI_sample$diff_AR)
library(bayestestR)
# 95% High Density Interval of the paired difference
diff_hdi <- hdi(HDI_sample$diff_AR , ci = 0.95)
print(diff_hdi, digits = 8)


library(ggplot2)
library(bayestestR)

# 1. Compute HDI
diff_hdi <- hdi(HDI_sample$diff_AR, ci = 0.95)
ci_low <- diff_hdi$CI_low
ci_high <- diff_hdi$CI_high

# 2. Build summary dataframe for HDI shading and reference lines
diff_summary <- data.frame(
  low = ci_low,
  high = ci_high,
  zero_ref = 0
)

# 3. Plot matching your custom template
ggplot(HDI_sample, aes(x = diff_AR)) +
  theme_minimal() +
  
  # Histogram & Density overlay using density scaling
  geom_histogram(aes(y = after_stat(density)), fill = "#D55E00", color = "#D55E00", position = "identity", alpha = 0.2, bins = 50) +
  geom_density(fill = "#D55E00", color = "#D55E00", alpha = 0.2, adjust = 2) +
  
  # Shaded 95% HDI Rectangle
  geom_rect(
    data = diff_summary,
    aes(xmin = low, xmax = high, ymin = 0, ymax = Inf),
    inherit.aes = FALSE,
    alpha = 0.1,
    fill = "grey20"
  ) +
  
  # Vertical reference line at 0 (matching your dashed red line style)
  geom_vline(
    data = diff_summary,
    aes(xintercept = zero_ref),
    color = "red",
    linewidth = 0.75,
    linetype = 2
  ) +
  
  # Labels and Axis Formatting
  labs(
    title = "A",
    x = expression(Delta ~ "Attack Rate" ~ (H[A]-H[0])),
    y = "Density"
  ) +
  scale_x_continuous(n.breaks = 6) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)), n.breaks = 6) +
  
  # Theme and Font Sizing matching your original setup
  theme(
    title = element_text(size = 20, face = "bold"),
    strip.text = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = 15, color = "black"),
    panel.spacing = unit(0, "lines"),
    legend.position = "none"
  )
# 
# # 1. Reshape data from wide to long format
# df_long <- HDI_sample %>%
#   select(attack_rate_alt, attack_rate_null) %>%
#   pivot_longer(
#     cols = c(attack_rate_alt, attack_rate_null),
#     names_to = "variable",
#     values_to = "attack_rate"
#   )
# 
# # 2. Plot overlapping histograms
# ggplot(df_long, aes(x = attack_rate, fill = variable)) +
#   geom_histogram(alpha = 0.5, position = "identity", bins = 30, color = "black") +
#   scale_fill_manual(
#     values = c("attack_rate_alt" = "#0072B2", "attack_rate_null" = "#D55E00"),
#     labels = c("Alternative", "Null")
#   ) +
#   scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
#   labs(
#     x = "Attack Rate",
#     y = "Frequency",
#     fill = "Model Scenario",
#     title = "Distribution of Attack Rates (Alt vs. Null)"
#   ) +
#   theme_minimal(base_size = 14)
# 
# 
# # 1. Reshape data from wide to long format
# df_long <- HDI_sample %>%
#   select(Cases_alt, Cases_null) %>%
#   pivot_longer(
#     cols = c(Cases_alt, Cases_null),
#     names_to = "variable",
#     values_to = "cases"
#   )
# 
# # 2. Plot overlapping histograms
# ggplot(df_long, aes(x = cases, fill = variable)) +
#   geom_histogram(alpha = 0.5, position = "identity", bins = 30, color = "black") +
#   scale_fill_manual(
#     values = c("Cases_alt" = "#0072B2", "Cases_null" = "#D55E00"),
#     labels = c("Alternative", "Null")
#   ) +
#   scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
#   labs(
#     x = "Total Cases",
#     y = "Frequency",
#     fill = "Model Scenario",
#     title = "Total Cases (Alt vs. Null)"
#   ) +
#   theme_minimal(base_size = 14)
# 
