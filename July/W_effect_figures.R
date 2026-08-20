setwd("~/BRUCE-CHAIN-main/Network MS/July")
source("SETUP Mech Model_July.R")
p=0.3
A_mat <- matrix(rbinom(25,1,p), nrow=5)
diag(A_mat) <- 1
A_mat[lower.tri(A_mat)] <- t(A_mat)[lower.tri(A_mat)]
#A_mat == t(A_mat)   # should be all TRUE
A_vec <- as.vector(A_mat)

gamma=1/10

w_per <- seq(0, 1, length.out=1000)
R0_pop <- vector()
for (i in 1:length(w_per)){
  Cmat <- A_mat * w_per[i] + (1 - w_per[i])   # off-diagonal formula applied everywhere
  diag(Cmat) <- 1 
  # Compute eigen decomposition
  eig <- eigen(Cmat)
  # Eigenvalues
  dom_eig <- max(eig$values)
  R0_pop[i] <- 5000*0.00005 / gamma * dom_eig
}
# Combine into a data frame
df <- data.frame(w_per = w_per, R0_pop = R0_pop)

# Plot with ggplot2
ggplot(df, aes(x = w_per, y = R0_pop)) +
  geom_line(color = "#0072B2", linewidth = 1) +
  labs(
    x = expression(w[per]),
    y = expression(R[0] ~ "(pop)"),
    title = expression("Dominant Eigenvalue Impact on " ~ R[0])
  ) +
  theme_minimal(base_size = 14)

library(ggplot2)

w_per <- seq(0, 1, length.out = 1000)
beta_par <- c(0.000005, 0.000015, 0.000025, 0.00005, 0.00001)

# 1. Compute dominant eigenvalues once across w_per
dom_eig <- numeric(length(w_per))
for (i in seq_along(w_per)) {
  Cmat <- A_mat * w_per[i] + (1 - w_per[i])
  diag(Cmat) <- 1
  dom_eig[i] <- max(eigen(Cmat)$values)
}

# 2. Build grid for all combinations of w_per and beta
df <- expand.grid(w_per = w_per, beta = beta_par)
df$dom_eig <- rep(dom_eig, times = length(beta_par))

# 3. Calculate R0_pop across all combinations
df$R0_pop <- (5000 * df$beta / gamma) * df$dom_eig

# Install if needed: install.packages("wesanderson")
library(wesanderson)
library(ggsci)

# Warm cinematic palette
scale_color_manual(values = wes_palette("Darjeeling1", n = 5, type = "continuous"))

# Cool vintage palette
scale_color_manual(values = wes_palette("Zissou1", n = 5, type = "continuous"))

# 4. Plot with ggplot2
ggplot(df, aes(x = w_per, y = R0_pop, color = factor(beta))) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_y_continuous(
    limits = c(0, NA),
    breaks = seq(0, max(df$R0_pop, na.rm = TRUE), by = 2)
  ) +
  scale_color_manual(values = wes_palette("Zissou1", n = 5, type = "continuous")) +
  labs(
    x = expression(w),
    y = expression(R[0] ~ "(pop)"),
    color = expression(beta),
    title = expression("Impact of " ~ w ~ " on " ~ R[0] ~ " across " ~ beta)
  ) +
  theme_minimal(base_size = 14)



day <- seq(1, 200, 1)
Data <- as.data.frame(day)
obs_names <- paste0("reports", 1:5)
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
beta_par=0.000025*10000
gamma=1/10
rho=0.5
k=3

w_per <- seq(0, 1, length.out=1000)
total_cases <- vector()
total_incidence <- vector()
for (i in 1:length(w_per)){
  pars <- c(w_per[i], beta_par, gamma, rho, k, A_vec)
  names(pars) <- paramnames
  
  repeat {
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
  total_incidence[i] <- sum(sims[, c("H1","H2","H3","H4","H5")])
  total_cases[i] <- sum(sims[, c("reports1","reports2","reports3","reports4","reports5")])
}


df_inc <- data.frame(w_per, total_incidence, total_cases)
df_inc[1,2]*0.95

plot(w_per, total_incidence/(5000*5))  


library(ggplot2)
library(dplyr)

w_per <- seq(0, 1, length.out = 100)
beta_par <- c(0.000005, 0.000015, 0.000025, 0.00005, 0.00001)*10000

# Create grid of all combinations
df <- expand.grid(w_per = w_per, beta = beta_par)
df$total_incidence <- numeric(nrow(df))
df$total_cases <- numeric(nrow(df))

# Run simulation for each combination of w_per and beta
for (i in 1:nrow(df)) {
  current_w <- df$w_per[i]
  current_beta <- df$beta[i]
  
  # Format parameters
  pars <- c(current_w, current_beta, gamma, rho, k, A_vec)
  names(pars) <- paramnames
  
  repeat {
    sims <- gen_data %>%
      simulate(
        params = pars,
        nsim = 1,
        format = "data.frame",
        include.data = FALSE
      )
    
    # Check if all reports1–reports5 have at least one non-zero value
    nonzero_check <- colSums(sims[, paste0("reports", 1:5)]) > 0
    
    if (all(nonzero_check)) {
      break
    }
  }
  
  df$total_incidence[i] <- sum(sims[, c("H1", "H2", "H3", "H4", "H5")])
  df$total_cases[i] <- sum(sims[, paste0("reports", 1:5)])
}

# Calculate per-capita incidence
df$incidence_per_capita <- df$total_incidence / (5000 * 5)

# Plot in ggplot2 with grey scale lines
ggplot(df, aes(x = w_per, y = incidence_per_capita, color = factor(beta/10000))) +
  geom_line(linewidth = 1) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_color_manual(values = wes_palette("Zissou1", n = 5, type = "continuous")) +
  labs(
    x = expression(w),
    y = "Total Incidence per Capita",
    color = expression(beta),
    title = expression("Per Capita Incidence across " ~ w ~ " and " ~ beta)
  ) +
  theme_minimal(base_size = 14) 
