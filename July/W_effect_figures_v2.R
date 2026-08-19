setwd("~/BRUCE-CHAIN-main/Network MS/July")
source("SETUP Mech Model_July.R")
registerDoParallel()
registerDoRNG(2488220)
p=0.3
A_mat <- matrix(rbinom(25,1,p), nrow=5)
diag(A_mat) <- 1
A_mat[lower.tri(A_mat)] <- t(A_mat)[lower.tri(A_mat)]
#A_mat == t(A_mat)   # should be all TRUE
A_vec <- as.vector(A_mat)
A_mat 
gamma=1/10
rho=0.5
k=3



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

ggplot(df, aes(x = w_per, y = R0_pop, color = factor(as.numeric(factor(beta))))) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_y_continuous(
    limits = c(0, NA),
    breaks = seq(0, max(df$R0_pop, na.rm = TRUE), by = 2)
  ) +
  scale_color_manual(
    values = c(
      "1" = "#999999",  # Grey
      "2" = "#CC79A7", # Pink
      "3" = "#009E73", # Green
      "4" = "#E69F00", # Orange
      "5" = "#0072B2" # Blue
    ),
    labels = unique(df$beta) # Restores original beta values in the legend
  ) +
  labs(
    x = expression(w),
    y = expression("Total Population" ~ R[0]),
    color = expression(beta),
    title = "A"
  ) +
  theme_minimal(base_size = 20)



day <- seq(1, 500, 1)
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


w_per <- seq(0, 1, length.out = 500)
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

ggplot(df, aes(x = w_per, y = incidence_per_capita, color = factor(as.numeric(factor(beta))))) +
  geom_line(linewidth = 1) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_color_manual(
    values = c(
      "1" = "#999999",  # Grey
      "2" = "#CC79A7", # Pink
      "3" = "#009E73", # Green
      "4" = "#E69F00", # Orange
      "5" = "#0072B2"  # Blue
    ),
    labels = sort(unique(df$beta)) / 10000  # Divides beta by 10000 for legend labels
  ) +
  labs(
    x = expression(w),
    y = "Attack Rate",
    color = expression(beta),
    title = "B"
  ) +
  theme_minimal(base_size = 20)
