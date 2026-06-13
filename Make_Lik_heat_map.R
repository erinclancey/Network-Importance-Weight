setwd("~/BRUCE-CHAIN-main/Network MS")
source("Network Mech Model Pomp_SETUP (1).R")

# 1. Define the parameter grid arrays
w_vals <- seq(0, 1, by = 0.01)   # Adjust range/grain as needed
p_vals <- seq(0.2, 1, by = 0.01)   # Adjust range/grain as needed

# Create all combinations
grid_results <- expand.grid(w = w_vals, p = p_vals)
grid_results$LL_diff <- NA  # Placeholder for results

# Fixed Simulation Framework
day <- seq(1, 200, 1)
Data <- as.data.frame(day)
obs_names <- paste0("reports", 1:5)
max_attempts <- 50  # Lowered safety break slightly for speed over large loops

# 2. Iterate through the parameter grid
for (i in 1:nrow(grid_results)) {
  
  current_w <- grid_results$w[i]
  current_p <- grid_results$p[i]
  
  attempt <- 0
  success <- FALSE
  
  # Stochastic simulation loop based on current p
  while (attempt < max_attempts) {
    attempt <- attempt + 1
    
    # Generate network matrix using current p
    A_mat <- matrix(rbinom(25, 1, current_p), nrow = 5)
    diag(A_mat) <- 1
    A_mat[lower.tri(A_mat)] <- t(A_mat)[lower.tri(A_mat)]
    A_vec <- as.vector(A_mat)
    
    # Baseline Parameters
    beta_par <- 0.000025 * 10000
    gamma <- 1 / 10
    rho <- 0.5
    k <- 3
    
    pars <- c(current_w, beta_par, gamma, rho, k, A_vec)
    names(pars) <- paramnames
    
    # Build simulation container
    gen_data <- Data %>% 
      pomp(
        times = "day", t0 = 0,
        rprocess = euler(step, delta.t = tau),
        rinit = rinit, accumvars = accumvars,
        statenames = statenames, paramnames = paramnames,
        rmeasure = rmeas, dmeasure = dmeas, obsnames = obs_names
      )
    
    # Simulate
    sims <- gen_data %>% 
      simulate(params = pars, nsim = 1, format = "data.frame", include.data = FALSE)
    
    # Verify non-zero signals across reports
    nonzero_check <- colSums(sims[, obs_names]) > 0
    if (all(nonzero_check)) {
      success <- TRUE
      break
    }
  }
  
  # Skip grid coordinate if data simulation persistently failed to seed
  if (!success) {
    next
  }
  
  # 3. Fit data back into pomp object to run particle filters
  sim_dat <- sims %>%
    select(day, all_of(obs_names)) %>%
    pomp(
      times = "day", t0 = 0,
      rprocess = euler(step, delta.t = tau),
      rinit = rinit, statenames = statenames,
      accumvars = accumvars, paramnames = paramnames,
      rmeasure = rmeas, dmeasure = dmeas, obsnames = obs_names
    )
  
  make_prof <- sim_dat %>% pomp(params = pars)
  
  # True vs Null Evaluation
  theta <- pars
  theta_null <- replace(theta, 1, 0)
  
  # Run particle filters (consider Np = 1000 or higher if surface is too noisy)
  LL_true <- logLik(make_prof %>% pfilter(params = theta, Np = 500))
  LL_null <- logLik(make_prof %>% pfilter(params = theta_null, Np = 500))
  
  # Save the test statistic into our grid tracker
  grid_results$LL_diff[i] <- LL_true - LL_null
}

grid_results$LL_diff[grid_results$LL_diff < 0] <- 0
grid_results$LL_diff_log <- log(grid_results$LL_diff, base=10)
View(grid_results)
# Generate Heat Map Plot
heatmap_plot <- ggplot(grid_results, aes(x = w, y = p, fill = LL_diff_log)) +
  geom_tile(color = "white", lwd = 0.2, linetype = 1) +
  scale_fill_viridis_c(name = "LL Diff", option = "viridis", na.value = "grey90") +
  labs(
    title = "Likelihood Ratio Profile Surface (LL_diff)",
    subtitle = "Varying Entry Weight (w) and Connectivity Probability (p)",
    x = "Network Entry Weight (w)",
    y = "Matrix Connection Probability (p)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold")
  )

# Render Plot
print(heatmap_plot)

