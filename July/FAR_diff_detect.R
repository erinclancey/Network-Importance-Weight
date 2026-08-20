setwd("~/BRUCE-CHAIN-main/Network MS/July")
source("SETUP Mech Model_July.R")
registerDoParallel()
registerDoRNG(2488420)
p=0.3
gamma=1/10
rho=0.5
k=3



#####################
A_vecs <- matrix(NA, nrow = 100, ncol = 25) 

for(i in 1:nrow(A_vecs)){ 
  A_mat <- matrix(rbinom(25, 1, p), nrow = 5) 
  diag(A_mat) <- 1 
  A_mat[lower.tri(A_mat)] <- t(A_mat)[lower.tri(A_mat)] 
  A_vecs[i, ] <- as.vector(A_mat) 
}
n_sims <- nrow(A_vecs)
w_per <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
beta_par <- 0.00005*10000
day <- seq(1, 200, 1)
Data <- as.data.frame(day)
obs_names <- paste0("reports", 1:5)
max_attempts <- 500

# Initialize base pomp object
Data %>% pomp(
  times = "day", t0 = 0, 
  rprocess = euler(step, delta.t = tau), 
  rinit = rinit, 
  accumvars = accumvars, 
  statenames = statenames, 
  paramnames = paramnames, 
  rmeasure = rmeas, 
  dmeasure = dmeas, 
  rprior = prior_sampler, 
  dprior = prior_density, 
  obsnames = obs_names
) -> gen_data

# --- 2. Build the Expanded Grid ---
# Include sim_id so every parameter pairing maps to its unique matrix connection row
df <- expand.grid(
  sim_id = 1:n_sims,
  w_per = w_per
)

# Pre-allocate tracking columns
df$total_incidence <- rep(NA, nrow(df) )
df$total_cases <- rep(NA, nrow(df) )

# --- 3. Execute Simulation Loop Across Matrices ---
for (i in 1:nrow(df)){
  current_w <- df$w_per[i]
  current_sim_id <- df$sim_id[i]
  
  # Extract the specific matrix vector for this simulation row
  A_vec <- A_vecs[current_sim_id, ]
  
  # Format parameters matching your paramnames definition
  pars <- c(current_w, beta_par, gamma, rho, k, A_vec)
  names(pars) <- paramnames
  
  attempts <- 0
  repeat {
    attempts <- attempts + 1
    
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
    
    if (attempts >= max_attempts) {
      warning(paste("reached max attempts at row", i))
      break
    }
  }
  
  # Extract aggregate simulation counts
  df$total_incidence[i] <- sum(sims[, c("H1", "H2", "H3", "H4", "H5")])
  df$total_cases[i] <- sum(sims[, paste0("reports", 1:5)])
}

# --- 4. Post-processing Calculations ---
df$incidence_per_capita <- df$total_incidence / (5000 * 5)

library(dplyr)

df_transformed <- df %>%
  group_by(sim_id) %>%
  mutate(
    # Subtracts the incidence_per_capita where w_per is 0 from every row in the group
    incidence_diff = incidence_per_capita - incidence_per_capita[w_per == 0]
  ) %>%
  ungroup()

library(ggplot2)

ggplot(df_transformed, aes(x = factor(w_per), y = incidence_diff)) +
  # Adds the box and whisker layers
  geom_boxplot(
    fill = "#0072B2",      # A clean, professional blue fill
    color = "#333333",     # Dark grey outlines
    alpha = 0.5,           # Slight transparency
    outlier.alpha = 0.5,   # Makes outlier dots semi-transparent
    outlier.color = "red"  # Highlights outliers in red
  ) +
  # Adds a horizontal reference line at 0 (the baseline)
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.8) +
  labs(
    x = expression(paste("Importance Weight (", w, ")")),
    y = expression(Delta ~ "Attack Rate" ~ (H[0]-H[A])),
    title = expression(beta~"=5e-05"),
    subtitle = ""
  ) +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black")
  )

# 
# library(dplyr)
# 
# # Run a paired test for each w_per level against the w_per = 0 baseline
# significance_summary <- df_transformed %>%
#   filter(w_per != 0) %>% # Drop 0 since it's the comparison group
#   group_by(w_per) %>%
#   summarize(
#     # Tests if the mean difference is significantly different from 0
#     p_value = t.test(incidence_diff)[["p.value"]],
#     mean_change = mean(incidence_diff)
#   ) %>%
#   mutate(
#     # Apply a Bonferroni correction because you are doing multiple comparisons
#     p_adjusted = p.adjust(p_value, method = "bonferroni"),
#     # Mark significance with standard asterisks
#     significant = ifelse(p_adjusted < 0.05, "Yes", "No")
#   )
# 
# print(significance_summary)





