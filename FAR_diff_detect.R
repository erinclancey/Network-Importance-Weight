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
w_per <- c(0, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
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
View(df)



# Extract the incidence_per_capita of the first 100 rows
baseline_vals <- df$incidence_per_capita[1:100]
baseline_ince <- df$total_incidence[1:100]

df_transformed <- df[-(1:100), ] %>%
  mutate(
    # Use modulo or row indexing to map rows back to the 1-100 baseline vector
    # Assuming your grid cycles through the 100 sims or you want to match by row index
    incidence_diff = incidence_per_capita - baseline_vals[(row_number() - 1) %% 100 + 1]
  )
View(df_transformed)

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
    title = expression(beta~"=1.5e-05"),
    subtitle = ""
  ) +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black")
  )

expected_diff <- mean(df_transformed[df_transformed$w_per == 0.2, ]$incidence_diff)
t.test(df_transformed[df_transformed$w_per == 0.2, ]$incidence_diff, mu = 0)


library(dplyr)
library(purrr)

# Example paired vectors
baseline <- baseline_ince
treatment <- df_transformed[df_transformed$w_per == 0.2, ]$total_incidence

# Run paired t-test
paired_result <- t.test(treatment, baseline, paired = TRUE)

print(paired_result)

# Example data: 
# Group 1: 45 successes out of 500 trials
# Group 2: 70 successes out of 500 trials

successes <- c(25000, 25000)
trials    <- c(25000, 25000)

# Run the two-proportion test without continuity correction
prop_test_result <- prop.test(x = successes, n = trials, correct = FALSE)

print(prop_test_result)

library(dplyr)
library(broom) # For clean tidying of test outputs

# 1. Run the hypothesis tests for each w_per group
test_results <- df_transformed %>%
  group_by(w_per) %>%
  do(tidy(t.test(.$incidence_diff, mu = 0))) %>%
  ungroup() %>%
  mutate(
    # Apply Benjamini-Hochberg (FDR) correction for multiple testing
    p_adj = p.value * n(), # or use p.adjust(p.value, method = "BH")
    p_adj_BH = p.adjust(p.value, method = "BH"),
    # Create significance stars for plotting
    sig_label = case_when(
      p_adj_BH < 0.001 ~ "***",
      p_adj_BH < 0.01  ~ "**",
      p_adj_BH < 0.05  ~ "*",
      TRUE             ~ "ns"
    )
  )

# View the full statistical table
print(test_results %>% select(w_per, estimate, statistic, p.value, p_adj_BH, sig_label))

# 2. (Optional) Find the maximum y-value per group to position significance labels on your plot
max_y_vals <- df_transformed %>%
  group_by(w_per) %>%
  summarize(y_max = max(incidence_diff, na.rm = TRUE), .groups = "drop")

plot_annotations <- left_join(test_results, max_y_vals, by = "w_per")

