setwd("~/BRUCE-CHAIN-main/Network MS/July")
source("SETUP Mech Model_July.R")
registerDoParallel()
registerDoRNG(2488220)


gamma=1/10
rho=0.5
k=3
p=0.3

A_vecs <- matrix(NA, nrow = 100, ncol = 25) 

for(i in 1:nrow(A_vecs)){ 
  A_mat <- matrix(rbinom(25, 1, p), nrow = 5) 
  diag(A_mat) <- 1 
  A_mat[lower.tri(A_mat)] <- t(A_mat)[lower.tri(A_mat)] 
  A_vecs[i, ] <- as.vector(A_mat) 
}

# --- Setup Parameters ---
w_per <- seq(0, 1, length.out = 100) # Reduced to 100 for speed; change back to 1000 if needed
beta_par <- c(0.000005,0.00001, 0.000015, 0.000025, 0.00005)
n_sims <- nrow(A_vecs) # Number of rows in A_vecs

# --- 1. Compute dominant eigenvalues for ALL 1,000 matrices ---
# We create a long dataframe storing the eigenvalues for every simulation and w_per
eig_list <- list()

for (s in 1:n_sims) {
  # Reconstruct the 5x5 A_mat from the s-th row of A_vecs
  A_mat <- matrix(A_vecs[s, ], nrow = 5, ncol = 5)
  
  dom_eig <- numeric(length(w_per))
  for (i in seq_along(w_per)) {
    Cmat <- A_mat * w_per[i] + (1 - w_per[i])
    diag(Cmat) <- 1
    dom_eig[i] <- max(Re(eigen(Cmat)$values)) # real() handles tiny numerical imaginary parts
  }
  
  # Store results for this matrix simulation
  eig_list[[s]] <- data.frame(
    sim_id = s,
    w_per = w_per,
    dom_eig = dom_eig
  )
}

# Combine all simulations into one master eigenvalue dataframe
all_eigs <- do.call(rbind, eig_list)

# --- 2. Build the full grid incorporating beta ---
# Create a grid combining all simulation/w_per data with all beta values
df_beta <- data.frame(beta = beta_par)
df <- merge(all_eigs, df_beta, all = TRUE)

# --- 3. Calculate R0_pop ---
df$R0_pop <- (5000 * df$beta / gamma) * df$dom_eig

# Create an explicit group identifier combining sim_id and beta
df$group_id <- paste0(df$sim_id, "_", df$beta)

# --- 4. Plot the Grid ---
ggplot(df, aes(x = w_per, y = R0_pop, 
               color = factor(as.numeric(factor(beta))), 
               group = group_id)) + # Grouping ensures 1000 distinct lines per beta
  geom_line(linewidth = 0.5, alpha = 0.15) + # Thin lines and low alpha to handle 1000 lines
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
  scale_y_continuous( 
    limits = c(0, NA), 
    breaks = seq(0, max(df$R0_pop, na.rm = TRUE) + 1, by = 2) 
  ) + 
  scale_color_manual( 
    values = c( 
      "1" = "#999999", # Grey 
      "2" = "#CC79A7", # Pink 
      "3" = "#009E73", # Green 
      "4" = "#E69F00", # Orange 
      "5" = "#0072B2"  # Blue 
    ), 
    labels = unique(df$beta) 
  ) + 
  labs( 
    x = "Importance Weight (w)", 
    y = expression("Total Population" ~ R[0]), 
    color = expression(beta), 
    title = "A" 
  ) + 
  theme_minimal(base_size = 20) + 
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",       # Moves the legend below the plot
    legend.direction = "horizontal"    # Aligns the legend items side-by-side
  ) + 
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) # Fixes legend fading







#####################
A_vecs <- matrix(NA, nrow = 100, ncol = 25) 

for(i in 1:nrow(A_vecs)){ 
  A_mat <- matrix(rbinom(25, 1, p), nrow = 5) 
  diag(A_mat) <- 1 
  A_mat[lower.tri(A_mat)] <- t(A_mat)[lower.tri(A_mat)] 
  A_vecs[i, ] <- as.vector(A_mat) 
}
n_sims <- nrow(A_vecs)
w_per <- seq(0, 1, length.out = 100)
beta_par <- c(0.000005, 0.000015, 0.000025, 0.00005, 0.00001)*10000
day <- seq(1, 200, 1)
Data <- as.data.frame(day)
obs_names <- paste0("reports", 1:5)
max_attempts <- 100

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
  w_per = w_per,
  beta = beta_par
)

# Pre-allocate tracking columns
df$total_incidence <- rep(NA, nrow(df) )
df$total_cases <- rep(NA, nrow(df) )

# --- 3. Execute Simulation Loop Across Matrices ---
for (i in 1:nrow(df)){
  current_w <- df$w_per[i]
  current_beta <- df$beta[i]
  current_sim_id <- df$sim_id[i]
  
  # Extract the specific matrix vector for this simulation row
  A_vec <- A_vecs[current_sim_id, ]
  
  # Format parameters matching your paramnames definition
  pars <- c(current_w, current_beta, gamma, rho, k, A_vec)
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

# Create a unique group string for every single matrix-line series
df$group_id <- paste0(df$sim_id, "_", df$beta)

# --- 5. Generate Transparency Grid Plot ---
ggplot(df, aes(x = w_per, y = incidence_per_capita, 
               color = factor(as.numeric(factor(beta))), 
               group = group_id)) + # Grouping ensures 1000 lines draw separately
  geom_line(linewidth = 0.5, alpha = 0.15) + # Transparent lines to overlay density distributions
  scale_y_continuous(
    limits = c(0, NA), 
    expand = expansion(mult = c(0, 0.05))
  ) + 
  scale_color_manual(
    values = c(
      "1" = "#999999", # Grey
      "2" = "#CC79A7", # Pink
      "3" = "#009E73", # Green
      "4" = "#E69F00", # Orange
      "5" = "#0072B2"  # Blue
    ),
    labels = sort(unique(df$beta)) / 10000 
  ) + 
  labs(
    x = "Importance Weight (w)",
    y = "Attack Rate",
    color = expression(beta),
    title = "B"
  ) + 
  theme_minimal(base_size = 20) + 
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",       # Moves the legend below the plot
    legend.direction = "horizontal"    # Aligns the legend items side-by-side
  ) + 
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) # Solidifies side legend








