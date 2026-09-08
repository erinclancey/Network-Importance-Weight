setwd("~/BRUCE-CHAIN-main/Network MS/July")
source("SETUP Mech Model_July.R")
library(ggplot2)
library(dplyr)
library(purrr)
library(HDInterval)
library(grid)
library(ggridges)
registerDoParallel()
registerDoRNG(2488620)

post_0 <- read.csv(file="estim_0_Aug11.csv", header=TRUE)
post_01 <- read.csv(file="estim_01_July13.csv", header=TRUE)
post_02 <- read.csv(file="estim_02_Sep1.csv", header=TRUE)
post_03 <- read.csv(file="estim_03_Sep1.csv", header=TRUE)
post_04 <- read.csv(file="estim_04_Sep1.csv", header=TRUE)
post_05 <- read.csv(file="estim_05_July16.csv", header=TRUE)
post_06 <- read.csv(file="estim_06_Sep8.csv", header=TRUE)
post_07 <- read.csv(file="estim_07_Sep8.csv", header=TRUE)
# post_08 <- read.csv(file="estim_08.csv", header=TRUE)
post_09 <- read.csv(file="estim_09_July16.csv", header=TRUE)
post_1 <- read.csv(file="estim_1_Aug19.csv", header=TRUE)

w_vec <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 1)
post_list <- list(post_0, post_01, post_02, post_03, post_04, post_05,
                  post_06, post_07, post_09, post_1)

AR_post_list <- list()

# --- Run this once to establish your base pomp object ---
day <- seq(1, 200, 1)
Data <- as.data.frame(day)
obs_names <- paste0("reports", 1:5)

p <- 0.3
A_mat <- matrix(rbinom(25, 1, p), nrow=5)
diag(A_mat) <- 1
A_mat[lower.tri(A_mat)] <- t(A_mat)[lower.tri(A_mat)]
A_vec <- as.vector(A_mat)

# Define generic baseline parameters for the pomp structure setup
gamma <- 1/10
rho <- 0.5

# Build the structural pomp object skeleton
Data %>% pomp(
  times="day", t0=0, rprocess=euler(step, delta.t=tau), rinit=rinit,
  accumvars=accumvars, statenames=statenames, paramnames=paramnames,
  rmeasure=rmeas, dmeasure=dmeas, rprior = prior_sampler,
  dprior = prior_density, obsnames=obs_names
) -> gen_data



for (j in 1:length(w_vec)) {
  # Grab the fresh, unmodified posterior from your list
  posterior <- post_list[[j]]
  
  # Sample your 1000 evaluation rows
  HDI_sample <- posterior %>% slice_sample(n = 1000)
  
  # Initialize outcome containers
  HDI_sample$Incidence_alt  <- rep(NA, nrow(HDI_sample))
  HDI_sample$Incidence_null <- rep(NA, nrow(HDI_sample))
  
  max_attempts <- 100 
  
  for (i in 1:nrow(HDI_sample)) {
    w_alt  <- HDI_sample$w[i]
    w_null <- 0
    
    # Check your scale here: if raw posterior beta is already 0.00025, remove '* 10000'
    beta_post <- HDI_sample$beta_par[i]  
    k_post    <- HDI_sample$k[i]
    
    # --- ALTERNATIVE SCENARIO SIMULATION ---
    pars_alt <- c(w_alt, beta_post, gamma, rho, k_post, A_vec)
    names(pars_alt) <- paramnames
    
    attempts <- 0
    sims_alt <- NULL
    repeat {
      attempts <- attempts + 1
      sim <- gen_data %>% simulate(params = pars_alt, nsim = 1, format = "data.frame", include.data = FALSE)
      if (all(colSums(sim[, paste0("reports", 1:5)]) > 0)) {
        sims_alt <- sim
        break
      }
      if (attempts >= max_attempts) { break }
    }
    
    if (!is.null(sims_alt)) {
      HDI_sample$Incidence_alt[i] <- sum(sims_alt[, c("H1", "H2", "H3", "H4", "H5")])
    }
    
    # --- NULL SCENARIO SIMULATION ---
    pars_null <- c(w_null, beta_post, gamma, rho, k_post, A_vec)
    names(pars_null) <- paramnames
    
    attempts <- 0
    sims_null <- NULL
    repeat {
      attempts <- attempts + 1
      sim <- gen_data %>% simulate(params = pars_null, nsim = 1, format = "data.frame", include.data = FALSE)
      if (all(colSums(sim[, paste0("reports", 1:5)]) > 0)) {
        sims_null <- sim
        break
      }
      if (attempts >= max_attempts) { break }
    }
    
    if (!is.null(sims_null)) {
      HDI_sample$Incidence_null[i] <- sum(sims_null[, c("H1", "H2", "H3", "H4", "H5")])
    }
  }
  
  # Calculate final metrics safely using the internal loop dataframe copies
  HDI_sample$attack_rate_alt  <- HDI_sample$Incidence_alt / (5000 * 5)
  HDI_sample$attack_rate_null <- HDI_sample$Incidence_null / (5000 * 5)
  HDI_sample$diff_AR          <- HDI_sample$attack_rate_null - HDI_sample$attack_rate_alt
  
  # Save to the output list
  AR_post_list[[j]] <- HDI_sample
}

########################################

# Combine all posterior samples into one dataframe
plot_dat <- map_dfr(seq_along(AR_post_list), function(i) {
  data.frame(diff_AR = AR_post_list[[i]]$diff_AR,
             k = i,
             w = w_vec[i])
})

# Compute 95% HDI for each w
hdi_dat <- map_dfr(seq_along(AR_post_list), function(i) {
  h <- hdi(AR_post_list[[i]]$diff_AR, ci = 0.95)
  data.frame(k = i,
             w = w_vec[i],
             low = h[1],
             high = h[2])
})

# Faceted plot
ggplot(plot_dat, aes(x = diff_AR)) +
  theme_minimal() +
  geom_histogram(aes(y = after_stat(density)),
                 fill = "#D55E00",
                 color = "#D55E00",
                 alpha = 0.2,
                 bins = 50) +
  geom_density(fill = "#D55E00",
               color = "#D55E00",
               alpha = 0.2,
               adjust = 2) +
  geom_rect(data = hdi_dat,
            aes(xmin = low,
                xmax = high,
                ymin = 0,
                ymax = Inf),
            inherit.aes = FALSE,
            alpha = 0.1,
            fill = "grey20") +
  geom_vline(xintercept = 0,
             color = "red",
             linewidth = 0.75,
             linetype = 2) +
  facet_wrap(~paste0("w = ", w), ncol = 4) +
  labs(x = expression(Delta ~ "Attack Rate" ~ (H[0] - H[A])),
       y = "Density") +
  scale_x_continuous(n.breaks = 6) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)),
                     n.breaks = 6) +
  theme(strip.text = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 15, color = "black"),
        panel.spacing = unit(0.5, "lines"),
        legend.position = "none")


# Combine the list of dataframes into one master plotting dataframe
plot_data <- imap_dfr(AR_post_list, function(df, idx) {
  df %>%
    # Add explicit column for the true w vector matching this loop iteration
    mutate(w_true = w_vec[idx]) %>%
    # Keep only the columns we need to save memory
    select(w_true, diff_AR) %>%
    # Remove rows where simulations failed (NA)
    filter(!is.na(diff_AR))
})


ggplot(plot_data, aes(x = diff_AR, y = as.factor(w_true), fill = w_true, color = w_true)) + 
  geom_density_ridges(alpha = 0.5, scale = 1.5, bandwidth = 0.0005) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1.0) +
  scale_fill_viridis_c(option = "viridis") + 
  scale_color_viridis_c(option = "viridis") + # Ensures the outlines match the fill color scale
  theme_ridges() + 
  labs(x = "Δ Attack Rate (H0 - HA)", y = "Weight (w)") + 
  theme(
    legend.position = "none",
    axis.title.x = element_text(hjust = 0.5),
    axis.title.y = element_text(hjust = 0.5)
  )

################################

posterior_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

sum_df <- plot_data %>%
  group_by(w_true) %>%
  summarise(
    mode = posterior_mode(diff_AR),
    lower = hdi(diff_AR, ci = 0.95)["lower"],
    upper = hdi(diff_AR, ci = 0.95)["upper"],
    .groups = "drop"
  ) 

latex_tab <- sum_df %>%
  mutate(
    Weight = paste0(w_true),
    HDI = paste0(
      "[",
      signif(lower, 2),
      ", ",
      signif(upper, 2),
      "]"
    ),
    `MAP Estimate` = signif(mode, 2)
  ) %>%
  select(
    Weight,
    `MAP Estimate`,
    `95\\% HDI` = HDI
  )

kable(
  latex_tab,
  format = "latex",
  booktabs = TRUE,
  caption = "Parameter estimation results including MAP estimates, 95\\% HDIs and true values across all simulation scenarios.",
  align = c("l","l","l","l","l")
)

hdi_df <- sum_df %>%
  mutate(
    y = seq_along(w_true),
    ymin = y - 0,
    ymax = y + 1
  )

ggplot(plot_data,
       aes(x = diff_AR,
           y = as.factor(w_true),
           fill = w_true,
           color = w_true)) +
  geom_rect(
    data = hdi_df,
    aes(
      xmin = lower,
      xmax = upper,
      ymin = ymin,
      ymax = ymax
    ),
    inherit.aes = FALSE,
    fill = "grey50",
    alpha = 0.4
  ) +
  geom_density_ridges(
    alpha = 0.35,
    scale = 1.5,
    bandwidth = 0.0005
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 1
  ) +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  scale_color_viridis_c(option = "plasma", direction = -1)+
  theme_ridges() +
  labs(
    x = expression(Delta*" Attack Rate ("*H[0]-H[A]*")"),
    y = "Importance Weight (w)"
  ) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 24, hjust = 0.5),
    axis.title.y = element_text(size = 24, hjust = 0.5),
    axis.text.x  = element_text(size = 22),
    axis.text.y  = element_text(size = 22)
  )

