setwd("~/BRUCE-CHAIN-main/Network MS/Sim_Study_FINAL_Sep2026")
library(tidyverse)
library(ggplot2)
library(latex2exp)
library(patchwork)
library(scales)
library(dplyr)

sim_df_original <- read.csv(file="simstudy_original_50.csv", header=TRUE)
sim_new_rows <- read.csv(file="simstudy_new_22rows.csv", header=TRUE)
sim_new_rows <- sim_new_rows %>% select(-logLik)

sim_df_updated <- sim_df_original %>%
  rows_update(sim_new_rows, by = "X")

sim_df <- sim_df_updated
sim_df <- na.omit(sim_df)
nrow(sim_df)


# Percent convergence
sim_df <- subset(sim_df, sim_df$Gel_w<1.05 & sim_df$Gel_beta <1.05 & sim_df$Gel_k<1.05)
nrow(sim_df)

# Coverage
sim_df <- sim_df %>%
  mutate(
    w_inHDI = ifelse(w_true >= w_low_H & w_true <= w_hi_H, 1, 0),
    beta_inHDI = ifelse(beta_par_true >= beta_par_low_H & beta_par_true <= beta_par_hi_H, 1, 0),
    k_inHDI = ifelse(k_true >= k_low_H & k_true <= k_hi_H, 1, 0)
  )

mean(sim_df$w_true >= sim_df$w_low_H & sim_df$w_true <= sim_df$w_hi_H)
mean(sim_df$beta_par_true >= sim_df$beta_par_low_H & sim_df$beta_par_true <= sim_df$beta_par_hi_H)
mean(sim_df$k_true >= sim_df$k_low_H & sim_df$k_true <= sim_df$k_hi_H)
sum(sim_df$w_true >= sim_df$w_low_H & sim_df$w_true <= sim_df$w_hi_H)
sum(sim_df$beta_par_true >= sim_df$beta_par_low_H & sim_df$beta_par_true <= sim_df$beta_par_hi_H)
sum(sim_df$k_true >= sim_df$k_low_H & sim_df$k_true <= sim_df$k_hi_H)

# Relative Bias, RMSE and COR 

sim_df %>%
  summarize(
    w_rel_bias = 100 * mean((w_mode - w_true) / w_true),
    w_rmse     = sqrt(mean((w_mode - w_true)^2)),
    w_cor      = cor(w_true, w_mode),
    
    beta_rel_bias = 100 * mean((beta_par_mode - beta_par_true) / beta_par_true),
    beta_rmse     = sqrt(mean((beta_par_mode - beta_par_true)^2)),
    beta_cor      = cor(beta_par_true, beta_par_mode),
    
    k_rel_bias = 100 * mean((k_mode - k_true) / k_true),
    k_rmse     = sqrt(mean((k_mode - k_true)^2)),
    k_cor      = cor(k_true, k_mode)
  )

#####SET 1
# Make linear prediction for plots
wmod <- summary(lm(w_mode ~ w_true , data = sim_df))
betamod <- summary(lm(beta_par_mode ~ beta_par_true , data = sim_df))
kmod <- summary(lm(k_mode ~ k_true , data = sim_df))
nrow(sim_df)



## ---- Plot 1: w recovery ----
p1 <- ggplot(sim_df, aes(x = w_true, y = w_mode)) +
  geom_point(color = "#0072B2", shape = 20, size = 3, alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, 
              linetype = 2, color = "black", size = 1) +
  geom_abline(intercept = wmod$coefficients[1,1], 
              slope = wmod$coefficients[2,1],
              linetype = 1, color = "#0072B2", size = 1) +
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
  geom_point(color = "#E69F00", shape = 20, size = 3, alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1,
              linetype = 2, color = "black", size = 1) +
  geom_abline(intercept = betamod$coefficients[1,1] / 10000,
              slope = betamod$coefficients[2,1],
              linetype = 1, color = "#E69F00", size = 1) +
  scale_x_continuous(
    limits = c(0.2/10000, 0.6/10000),
    labels = sci_hybrid
  ) +
  scale_y_continuous(
    limits = c(0.2/10000, 0.6/10000),
    labels = sci_hybrid
  ) +
  labs(x = TeX("$\\beta$"),
       y = TeX("$\\hat{\\beta}$")) +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 16))

## ---- Plot 1: w recovery ----
p3 <- ggplot(sim_df, aes(x = k_true, y = k_mode)) +
  geom_point(color =  "#CC79A7", shape = 20, size = 3, alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, 
              linetype = 2, color = "black", size = 1) +
  geom_abline(intercept = kmod$coefficients[1,1], 
              slope = kmod$coefficients[2,1],
              linetype = 1, color =  "#CC79A7", size = 1) +
  xlim(0, 8) + ylim(0, 8) +
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











