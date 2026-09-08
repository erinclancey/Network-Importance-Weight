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

###########################
chains.long <- map_dfr(seq_along(post_list), \(i)
                       post_list[[i]] %>%
                         select(w, beta_par, k, chain, iter) %>%
                         mutate(
                           beta_par = beta_par / 10000,
                           dataset = paste0("w = ", signif(w_vec[i], 3)),
                           chain = as.character(chain)
                         ) %>%
                         rename(k_par = k) %>%
                         pivot_longer(c(w, beta_par, k_par),
                                      names_to = "variable",
                                      values_to = "value")
) %>%
  mutate(variable = factor(variable, levels = c("w", "beta_par", "k_par")))

plot_names <- as_labeller(
  c(w = "paste(hat(w))",
    beta_par = "paste(hat(beta))",
    k_par = "paste(hat(k))"),
  label_parsed
)

ref_lines <- bind_rows(
  data.frame(variable = "w",
             dataset = paste0("w = ", signif(w_vec, 3)),
             yintercept = w_vec),
  data.frame(variable = "beta_par",
             dataset = unique(chains.long$dataset),
             yintercept = 0.000025),
  data.frame(variable = "k_par",
             dataset = unique(chains.long$dataset),
             yintercept = 3)
) %>%
  mutate(variable = factor(variable, levels = c("w", "beta_par", "k_par")))

ref_lines$dataset <- factor(
  ref_lines$dataset,
  levels = levels(factor(chains.long$dataset))
)

ggplot(chains.long,
       aes(iter, value, group = chain, color = chain)) +
  geom_line(alpha = 0.8) +
  geom_hline(data = ref_lines,
             aes(yintercept = yintercept),
             color = "black",
             linetype = "dashed",
             linewidth = 0.9,
             inherit.aes = FALSE) +
  facet_grid(variable ~ dataset,
             labeller = labeller(variable = plot_names),
             scales = "free_y") +
  scale_color_manual(values = c("#F0E442", "#000000", "#009E73",
                                "#999999", "#0072B2")) +
  labs(x = "Iteration", y = "Parameter Value") +
  theme_minimal(base_size = 20) +
  scale_x_continuous(
    breaks = c(0, 5000, 10000)
  )+
  theme(
      strip.text.x = element_text(size = 24),
      strip.text.y.right = element_text(angle = 0, size = 24, face = "bold"),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 20),
    legend.key.size = unit(1.2, "cm"),
    panel.spacing = unit(0.2, "lines"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    legend.position = "bottom"
  )

