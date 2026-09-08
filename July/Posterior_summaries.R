setwd("~/BRUCE-CHAIN-main/Network MS/July")
source("SETUP Mech Model_July.R")
library(ggplot2)
library(dplyr)
library(purrr)
library(HDInterval)
library(grid)
library(ggridges)
library(patchwork)
library(knitr)
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
M=10000
plot_names <- as_labeller(c(
  'w' = "paste(hat(w))",
  'beta_par' = "paste(hat(beta))",
  'k' = "paste(hat(k))"
), label_parsed)

sum_list <- list()
plot_list <- vector("list", length(w_vec))
beta_plot_list <- vector("list", length(w_vec))
k_plot_list <- vector("list", length(w_vec))

for(j in 1:length(w_vec)){
  
  posterior <- post_list[[j]]
  posterior$beta_par <- posterior$beta_par/10000
  
  post <- posterior %>%
    dplyr::select(w, beta_par, k, chain, iter)
  
  post.list <- split(post, f = post$chain)
  
  mcmc.list <- mcmc.list(list())
  for(i in seq_along(post.list)){
    mcmc.list[[i]] <- mcmc(post.list[[i]])
  }
  
  processed <- window(mcmc.list, start = 2000, end = M + 1, thin = 1)
  processed <- data.frame(do.call(rbind, processed))
  
  processed.long <- processed %>%
    pivot_longer(
      cols = c(w, beta_par, k),
      names_to = "variable",
      values_to = "value"
    )
  
  w_mode <- as.vector(posterior.mode(mcmc(processed$w), adjust = 2))
  w_low  <- ci(processed$w, ci = 0.95, method = "HDI")$CI_low
  w_hi   <- ci(processed$w, ci = 0.95, method = "HDI")$CI_high
  
  beta_par_mode <- as.vector(posterior.mode(mcmc(processed$beta_par), adjust = 2))
  beta_par_low  <- ci(processed$beta_par, ci = 0.95, method = "HDI")$CI_low
  beta_par_hi   <- ci(processed$beta_par, ci = 0.95, method = "HDI")$CI_high
  
  k_mode <- as.vector(posterior.mode(mcmc(processed$k), adjust = 2))
  k_low  <- ci(processed$k, ci = 0.95, method = "HDI")$CI_low
  k_hi   <- ci(processed$k, ci = 0.95, method = "HDI")$CI_high
  
  summaries <- data.frame(
    variable = c("w", "beta_par", "k"),
    true = c(
      w = w_vec[j],
      beta_par = 0.000025,
      k = 3
    ),
    mode = c(w_mode, beta_par_mode, k_mode),
    low = c(w_low, beta_par_low, k_low),
    high = c(w_hi, beta_par_hi, k_hi)
  )
  
  summaries <- summaries %>%
    mutate(across(where(is.numeric), ~signif(.x, digits = 3)))
  
  processed.long_w <- processed.long %>%
    filter(variable == "w")
  
  processed.long_beta <- processed.long %>%
    filter(variable == "beta_par")
  
  processed.long_k <- processed.long %>%
    filter(variable == "k")
  
  summaries_w <- summaries %>%
    filter(variable == "w")
  
  summaries_beta <- summaries %>%
    filter(variable == "beta_par")
  
  summaries_k <- summaries %>%
    filter(variable == "k")
  
  blank_data <- data.frame(
    variable = factor("w", levels = "w"),
    value = c(0, 1)
  )
  
  plot_list[[j]] <-
    ggplot(processed.long_w,
           aes(x = value,
               fill = variable,
               color = variable)) +
    theme_minimal() +
    geom_blank(data = blank_data) +
    geom_histogram(aes(y = after_stat(density)),
                   position = "identity",
                   alpha = 0.2,
                   bins = 50) +
    geom_density(alpha = .2, adjust = 2) +
    geom_vline(data = summaries_w,
               aes(xintercept = mode),
               color = "black",
               linewidth = 0.75,
               linetype = 2) +
    geom_vline(data = summaries_w,
               aes(xintercept = true),
               color = "red",
               linewidth = 0.75,
               linetype = 2) +
    geom_rect(data = summaries_w,
              aes(xmin = low,
                  xmax = high,
                  ymin = 0,
                  ymax = Inf),
              inherit.aes = FALSE,
              alpha = 0.1,
              fill = "grey20") +
    scale_fill_manual(values = c("w" = "#0072B2")) +
    scale_color_manual(values = c("w" = "#0072B2")) +
    labs(
      title = paste0("w = ", w_vec[j]),
      x = expression(hat(w)),
      y = "Density"
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      n.breaks = 6
    ) +
    scale_y_continuous(
      expand = expansion(mult = 0.2),
      n.breaks = 6
    ) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text = element_text(size = 12, color = "black"),
      legend.position = "none"
    )
  
  beta_plot_list[[j]] <-
    ggplot(processed.long_beta,
           aes(x = value,
               fill = variable,
               color = variable)) +
    theme_minimal() +
    geom_histogram(aes(y = after_stat(density)),
                   position = "identity",
                   alpha = 0.2,
                   bins = 50) +
    geom_density(alpha = .2, adjust = 2) +
    geom_vline(data = summaries_beta,
               aes(xintercept = mode),
               color = "black",
               linewidth = 0.75,
               linetype = 2) +
    geom_vline(xintercept = 0.000025,
               color = "red",
               linewidth = 0.75,
               linetype = 2) +
    geom_rect(data = summaries_beta,
              aes(xmin = low,
                  xmax = high,
                  ymin = 0,
                  ymax = Inf),
              inherit.aes = FALSE,
              alpha = 0.1,
              fill = "grey20") +
    scale_fill_manual(values = c("beta_par" = "grey40")) +
    scale_color_manual(values = c("beta_par" = "grey40")) +
    labs(
      title = paste0("w = ", w_vec[j]),
      x = expression(hat(beta)),
      y = "Density"
    ) +
    scale_x_continuous(n.breaks = 6) +
    scale_y_continuous(
      expand = expansion(mult = 0.2),
      n.breaks = 6
    ) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text = element_text(size = 12, color = "black"),
      legend.position = "none"
    )
  
  k_plot_list[[j]] <-
    ggplot(processed.long_k,
           aes(x = value,
               fill = variable,
               color = variable)) +
    theme_minimal() +
    geom_histogram(aes(y = after_stat(density)),
                   position = "identity",
                   alpha = 0.2,
                   bins = 50) +
    geom_density(alpha = .2, adjust = 2) +
    geom_vline(data = summaries_k,
               aes(xintercept = mode),
               color = "black",
               linewidth = 0.75,
               linetype = 2) +
    geom_vline(xintercept = 3,
               color = "red",
               linewidth = 0.75,
               linetype = 2) +
    geom_rect(data = summaries_k,
              aes(xmin = low,
                  xmax = high,
                  ymin = 0,
                  ymax = Inf),
              inherit.aes = FALSE,
              alpha = 0.1,
              fill = "grey20") +
    scale_fill_manual(values = c("k" = "grey40")) +
    scale_color_manual(values = c("k" = "grey40")) +
    labs(
      title = paste0("w = ", w_vec[j]),
      x = expression(hat(k)),
      y = "Density"
    ) +
    scale_x_continuous(n.breaks = 6) +
    scale_y_continuous(
      expand = expansion(mult = 0.2),
      n.breaks = 6
    ) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text = element_text(size = 12, color = "black"),
      legend.position = "none"
    )
  
  clean_summary_df <- function(summary_obj){
    df <- as.data.frame(summary_obj)
    
    if(!"Parameter" %in% colnames(df)){
      df <- cbind(Parameter = rownames(df), df)
    }
    
    rownames(df) <- NULL
    
    cols_to_convert <- c("true", "mode", "low", "high")
    
    for(col in cols_to_convert){
      if(col %in% colnames(df)){
        numeric_values <- as.numeric(as.character(df[[col]]))
        df[[col]] <- format(
          round(numeric_values, 10),
          scientific = FALSE,
          drop0trailing = TRUE
        )
      }
    }
    
    df
  }
  
  cleaned_summary <- clean_summary_df(summaries)[,-1]
  sum_list[[j]] <- cleaned_summary
}


wrap_plots(plot_list, ncol = 3) +
  plot_annotation(
    title = expression("Posterior Distributions of " * hat(w)),
    theme = theme(
      plot.title = element_text(
        size = 28,
        face = "bold",
        hjust = 0.5
      )
    )
  )

wrap_plots(beta_plot_list, ncol = 3) +
  plot_annotation(
    title = expression("Posterior Distributions of " * hat(beta)),
    theme = theme(
      plot.title = element_text(
        size = 28,
        face = "bold",
        hjust = 0.5
      )
    )
  )

wrap_plots(k_plot_list, ncol = 3) +
  plot_annotation(
    title = expression("Posterior Distributions of " * hat(k)),
    theme = theme(
      plot.title = element_text(
        size = 28,
        face = "bold",
        hjust = 0.5
      )
    )
  )

latex_tab <- bind_rows(
  lapply(seq_along(sum_list), function(i) {
    sum_list[[i]] %>%
      mutate(
        Weight = paste0(w_vec[i]),
        HDI = paste0("[", low, ", ", high, "]")
      ) %>%
      select(`Weight`, Parameter = variable,
             `MAP Estimate` = mode, `95\\% HDI` = HDI)
  })
)

kable(
  latex_tab,
  format = "latex",
  booktabs = TRUE,
  caption = "Parameter estimation results including MAP estimates, 95\\% HDIs and true values across all simulation scenarios.",
  align = c("l","l","l","l","l")
)

