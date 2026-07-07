setwd("~/BRUCE-CHAIN-main/Network MS/July 2026")
source("SETUP Mech Model_July.R")

iter_sim_study <- data.frame(matrix(ncol = 6, nrow = nrow(pars)))
x <- c("w_true","w_hat",
       "beta_par_true", "beta_par_hat",
       "k_true", "k_hat")
colnames(iter_sim_study) <- x
#############################################################
for(i in 1:nrow(pars)){
# 2. Build your pomp object
sims <- data_list[[i]] %>% select(day, paste0("reports", 1:5))

sims %>%
  select(day, all_of(obs_names)) %>%
  pomp(
    times      = "day",
    t0         = 0,
    rprocess   = euler(step, delta.t = tau),
    rinit      = rinit,
    statenames = statenames,
    accumvars  = accumvars,
    paramnames = paramnames,
    rmeasure   = rmeas,
    dmeasure   = dmeas,
    obsnames   = obs_names
  ) -> sim_dat


#set true parameter values
iter_sim_study$w_true[i] <- pars$w[i]
iter_sim_study$beta_par_true[i] <- pars$beta_par[i]
iter_sim_study$k_true[i] <- pars$k[i]

# Set Up Chains
nchains=1
w = runif(nchains, 0.4,0.6)
beta_par = runif(nchains, 0.00003*10000, 0.00004*10000)
gamma = rep(pars["gamma"][i,], nchains)
rho = rep(pars["rho"][i,], nchains)
k = runif(nchains, 3,5)
A_frame =  pars[i, 6:30]
rownames(A_frame) <- NULL
theta.start <- data.frame(cbind(w,beta_par,gamma,rho,k, A_frame))

  
mif2(
  data = sim_dat,           # Explicitly pass your pomp object here
  Nmif = 500,                  
  Np = 500,                  
  cooling.fraction.50 = 0.5,  
  rw.sd = rw_sd(
    w = 0.02,            
    beta_par = 0.02,
    k = 0.1
  ),
  params = theta.start[1, ]   # Your starting parameter vector
) -> mif_output

# Extract final MLE parameter estimates
mle_estimates <- coef(mif_output)

iter_sim_study$w_hat[i] <- mle_estimates[1]
iter_sim_study$beta_par_hat[i] <- mle_estimates[2]
iter_sim_study$k_hat[i] <- mle_estimates[5]
print(i)
}

write.csv(iter_sim_study, file = "iter_sim_study.csv", row.names = FALSE)

#####SET 1
# Make linear prediction for plots
wmod <- summary(lm(w_hat ~ w_true , data = iter_sim_study))
betamod <- summary(lm(beta_par_hat ~ beta_par_true , data = iter_sim_study))
kmod <- summary(lm(k_hat ~ k_true , data = iter_sim_study))




library(ggplot2)
library(latex2exp)
library(patchwork)
library(scales)

## ---- Plot 1: w recovery ----
p1 <- ggplot(iter_sim_study, aes(x = w_true, y = w_hat)) +
  geom_point(color = "darkorange2", shape = 20, size = 3, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, 
              linetype = 2, color = "black", size = 1) +
  geom_abline(intercept = wmod$coefficients[1,1], 
              slope = wmod$coefficients[2,1],
              linetype = 1, color = "darkorange4", size = 1) +
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

p2 <- ggplot(iter_sim_study, aes(x = beta_par_true / 10000,
                         y = beta_par_hat / 10000)) +
  geom_point(color = "steelblue", shape = 20, size = 3, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1,
              linetype = 2, color = "black", size = 1) +
  geom_abline(intercept = betamod$coefficients[1,1] / 10000,
              slope = betamod$coefficients[2,1],
              linetype = 1, color = "steelblue4", size = 1) +
  scale_x_continuous(
    limits = c(0, 0.8/10000),
    labels = sci_hybrid
  ) +
  scale_y_continuous(
    limits = c(0, 0.8/10000),
    labels = sci_hybrid
  ) +
  labs(x = TeX("$\\beta$"),
       y = TeX("$\\hat{\\beta}$")) +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 16))

## ---- Plot 1: w recovery ----
p3 <- ggplot(iter_sim_study, aes(x = k_true, y = k_hat)) +
  geom_point(color = "darkolivegreen4", shape = 20, size = 3, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, 
              linetype = 2, color = "black", size = 1) +
  geom_abline(intercept = kmod$coefficients[1,1], 
              slope = kmod$coefficients[2,1],
              linetype = 1, color = "darkolivegreen", size = 1) +
  xlim(0, 10) + ylim(0, 10) +
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












