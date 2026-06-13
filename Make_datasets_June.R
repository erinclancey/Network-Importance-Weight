setwd("~/BRUCE-CHAIN-main/Network MS")
source("Network Mech Model Pomp_SETUP_prior.R")
registerDoParallel()
registerDoRNG(2488420)


############# Make an empty pomp object ############# 
day <- seq(1, 200, 1)
empty_df <- as.data.frame(day)
obs_names <- paste0("reports", 1:5)
empty_data <- pomp(data=empty_df,
                   times="day",t0=0,
                   rprocess=euler(step,delta.t=tau),
                   rinit=rinit,
                   accumvars=accumvars,
                   statenames=statenames,
                   paramnames=paramnames,
                   rmeasure=rmeas,
                   dmeasure=dmeas,
                   obsnames=obs_names)
###############set up true  parameter grid for simulating data########
n=200
w = runif(n,0,1)
beta_par = runif(n, 0.000015*10000, 0.00005*10000)
gamma = runif(n, 1/10, 1/10)
rho = runif(n, 0.25,0.75)
k = runif(n, 1,5)
p = runif(n, 0.2,0.5)

params <- cbind(w, beta_par, gamma, rho, k)
colnames(params) <- paramnames[c(1:5)]
A_NA =  matrix(NA, nrow = n, ncol = 25)
colnames(A_NA) <- A_names
pars_NA <- matrix(NA, nrow = n, ncol = 5)
colnames(pars_NA) <- paramnames[c(1:5)]
pars <- cbind(pars_NA,A_NA)
R0_local <- vector()
R0_total <- vector()
data_list <- list()
end_day <- vector()
total_incidence <- vector()
total_cases <- vector()
##########################Make parameter dataframe
for (i in 1:nrow(pars)){
  repeat {
    # Matrix A generated then hard coded into step
    A_mat <- matrix(rbinom(25,1,p[i]), nrow=5)
    diag(A_mat) <- 1
    A_mat[lower.tri(A_mat)] <- t(A_mat)[lower.tri(A_mat)]
    A_vec <- as.vector(A_mat)
    names(A_vec) <- A_names
    A_row <- t(as.matrix(A_vec))
    pars[i,] <- c(params[i,],A_row)
    
    empty_data %>%
      simulate(
        params=pars[i,],
        nsim=1,
        format="data.frame",include.data=FALSE) -> sim_dat
    
    # Check if all reports1–reports5 have at least one non-zero value
    nonzero_check <- colSums(sim_dat[ , paste0("reports", 1:5)]) > 0
    
    if (all(nonzero_check)) {
      break
    }
  }
  end_day_row <- sim_dat %>%
    filter(day > 50) %>%                 # only rows after day 50
    filter(I1 == 0, I2 == 0, I3 == 0, I4 == 0, I5 == 0) %>%  # all I’s zero
    slice(1)
  val <- end_day_row %>% dplyr::pull(day) %>% .[1]
  end_day[i] <- ifelse(length(val) == 0 || is.na(val), 200, val)
  
  
  total_incidence[i] <- sum(sim_dat[, c("H1","H2","H3","H4","H5")])
  total_cases[i] <- sum(sim_dat[, c("reports1","reports2","reports3","reports4","reports5")])

  data_list[[i]] <- sim_dat
  # Build Cmat
  #THis is working now because all populations are the same size. 
  R0_local[i] <- 5000*pars[i,"beta_par"]/10000 / pars[i,"gamma"]
  
  #######Calculate R0
  # Build Cmat
  Cmat <- A_mat * w[i] + (1 - w[i])   # off-diagonal formula applied everywhere
  diag(Cmat) <- 1 
  # Compute eigen decomposition
  eig <- eigen(Cmat)
  # Eigenvalues
  dom_eig <- max(eig$values)
  R0_total[i] <- 5000*pars[i,"beta_par"]/10000 / pars[i,"gamma"] * dom_eig
  
  
}

pars <- data.frame(pars) 
sim_dat_metadata <- as.data.frame(cbind(pars, p, R0_local, R0_total, end_day, total_incidence, total_cases))

write.csv(sim_dat_metadata, file="sim_dat_metadata_June13.csv")


  