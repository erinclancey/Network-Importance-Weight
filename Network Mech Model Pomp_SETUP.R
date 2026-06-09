## ----General Setup----------------------------------------------------------
library(foreach)
library(iterators)
library(parallel)
library(rngtools)
library(doParallel)
library(doRNG)
library(pomp)
library(tidyverse)
library(dplyr)
library(tidyr)
library(coda)
library(MCMCglmm)
library(bayestestR)
library(patchwork)
# runs the ci function to get credible intervals
setwd("~/BRUCE-CHAIN-main/Network MS")
#dir.create("tmp")
options(pomp_cdir="./tmp")



## ----Model Setup----------------------------------------------------------

A_names <- paste0("A", seq(1,25,1))

accumvars=c("H1","H2","H3","H4","H5")
statenames=c("S1","S2","S3","S4","S5",
             "I1","I2","I3","I4","I5",
             "R1","R2","R3","R4","R5",
             "H1","H2","H3","H4","H5")
paramnames=c("w","beta_par","gamma", "rho","k", A_names)

# time step
tau=1

## ----Process Model----------------------------------------------------------
step <- Csnippet("
  double S[5] = {S1,S2,S3,S4,S5};
  double I[5] = {I1,I2,I3,I4,I5};
  double R[5] = {R1,R2,R3,R4,R5};
  double H[5] = {H1,H2,H3,H4,H5};

  int n = 5;
  int i, j;

  // adjacency matrix built from parameters A1..A25
  int A[5][5] = {
    {A1, A2, A3, A4, A5},
    {A6, A7, A8, A9, A10},
    {A11, A12, A13, A14, A15},
    {A16, A17, A18, A19, A20},
    {A21, A22, A23, A24, A25}
  };

  double Cmat[5][5];

  // diagonal entries always 1, off-diagonals: Aij*w + (1-w)
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      if (i == j) {
        Cmat[i][j] = 1;
      } else {
        Cmat[i][j] = A[i][j]*w + (1-w);
      }
    }
  }

  double dN_SI[5], dN_IR[5];
  for (i=0; i<n; i++) {
    double force = 0;
    for (j=0; j<n; j++) {
      force += Cmat[i][j]*I[j];
    }
    dN_SI[i] = rbinom(S[i], 1-exp(-beta_par/10000*force*dt));
    dN_IR[i] = rbinom(I[i], 1-exp(-gamma*dt));
    S[i] -= dN_SI[i];
    I[i] += dN_SI[i] - dN_IR[i];
    R[i] += dN_IR[i];
    H[i] += dN_SI[i];
  }

  S1=S[0]; S2=S[1]; S3=S[2]; S4=S[3]; S5=S[4];
  I1=I[0]; I2=I[1]; I3=I[2]; I4=I[3]; I5=I[4];
  R1=R[0]; R2=R[1]; R3=R[2]; R4=R[3]; R5=R[4];
  H1=H[0]; H2=H[1]; H3=H[2]; H4=H[3]; H5=H[4];
")



## ----Initial Conditions----------------------------------------------------------
rinit <- Csnippet("
  S1=5000; S2=5000; S3=5000; S4=5000; S5=5000;
  I1=0; I2=1; I3=0; I4=0; I5=0;
  R1=0; R2=0; R3=0; R4=0; R5=0;
  H1=0; H2=0; H3=0; H4=0; H5=0;
")

## ----Measurement Model Functions (Negative Binomial)--------------------------------------

rmeas <- Csnippet("
  double H[5] = {H1,H2,H3,H4,H5};
  reports1 = rnbinom_mu(k, H[0]*rho);
  reports2 = rnbinom_mu(k, H[1]*rho);
  reports3 = rnbinom_mu(k, H[2]*rho);
  reports4 = rnbinom_mu(k, H[3]*rho);
  reports5 = rnbinom_mu(k, H[4]*rho);
")

dmeas <- Csnippet("
  double H[5] = {H1,H2,H3,H4,H5};
  double reports[5] = {reports1,reports2,reports3,reports4,reports5};

  if (w < -1e-06 || w > 1+1e-06 || beta_par < 0 || k <= 0) {
    lik = R_NegInf;   // invalid parameters → -Inf
  } else {
    lik = 0;
    for (int i=0; i<5; i++) {
      // accumulate log-likelihood contributions
      lik += dnbinom_mu(reports[i], k, fmax(H[i]*rho, 1e-06), 1);
    }
  }
")





