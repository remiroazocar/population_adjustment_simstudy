# This file specifies the simulation setup

rm(list=ls())

# setwd("C:/Users/Antonio/Desktop/population_adjustment_simstudy") 

replicates <- 1000 # Monte Carlo replicates
allocation <- 1/2 # proportion of patients in active treatment (1:1 allocation ratio, vs. placebo)

no.chars <- 4 # number of baseline characteristics, these are prognostic variables
no.ems <- 2 # number of effect modifiers (subset of baseline characteristics)
pvs <- 1:4 # indices of prognostic variables
ems <- 1:2 # indices of effect modifiers 
N_AC <- c(150,300,600) # number of subjects in the AC trial
N_BC <- 600 # number of subjects in the BC trial 
b_trt <- log(0.25) # baseline effect (log HR) of active treatment vs. common comparator
b_X <- c(-log(0.67), -log(0.5), -log(0.33)) # effect of each prognostic variable
b_EM <- c(-log(0.67), -log(0.5), -log(0.33)) # interaction effect of each effect modifier
cens_prob <- 0.35 # censoring probability of time-to-event outcome
meanX_AC <- c(0.45,0.3,0.15) # mean of each normally-distributed covariate in AC trial
meanX_BC <- 0.6 # mean of each normally-distributed covariate in BC
sdX <- 0.2 # standard deviation of each covariate (same for AC and BC)
corX <- c(0,0.35) # covariate correlation coefficient  

# Weibull distribution parameters
weib_inv_scale <- 8.5 
weib_shape <- 1.3  
  
# parameter combinations for each scenario
param.combinations <- expand.grid(N_AC=N_AC, b_X=b_X, b_EM=b_EM, corX=corX, meanX_AC=meanX_AC)
pc <- param.combinations
scenarios <- nrow(param.combinations) # number of scenarios

# the rate parameter of the exponential distribution (from which censoring times are drawn from) 
# is selected to achieve specific censoring rate under active treatment at baseline. we select it 
# by simulating survival times and use optim to minimize difference between observed and target censoring rate
optim.function <- function(param, inv_scale, shape, b_trt, cens_prob, N) {
  cens_rate <- param # targeted censoring rate
  U <- runif(N)
  Tlat <- -log(U)/(inv_scale*exp(b_trt))^(1/shape) # latent survival time
  C <- rexp(n=N, rate=cens_rate) # censoring time
  prop_cens <- sum(Tlat>C)/N # observed censoring rate
  fit <- sum((prop_cens - cens_prob)^2) # minimize difference
  return(fit)
} 

cens_rate <- optim(par=1, fn=optim.function, inv_scale=weib_inv_scale, shape=weib_shape, b_trt=b_trt,
                   cens_prob=cens_prob, N=1000000, method="Brent", lower=0, upper=10)$par
save.image(file="survival_settings.RData") # save simulation settings