# This file performs the indirect comparison methods on the simulated data

rm(list=ls())

# setwd("C:/Users/Antonio/Desktop/population_adjustment_simstudy") 
load(file="survival_settings.RData") # load simulation settings
source('functions.R') # load MAIC functions

set.seed(444)
scenarios <- nrow(pc) # number of simulation scenarios

# for data manipulation 
if(!require(dplyr)) {install.packages("dplyr"); library(dplyr)}
if(!require(tidyr)) {install.packages("tidyr"); library(tidyr)}
# copula simulation
if(!require(simstudy)) {install.packages("simstudy"); library(simstudy)}
# for detect cores
if(!require(parallel)) {install.packages("parallel"); library(parallel)}
# for parallel cluster
if(!require(doSNOW)) {install.packages("doSNOW"); library(doSNOW)}
# to conduct standard cox regressions, e.g. for STC
if(!require(survival)) {install.packages("survival"); library(survival)}
# survey package for weighted survival 
if(!require(survey)) {install.packages("survey"); library(survey)}

# load simulated patient-level (AC) and aggregate-level (BC) datasets for all scenarios
IPD.AC.all <- vector(mode="list", scenarios)
ALD.BC.all <- vector(mode="list", scenarios)
for (i in 1:scenarios) {
  file.id <- paste0("N_AC", pc$N_AC[i], "b_X", round(pc$b_X[i], digits=2), 
                    "b_EM", round(pc$b_EM[i], digits=2),  
                    "prop_diff", pc$prop.diff[i], "corX", pc$corX[i]) 
  load(paste0("Data/IPD_AC_", file.id, ".RData"))
  load(paste0("Data/ALD_BC_", file.id, ".RData"))  
  IPD.AC.all[[i]] <- IPD.AC
  ALD.BC.all[[i]] <- ALD.BC
}

### Standard indirect treatment comparison (Bucher method)
bucher.wrapper <- function(data.AC, data.BC) {
  data.AC$trt <- factor(data.AC$trt, levels=c("C","A"))
  d.AC.bucher <- summary(coxph(Surv(time, status)~trt, data=data.AC))$coef[1]
  var.d.AC.bucher <- vcov(coxph(Surv(time, status)~trt, data=data.AC))[[1]]
  d.BC.bucher <- with(data.BC,logHR_B) 
  var.d.BC.bucher <- with(data.BC, var_logHR_B)
  d.AB.bucher <- unname(d.AC.bucher - d.BC.bucher) # average treatment effect
  var.d.AB.bucher <- unname(var.d.AC.bucher+var.d.BC.bucher) # variance of treatment effect
  list(d.AB.bucher, var.d.AB.bucher)
}  

### Matching-adjusted indirect comparison (MAIC)
maic.wrapper <- function(data.AC, data.BC, ems) { # ems indexes the position (columns) of the effect modifiers
  data.AC$trt <- factor(data.AC$trt, levels=c("C","A"))
  AC.ems <- data.AC[,1+ems] # column 1 of IPD is treatment indicator 
  maic.weights <- maic(A.X=AC.ems, B.summary=data.BC[ems]) # maic weights through method of moments
  maic.aess <- approx.ess(maic.weights) # approximate effective sample size
  # MAIC incorporated as weights in survey design 
  svy.design <- data.AC %>% mutate(wt=maic.weights) %>% svydesign(id=~1, weights=maic.weights, data=.)
  # Survey-weighted survival model
  outcome.fit.maic <- data.AC %>% svycoxph(Surv(time,status)~trt, design=svy.design, data=.)
  d.AC.maic <- summary(outcome.fit.maic)$coef[1]
  var.d.AC.maic <- vcov(outcome.fit.maic)[[1]] # sandwich-like variance estimator for A vs. C
  d.BC.maic <- with(data.BC,logHR_B)
  var.d.BC.maic <- with(data.BC, var_logHR_B)
  d.AB.maic <- d.AC.maic - d.BC.maic # ATE for A vs. B
  var.d.AB.maic <- var.d.AC.maic + var.d.BC.maic # Variance of A vs. B treatment effect
  list(d.AB.maic, var.d.AB.maic, maic.aess)
}

### Simulated treatment comparison (STC) - this is the "covariate simulation" approach
stc.wrapper <- function(data.AC, data.BC, pvs, ems) {
  data.AC$trt <- factor(data.AC$trt, levels=c("C","A"))
  # baseline characteristics/prognostic variables
  AC.chars <- data.AC[,1+pvs] # columns 1 of IPD is treatment indicator
  # cox regression on survival time performed in AC IPD is predictive model for outcome as 
  # a function of prognostic variables/effect modifiers. IPD covariates not centered in BC means
  stc.coxph <- coxph(as.formula(paste0("Surv(time,status)~",
                                       paste0(paste0(colnames(AC.chars)[pvs],
                                                     collapse="+"),"+",
                                              paste0(colnames(AC.chars)[ems],
                                                     "*trt",collapse="+")))),
                     data=data.AC)  
  S <- 100000 # Number of simulated patients in fake trial
  BC.sim <- data.frame(row.names=1:S) # data frame for simulated BC covariates 
  names.col <- NULL
  for (k in 1:length(pvs)) {
    name.col <- paste0("X", pvs[k]) 
    names.col <- c(names.col, name.col)
  }
  # simulate BC covariates from their approximate joint distribution using copula
  cor.AC <- cor(AC.chars[,pvs]) # AC trial correlation structure
  BC.sim <- as.data.frame(genCorGen(n=S, nvar=length(pvs), # copula
                          params1=as.numeric(data.BC[,pvs]), dist="binary", corstr="cs", 
                          corMatrix=cor.AC, wide=TRUE)[,-1])
  colnames(BC.sim) <- names.col
  BC.sim$trt <- NA
  # treatment assignment keeps A vs. C allocation ratio
  N_A <- floor(S*allocation) # number of patients under A in simulated dataset
  BC.sim$trt[1:N_A] <- "A" 
  BC.sim$trt[(N_A+1):S] <- "C" # under C in simulated
  BC.sim.trtA <- BC.sim[BC.sim$trt=="A",]
  BC.sim.trtC <- BC.sim[BC.sim$trt=="C",]
  # average linear predictions for each simulated patient under A and C
  predict.trtA <- predict(object=stc.coxph, newdata=BC.sim.trtA, type="lp", se.fit=TRUE) # A log HRs
  predict.trtC <- predict(object=stc.coxph, newdata=BC.sim.trtC, type="lp", se.fit=TRUE) # C log HRs
  logHR.AC <- mean(predict.trtA$fit)-mean(predict.trtC$fit) # A vs. C treatment effect
  var.logHR.AC <- mean(predict.trtA$se.fit^2)+mean(predict.trtC$se.fit^2) # A vs. C variance
  d.AC.stc <- logHR.AC
  var.d.AC.stc <- var.logHR.AC
  d.BC.stc <- with(data.BC,logHR_B)
  var.d.BC.stc <- with(data.BC, var_logHR_B)
  d.AB.stc <- d.AC.stc - d.BC.stc # A vs. B treatment effect
  var.d.AB.stc <- var.d.AC.stc + var.d.BC.stc # A vs. B variance
  list(d.AB.stc, var.d.AB.stc)
}  
  
### Simulated treatment comparison (STC-II) - this is the original "plug-in" approach
stc.ii.wrapper <- function(data.AC, data.BC, pvs, ems) {
  data.AC$trt <- factor(data.AC$trt, levels=c("C","A"))
  AC.chars <- data.AC[,1+pvs] # column 1 of IPD is treatment indicator
  # fit outcome regresion model with IPD effect modifiers centered at the mean BC values
  stc.coxph.ii <- coxph(as.formula(paste0("Surv(time,status)~", paste0(colnames(AC.chars)[pvs],
                                                                    collapse="+"),
                                       "+",  paste0("trt*I(", colnames(AC.chars)[ems],
                                                    "-", deparse(substitute(data.BC)),
                                                    "$prop.", colnames(AC.chars)[ems], ")",
                                                    collapse="+"))), data=data.AC)
  d.AC.stc.ii <- coef(stc.coxph.ii)["trtA"]
  var.d.AC.stc.ii <- vcov(stc.coxph.ii)["trtA", "trtA"]
  d.BC.stc.ii <- with(data.BC,logHR_B)
  var.d.BC.stc.ii <- with(data.BC, var_logHR_B)
  d.AB.stc.ii <- d.AC.stc.ii - d.BC.stc.ii # A vs. B treatment effect
  var.d.AB.stc.ii <- var.d.AC.stc.ii + var.d.BC.stc.ii # A vs. B variance
  list(d.AB.stc.ii, var.d.AB.stc.ii)  
}  
  
# set up cluster for parallel computing
num.cores <- detectCores()
cluster <- makeCluster(num.cores, type="SOCK", outfile="")
registerDoSNOW(cluster)
# progress bar
pb <- txtProgressBar(max=replicates, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

# combine lists in parallelisation
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

# run population adjustment methods for all replicates/scenarios in parallel
for(i in 1:scenarios) {
  IPD.AC <- IPD.AC.all[[i]]
  ALD.BC <- ALD.BC.all[[i]]
  pvs_i <- pvs
  ems_i <- ems
  file.id <- paste0("N_AC", pc$N_AC[i], "b_X", round(pc$b_X[i], digits=2), 
                    "b_EM", round(pc$b_EM[i], digits=2), 
                    "prop_diff", pc$prop.diff[i], "corX", pc$corX[i]) 
  bucher.results <- foreach(j=1:replicates, .combine='comb', .multicombine=TRUE,
                            .init=list(list(), list()), .options.snow=opts,
                            .packages=c("dplyr","tidyr", "survival")) %dopar% {
                              results <- bucher.wrapper(IPD.AC[[j]], ALD.BC[[j]])
                              return(results)
                            }
  close(pb)
  means <- unlist(bucher.results[[1]])
  variances <- unlist(bucher.results[[2]])
  save(means, file=paste0("Results/Bucher/means_", file.id, ".RData"))
  save(variances, file=paste0("Results/Bucher/variances_", file.id, ".RData"))  
  maic.results <- foreach(j=1:replicates, .combine='comb', .multicombine=TRUE,
                          .init=list(list(), list(), list()), .options.snow=opts,
                          .packages=c("dplyr","tidyr", "survey")) %dopar% {
                            results <- maic.wrapper(IPD.AC[[j]], ALD.BC[[j]], ems=ems_i)
                            return(results)
                          }
  close(pb)
  means <- unlist(maic.results[[1]])
  variances <- unlist(maic.results[[2]])
  approx.ess.maic <- unlist(maic.results[[3]])
  save(means, file=paste0("Results/MAIC/means_", file.id, ".RData"))
  save(variances, file=paste0("Results/MAIC/variances_", file.id, ".RData"))  
  save(approx.ess.maic, file=paste0("Results/MAIC/aess_", file.id, ".RData")) 
  stc.results <- foreach(j=1:replicates, .combine='comb', .multicombine=TRUE,
                         .init=list(list(), list()), .options.snow=opts,
                         .packages=c("survival", "simstudy")) %dopar% {
                           results <- stc.wrapper(IPD.AC[[j]], ALD.BC[[j]], pvs=pvs_i,
                                                  ems=ems_i)
                           return(results)
                         }
  close(pb)
  means <- unlist(stc.results[[1]])
  variances <- unlist(stc.results[[2]])
  save(means, file=paste0("Results/STC/means_", file.id, ".RData"))
  save(variances, file=paste0("Results/STC/variances_", file.id, ".RData"))
  stc.ii.results <- foreach(j=1:replicates, .combine='comb', .multicombine=TRUE,
                            .init=list(list(), list()), .options.snow=opts, 
                            .packages=c("survival")) %dopar% {
                            results <- stc.ii.wrapper(IPD.AC[[j]], ALD.BC[[j]], pvs=pvs_i, 
                                                      ems=ems_i)
                            return(results) 
                         }
  close(pb)
  means <- unlist(stc.ii.results[[1]])
  variances <- unlist(stc.ii.results[[2]])
  save(means, file=paste0("Results/STC_II/means_", file.id, ".RData"))
  save(variances, file=paste0("Results/STC_II/variances_", file.id, ".RData"))  
}  
 
stopCluster(cluster)
