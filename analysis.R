## This file processes the results of the simulation study and computes
## and graphs the relevant performance metrics. 

rm(list=ls())

# setwd("C:/Users/Antonio/Desktop/population_adjustment_simstudy")  
load(file="survival_settings.RData")
source("functions.R") # load functions to compute performance measures and for plotting

# container of means and variances for all replicates in each scenario  
maic.means.list <- vector(mode="list", scenarios)
maic.variances.list <- vector(mode="list", scenarios)
stc.means.list <- vector(mode="list", scenarios)
stc.variances.list <- vector(mode="list", scenarios)
stc.ii.means.list <- vector(mode="list", scenarios)
stc.ii.variances.list <- vector(mode="list", scenarios)
bucher.means.list <- vector(mode="list", scenarios)
bucher.variances.list <- vector(mode="list", scenarios)

# average treatment effect container 
maic.ate <- rep(NA, scenarios) 
stc.ate <- rep(NA, scenarios)
stc.ii.ate <- rep(NA, scenarios)
bucher.ate <- rep(NA, scenarios)
maic.ate.mcse <- rep(NA, scenarios) # monte carlo standard error
stc.ate.mcse <- rep(NA, scenarios)
stc.ii.ate.mcse <- rep(NA, scenarios) 
bucher.ate.mcse <- rep(NA, scenarios)

# lower bound of confidence interval container
maic.lci <- vector(mode="list", scenarios)
stc.lci <- vector(mode="list", scenarios)
stc.ii.lci <- vector(mode="list", scenarios)
bucher.lci <- vector(mode="list", scenarios)
maic.lci.mean <- rep(NA, scenarios)
stc.lci.mean <- rep(NA, scenarios)
stc.ii.lci.mean <- rep(NA, scenarios)
bucher.lci.mean <- rep(NA, scenarios)
maic.lci.mcse <- rep(NA, scenarios)
stc.lci.mcse <- rep(NA, scenarios)
stc.ii.lci.mcse <- rep(NA, scenarios)
bucher.lci.mcse <- rep(NA, scenarios)

# upper bound of confidence interval container
maic.uci <- vector(mode="list", scenarios)
stc.uci <- vector(mode="list", scenarios)
stc.ii.uci <- vector(mode="list", scenarios)
bucher.uci <- vector(mode="list", scenarios)
maic.uci.mean <- rep(NA, scenarios)
stc.uci.mean <- rep(NA, scenarios)
stc.ii.uci.mean <- rep(NA, scenarios)
bucher.uci.mean <- rep(NA, scenarios)
maic.uci.mcse <- rep(NA, scenarios)
stc.uci.mcse <- rep(NA, scenarios)
stc.ii.uci.mcse <- rep(NA, scenarios)
bucher.uci.mcse <- rep(NA, scenarios)

# Confidence interval widths
maic.ciwidth <- vector(mode="list", scenarios)  
stc.ciwidth <- vector(mode="list", scenarios)
stc.ii.ciwidth <- vector(mode="list", scenarios)
bucher.ciwidth <- vector(mode="list", scenarios)

# Bias containers
maic.bias <- rep(NA, scenarios) 
stc.bias <- rep(NA, scenarios)
stc.ii.bias <- rep(NA, scenarios)
bucher.bias <- rep(NA, scenarios)
maic.bias.mcse <- rep(NA, scenarios) 
stc.bias.mcse <- rep(NA, scenarios)
stc.ii.bias.mcse <- rep(NA, scenarios)
bucher.bias.mcse <- rep(NA, scenarios)

# mean absolute error (MAE) containers
maic.abs.err <- vector(mode="list", scenarios)
stc.abs.err <- vector(mode="list", scenarios)
stc.ii.abs.err <- vector(mode="list", scenarios)
bucher.abs.err <- vector(mode="list", scenarios)
maic.mae <- rep(NA, scenarios)
stc.mae <- rep(NA, scenarios)
stc.ii.mae <- rep(NA, scenarios)
bucher.mae <- rep(NA, scenarios)
maic.mae.mcse <- rep(NA, scenarios)
stc.mae.mcse <- rep(NA, scenarios)
stc.ii.mae.mcse <- rep(NA, scenarios)
bucher.mae.mcse <- rep(NA, scenarios)

# mean square error (MSE) containers
maic.mse <- rep(NA, scenarios)
stc.mse <- rep(NA, scenarios)
stc.ii.mse <- rep(NA, scenarios)
bucher.mse <- rep(NA, scenarios)
maic.mse.mcse <- rep(NA, scenarios)
stc.mse.mcse <- rep(NA, scenarios)
stc.ii.mse.mcse <- rep(NA, scenarios)
bucher.mse.mcse <- rep(NA, scenarios)

# variability ratio containers
maic.vr <- rep(NA, scenarios)
stc.vr <- rep(NA, scenarios)
stc.ii.vr <- rep(NA, scenarios)
bucher.vr <- rep(NA, scenarios)
maic.vr.mcse <- rep(NA, scenarios)
stc.vr.mcse <- rep(NA, scenarios)
stc.ii.vr.mcse <- rep(NA, scenarios)
bucher.vr.mcse <- rep(NA, scenarios)

# empirical standard error (EmpSE) containers
maic.empse <- rep(NA, scenarios)
stc.empse <- rep(NA, scenarios)
stc.ii.empse <- rep(NA, scenarios)
bucher.empse <- rep(NA, scenarios)
maic.empse.mcse <- rep(NA, scenarios)
stc.empse.mcse <- rep(NA, scenarios)
stc.ii.empse.mcse <- rep(NA, scenarios)
bucher.empse.mcse <- rep(NA, scenarios)

# coverage rate (%) of 95% confidence intervals
maic.cov <- rep(NA, scenarios)
stc.cov <- rep(NA, scenarios)
stc.ii.cov <- rep(NA, scenarios)
bucher.cov <- rep(NA, scenarios)
maic.cov.mcse <- rep(NA, scenarios)
stc.cov.mcse <- rep(NA, scenarios)
stc.ii.cov.mcse <- rep(NA, scenarios)
bucher.cov.mcse <- rep(NA, scenarios)

# % of replicates worse than Bucher for each scenario  
maic.error.worse.than <- rep(NA,scenarios)
maic.error.worse.than.mcse <- rep(NA,scenarios)
stc.error.worse.than <- rep(NA,scenarios)
stc.error.worse.than.mcse <- rep(NA,scenarios)
stc.ii.error.worse.than <- rep(NA,scenarios)
stc.ii.error.worse.than.mcse <- rep(NA,scenarios)

# standardized biases
maic.std.bias <- rep(NA,scenarios) 
stc.std.bias <- rep(NA,scenarios)
stc.ii.std.bias <- rep(NA,scenarios)
bucher.std.bias <- rep(NA,scenarios)

# approximate effective sample sizes
maic.aess.mean <- rep(NA,scenarios)
maic.aess.mcse <- rep(NA,scenarios)

# number of replicates for which STC does not converge 
stc.notconverges.list <- rep(NA, scenarios)
stc.ii.notconverges.list <- rep(NA, scenarios)
# assume that STC has not converged if variance greater than this number
max.variance <- 5  

# table storing parameter settings and performance measures for each scenario
scenarios.df <- data.frame()

for (i in 1:scenarios) {
  file.id <- paste0("N_AC", pc$N_AC[i], "b_X", round(pc$b_X[i], digits=2), 
                    "b_EM", round(pc$b_EM[i], digits=2),
                    "prop_diff", pc$prop.diff[i], "corX", pc$corX[i])   
  # MAIC
  load(paste0("Results/MAIC/means_", file.id, ".RData"))
  load(paste0("Results/MAIC/variances_", file.id, ".RData"))
  ### MATCHING-ADJUSTED INDIRECT COMPARISON
  maic.means.list[[i]] <- means
  maic.variances.list[[i]] <- variances
  maic.bias[i] <- bias(maic.means.list[[i]], b_trt-b_trt)  
  maic.bias.mcse[i] <- bias.mcse(maic.means.list[[i]])
  maic.mae[i] <- mae(maic.means.list[[i]], b_trt-b_trt)
  maic.abs.err[[i]] <- maic.means.list[[i]] - (b_trt-b_trt)
  maic.mae.mcse[i] <- mcse.estimate(maic.abs.err[[i]])
  maic.mse[i] <- mse(maic.means.list[[i]], b_trt-b_trt) 
  maic.mse.mcse[i] <- mse.mcse(maic.means.list[[i]], b_trt-b_trt) 
  maic.ate[i] <- mean(maic.means.list[[i]])
  maic.ate.mcse[i] <- mcse.estimate(maic.means.list[[i]])
  # construct confidence interval using normal distribution
  maic.lci[[i]] <- maic.means.list[[i]] + qnorm(0.025)*sqrt(maic.variances.list[[i]])
  maic.uci[[i]] <- maic.means.list[[i]] + qnorm(0.975)*sqrt(maic.variances.list[[i]])
  maic.ciwidth[[i]] <- maic.uci[[i]] - maic.lci[[i]]
  maic.lci.mean[i] <- mean(maic.lci[[i]])
  maic.lci.mcse[i] <- mcse.estimate(maic.lci[[i]])
  maic.uci.mean[i] <- mean(maic.uci[[i]])
  maic.uci.mcse[i] <- mcse.estimate(maic.uci[[i]])
  maic.cov[i] <- coverage(maic.lci[[i]], maic.uci[[i]], b_trt-b_trt)
  maic.cov.mcse[i] <- coverage.mcse(maic.cov[i], length(maic.lci[[i]]))
  maic.empse[i] <- empse(maic.means.list[[i]])
  maic.empse.mcse[i] <- empse.mcse(maic.empse[i], length(maic.means.list[[i]]))
  maic.vr[i] <- var.ratio(maic.means.list[[i]], sqrt(maic.variances.list[[i]]))
  maic.vr.mcse[i] <- var.ratio.mcse(avg.se=mean(sqrt(maic.variances.list[[i]])), 
                                    emp.se=maic.empse[i],
                                    var.avg.se=mcse.estimate(maic.variances.list[[i]])^2,
                                    var.emp.se=maic.empse.mcse[i]^2)
  maic.std.bias[i] <- (maic.bias[i]*100)/maic.empse[i]
  load(paste0("Results/MAIC/aess_", file.id, ".RData"))
  maic.aess.mean[i] <- mean(approx.ess.maic)
  maic.aess.mcse[i] <- mcse.estimate(approx.ess.maic) 
  ### SIMULATED TREATMENT COMPARISON ("covariate simulation" approach)
  load(paste0("Results/STC/means_", file.id, ".RData"))
  load(paste0("Results/STC/variances_", file.id, ".RData"))  
  stc.notconverges.list[i] <- sum(variances>max.variance)
  # discard replicates for which STC did not converge
  stc.converges <- variances<max.variance
  means <- means[stc.converges]
  variances <- variances[stc.converges]
  stc.means.list[[i]] <- means
  stc.variances.list[[i]] <- variances
  stc.bias[i] <- bias(stc.means.list[[i]], b_trt-b_trt)  
  stc.bias.mcse[i] <- bias.mcse(stc.means.list[[i]])
  stc.mae[i] <- mae(stc.means.list[[i]], b_trt-b_trt)
  stc.mae.mcse[i] <- mcse.estimate(stc.means.list[[i]])
  stc.abs.err[[i]] <- stc.means.list[[i]] - (b_trt-b_trt)
  stc.mae.mcse[i] <- mcse.estimate(stc.abs.err[[i]])
  stc.mse[i] <- mse(stc.means.list[[i]], b_trt-b_trt) 
  stc.mse.mcse[i] <- mse.mcse(stc.means.list[[i]], b_trt-b_trt) 
  stc.ate[i] <- mean(stc.means.list[[i]]) 
  stc.ate.mcse[i] <- mcse.estimate(stc.means.list[[i]])
  stc.lci[[i]] <- stc.means.list[[i]] + qnorm(0.025)*sqrt(stc.variances.list[[i]])
  stc.uci[[i]] <- stc.means.list[[i]] + qnorm(0.975)*sqrt(stc.variances.list[[i]])
  stc.ciwidth[[i]] <- stc.uci[[i]] - stc.lci[[i]]
  stc.lci.mean[i] <- mean(stc.lci[[i]])
  stc.lci.mcse[i] <- mcse.estimate(stc.lci[[i]])
  stc.uci.mean[i] <- mean(stc.uci[[i]])
  stc.uci.mcse[i] <- mcse.estimate(stc.uci[[i]])
  stc.cov[i] <- coverage(stc.lci[[i]], stc.uci[[i]], b_trt-b_trt)
  stc.cov.mcse[i] <- coverage.mcse(stc.cov[i], length(stc.lci[[i]]))
  stc.empse[i] <- empse(stc.means.list[[i]])
  stc.empse.mcse[i] <- empse.mcse(stc.empse[i], length(stc.means.list[[i]]))
  stc.vr[i] <- var.ratio(stc.means.list[[i]], sqrt(stc.variances.list[[i]]))
  stc.vr.mcse[i] <- var.ratio.mcse(avg.se=mean(sqrt(stc.variances.list[[i]])), 
                                   emp.se=stc.empse[i],
                                   var.avg.se=mcse.estimate(stc.variances.list[[i]])^2,
                                   var.emp.se=stc.empse.mcse[i]^2)  
  stc.std.bias[i] <- (stc.bias[i]*100)/stc.empse[i]
  ### SIMULATED TREATMENT COMPARISON ("plug-in" approach)
  load(paste0("Results/STC_II/means_", file.id, ".RData"))
  load(paste0("Results/STC_II/variances_", file.id, ".RData"))  
  stc.ii.notconverges.list[i] <- sum(variances>max.variance)
  # discard replicates for which STC did not converge
  stc.ii.converges <- variances<max.variance
  means <- means[stc.ii.converges]
  variances <- variances[stc.ii.converges]
  stc.ii.means.list[[i]] <- means
  stc.ii.variances.list[[i]] <- variances
  stc.ii.bias[i] <- bias(stc.ii.means.list[[i]], b_trt-b_trt)  
  stc.ii.bias.mcse[i] <- bias.mcse(stc.ii.means.list[[i]])
  stc.ii.mae[i] <- mae(stc.ii.means.list[[i]], b_trt-b_trt)
  stc.ii.mae.mcse[i] <- mcse.estimate(stc.ii.means.list[[i]])
  stc.ii.abs.err[[i]] <- stc.ii.means.list[[i]] - (b_trt-b_trt)
  stc.ii.mae.mcse[i] <- mcse.estimate(stc.ii.abs.err[[i]])
  stc.ii.mse[i] <- mse(stc.ii.means.list[[i]], b_trt-b_trt) 
  stc.ii.mse.mcse[i] <- mse.mcse(stc.ii.means.list[[i]], b_trt-b_trt) 
  stc.ii.ate[i] <- mean(stc.ii.means.list[[i]]) 
  stc.ii.ate.mcse[i] <- mcse.estimate(stc.ii.means.list[[i]])
  stc.ii.lci[[i]] <- stc.ii.means.list[[i]] + qnorm(0.025)*sqrt(stc.ii.variances.list[[i]])
  stc.ii.uci[[i]] <- stc.ii.means.list[[i]] + qnorm(0.975)*sqrt(stc.ii.variances.list[[i]])
  stc.ii.ciwidth[[i]] <- stc.ii.uci[[i]] - stc.ii.lci[[i]]
  stc.ii.lci.mean[i] <- mean(stc.ii.lci[[i]])
  stc.ii.lci.mcse[i] <- mcse.estimate(stc.ii.lci[[i]])
  stc.ii.uci.mean[i] <- mean(stc.ii.uci[[i]])
  stc.ii.uci.mcse[i] <- mcse.estimate(stc.ii.uci[[i]])
  stc.ii.cov[i] <- coverage(stc.ii.lci[[i]], stc.ii.uci[[i]], b_trt-b_trt)
  stc.ii.cov.mcse[i] <- coverage.mcse(stc.ii.cov[i], length(stc.ii.lci[[i]]))
  stc.ii.empse[i] <- empse(stc.ii.means.list[[i]])
  stc.ii.empse.mcse[i] <- empse.mcse(stc.ii.empse[i], length(stc.ii.means.list[[i]]))
  stc.ii.vr[i] <- var.ratio(stc.ii.means.list[[i]], sqrt(stc.ii.variances.list[[i]]))
  stc.ii.vr.mcse[i] <- var.ratio.mcse(avg.se=mean(sqrt(stc.ii.variances.list[[i]])), 
                                   emp.se=stc.ii.empse[i],
                                   var.avg.se=mcse.estimate(stc.ii.variances.list[[i]])^2,
                                   var.emp.se=stc.ii.empse.mcse[i]^2)  
  stc.ii.std.bias[i] <- (stc.ii.bias[i]*100)/stc.ii.empse[i]    
  ### BUCHER METHOD (STANDARD INDIRECT COMPARISON)
  load(paste0("Results/Bucher/means_", file.id, ".RData"))
  load(paste0("Results/Bucher/variances_", file.id, ".RData")) 
  bucher.means.list[[i]] <- means
  bucher.variances.list[[i]] <- variances
  bucher.bias[i] <- bias(bucher.means.list[[i]], b_trt-b_trt)  
  bucher.bias.mcse[i] <- bias.mcse(bucher.means.list[[i]])
  bucher.mae[i] <- mae(bucher.means.list[[i]], b_trt-b_trt)
  bucher.abs.err[[i]] <- bucher.means.list[[i]] - (b_trt-b_trt)
  bucher.mae.mcse[i] <- mcse.estimate(bucher.abs.err[[i]])
  bucher.mse[i] <- mse(bucher.means.list[[i]], b_trt-b_trt) 
  bucher.mse.mcse[i] <- mse.mcse(bucher.means.list[[i]], b_trt-b_trt) 
  bucher.ate[i] <- mean(bucher.means.list[[i]]) 
  bucher.ate.mcse[i] <- mcse.estimate(bucher.means.list[[i]])
  bucher.lci[[i]] <- bucher.means.list[[i]] + qnorm(0.025)*sqrt(bucher.variances.list[[i]])
  bucher.uci[[i]] <- bucher.means.list[[i]] + qnorm(0.975)*sqrt(bucher.variances.list[[i]])
  bucher.ciwidth[[i]] <- bucher.uci[[i]] - bucher.lci[[i]]
  bucher.lci.mean[i] <- mean(bucher.lci[[i]])
  bucher.lci.mcse[i] <- mcse.estimate(bucher.lci[[i]])
  bucher.uci.mean[i] <- mean(bucher.uci[[i]])
  bucher.uci.mcse[i] <- mcse.estimate(bucher.uci[[i]])
  bucher.cov[i] <- coverage(bucher.lci[[i]], bucher.uci[[i]], b_trt-b_trt)
  bucher.cov.mcse[i] <- coverage.mcse(bucher.cov[i], length(bucher.lci[[i]]))
  bucher.empse[i] <- empse(bucher.means.list[[i]])
  bucher.empse.mcse[i] <- empse.mcse(bucher.empse[i], replicates)
  bucher.vr[i] <- var.ratio(bucher.means.list[[i]], sqrt(bucher.variances.list[[i]]))
  bucher.vr.mcse[i] <- var.ratio.mcse(avg.se=mean(sqrt(bucher.variances.list[[i]])), 
                                      emp.se=bucher.empse[i],
                                      var.avg.se=mcse.estimate(bucher.variances.list[[i]])^2,
                                      var.emp.se=bucher.empse.mcse[i]^2) 
  bucher.std.bias[i] <- (bucher.bias[i]*100)/bucher.empse[i]
  truth <- b_trt-b_trt # true baseline A vs. B treatment effect
  maic.error.worse.than[i] <- sum(abs(bucher.means.list[[i]]-truth)<abs(maic.means.list[[i]]-truth))/replicates
  maic.error.worse.than.mcse[i] <- coverage.mcse(maic.error.worse.than[i], replicates)
  stc.error.worse.than[i] <- sum(abs(bucher.means.list[[i]][stc.converges]-truth)<abs(stc.means.list[[i]]-truth))/sum(stc.converges)
  stc.error.worse.than.mcse[i] <- coverage.mcse(stc.error.worse.than[i], sum(stc.converges))  
  stc.ii.error.worse.than[i] <- sum(abs(bucher.means.list[[i]][stc.ii.converges]-truth)<abs(stc.ii.means.list[[i]]-truth))/sum(stc.ii.converges)
  stc.ii.error.worse.than.mcse[i] <- coverage.mcse(stc.ii.error.worse.than[i], sum(stc.ii.converges)) 

  tmp.scenarios <- cbind(i, pc$N_AC[i], pc$b_X[i], pc$b_EM[i], pc$prop.diff[i], pc$corX[i])

  maic.tmp.metrics <- cbind(maic.ate[i], maic.ate.mcse[i], maic.lci.mean[i],
                            maic.lci.mcse[i], maic.uci.mean[i], maic.uci.mcse[i],
                            maic.bias[i], maic.bias.mcse[i], maic.mse[i], maic.mse.mcse[i],
                            maic.mae[i], maic.mae.mcse[i], maic.cov[i], maic.cov.mcse[i],
                            maic.empse[i], maic.empse.mcse[i], maic.vr[i], maic.vr.mcse[i],
                            maic.error.worse.than[i], maic.error.worse.than.mcse[i], 
                            maic.std.bias[i], maic.aess.mean[i], maic.aess.mcse[i])
  
  stc.tmp.metrics <- cbind(stc.ate[i], stc.ate.mcse[i], stc.lci.mean[i],
                           stc.lci.mcse[i], stc.uci.mean[i], stc.uci.mcse[i],
                           stc.bias[i], stc.bias.mcse[i], stc.mse[i], stc.mse.mcse[i],
                           stc.mae[i], stc.mae.mcse[i], stc.cov[i], stc.cov.mcse[i],
                           stc.empse[i], stc.empse.mcse[i], stc.vr[i], stc.vr.mcse[i],
                           stc.error.worse.than[i], stc.error.worse.than.mcse[i],
                           stc.std.bias[i], stc.notconverges.list[i])
  
  stc.ii.tmp.metrics <- cbind(stc.ii.ate[i], stc.ii.ate.mcse[i], stc.ii.lci.mean[i],
                              stc.ii.lci.mcse[i], stc.ii.uci.mean[i], stc.ii.uci.mcse[i],
                              stc.ii.bias[i], stc.ii.bias.mcse[i], stc.ii.mse[i], stc.ii.mse.mcse[i],
                              stc.ii.mae[i], stc.ii.mae.mcse[i], stc.ii.cov[i], stc.ii.cov.mcse[i],
                              stc.ii.empse[i], stc.ii.empse.mcse[i], stc.ii.vr[i], stc.ii.vr.mcse[i],
                              stc.ii.error.worse.than[i], stc.ii.error.worse.than.mcse[i],
                              stc.ii.std.bias[i], stc.ii.notconverges.list[i])

  bucher.tmp.metrics <- cbind(bucher.ate[i], bucher.ate.mcse[i], bucher.lci.mean[i],
                              bucher.lci.mcse[i], bucher.uci.mean[i], bucher.uci.mcse[i],
                              bucher.bias[i], bucher.bias.mcse[i], bucher.mse[i], bucher.mse.mcse[i],
                              bucher.mae[i], bucher.mae.mcse[i], bucher.cov[i], bucher.cov.mcse[i],
                              bucher.empse[i], bucher.empse.mcse[i], bucher.vr[i], bucher.vr.mcse[i], 
                              bucher.std.bias[i]) 
  
  tmp.scenarios<- cbind(tmp.scenarios, maic.tmp.metrics, stc.tmp.metrics, stc.ii.tmp.metrics,
                        bucher.tmp.metrics)
  
  scenarios.df <- rbind(scenarios.df, tmp.scenarios)
}  

colnames(scenarios.df) <- c("Scenario", "N_AC", "b_X","b_EM","prop_diff", "corX",
                            "maic.ate", "maic.ate.mcse", "maic.lci.mean",
                            "maic.lci.mcse", "maic.uci.mean", "maic.uci.mcse",
                            "maic.bias", "maic.bias.mcse", "maic.mse", "maic.mse.mcse",
                            "maic.mae", "maic.mae.mcse", "maic.cov", "maic.cov.mcse",
                            "maic.empse", "maic.empse.mcse", "maic.vr", "maic.vr.mcse",
                            "maic.error.worse.than", "maic.error.worse.than.mcse", 
                            "maic.std.bias", "maic.aess.mean", "maic.aess.mcse", 
                            "stc.ate", "stc.ate.mcse", 
                            "stc.lci.mean", "stc.lci.mcse", "stc.uci.mean", "stc.uci.mcse",
                            "stc.bias", "stc.bias.mcse", "stc.mse", "stc.mse.mcse",
                            "stc.mae", "stc.mae.mcse", "stc.cov", "stc.cov.mcse",
                            "stc.empse", "stc.empse.mcse", "stc.vr", "stc.vr.mcse",
                            "stc.error.worse.than", "stc.error.worse.than.mcse",
                            "stc.std.bias","stc.no.convergence", 
                            "stc.ii.ate", "stc.ii.ate.mcse", 
                            "stc.ii.lci.mean", "stc.ii.lci.mcse", "stc.ii.uci.mean", "stc.ii.uci.mcse",
                            "stc.ii.bias", "stc.ii.bias.mcse", "stc.ii.mse", "stc.ii.mse.mcse",
                            "stc.ii.mae", "stc.ii.mae.mcse", "stc.ii.cov", "stc.ii.cov.mcse",
                            "stc.ii.empse", "stc.ii.empse.mcse", "stc.ii.vr", "stc.ii.vr.mcse",
                            "stc.ii.error.worse.than", "stc.ii.error.worse.than.mcse",
                            "stc.ii.std.bias","stc.ii.no.convergence",
                            "bucher.ate", "bucher.ate.mcse", "bucher.lci.mean",
                            "bucher.lci.mcse", "bucher.uci.mean", "bucher.uci.mcse",
                            "bucher.bias", "bucher.bias.mcse", "bucher.mse", "bucher.mse.mcse",
                            "bucher.mae", "bucher.mae.mcse", "bucher.cov", "bucher.cov.mcse",
                            "bucher.empse", "bucher.empse.mcse", "bucher.vr", 
                            "bucher.vr.mcse", "bucher.std.bias")

write.csv(scenarios.df, "Analysis/scenarios.csv", row.names = FALSE)

scenarios.df <- read.csv(file="Analysis/scenarios.csv", header=TRUE, sep=",")

scenarios.df$b_EM <- round(scenarios.df$b_EM, 2)
scenarios.df$b_X <- round(scenarios.df$b_X, 2)

# reorder simulation results in order to be presented in nested loop plot (FOR BIAS).   
# (function is by Rücker, G., Schwarzer, G. Presenting simulation results in a nested loop plot. 
# BMC Med Res Methodol 14, 129 (2014) doi:10.1186/1471-2288-14-129)
nldata <- nestedloop(scenarios.df,
                     varnames=c("b_EM","b_X","N_AC","prop_diff", "corX"),
                     varlabels=c("Effect-modifying interaction",
                                 "Prognostic variable effect", "Subjects in trial with patient-level data",
                                 "Covariate overlap","Covariate correlation"),
                     sign=c(1, 1, 1, 1, 1))
pd.bias <- nldata # object to be used in nested loop plot
# use labels instead of numeric values for the following factors
pd.bias$prop_diff <- factor(pd.bias$prop_diff, levels=c(0.2,0.3,0.4), 
                            labels=c("strong","moderate","poor"))
pd.bias$corX <- factor(pd.bias$corX, levels=c(0,0.35),labels=c("none","moderate"))
pd.bias$b_EM <- factor(pd.bias$b_EM,levels=c(0.40,0.69,1.11), labels=c("moderate","strong",
                                                                       "very strong"))
pd.bias$b_X <- factor(pd.bias$b_X, levels=c(0.40,0.69,1.11), labels=c("moderate","strong",
                                                                      "very strong"))
# reorder simulation results in nested loop plot presentation order (FOR OTHER PERFORMANCE MEASURES)
nldata <- nestedloop(scenarios.df,
                     varnames=c("b_X","N_AC","b_EM","prop_diff", "corX"),
                     varlabels=
                       c("Prognostic variable effect", "Subjects in trial with patient-level data",
                         "Effect-modifiying interaction", "Covariate overlap",
                         "Covariate correlation"),
                     sign=c(1, 1, 1, 1, 1))
pd.other <- nldata
pd.other$prop_diff <- factor(pd.other$prop_diff,levels=c(0.2,0.3,0.4),
                             labels=c("strong","moderate","poor"))
pd.other$corX <- factor(pd.other$corX, levels=c(0,0.35),labels=c("none","moderate"))
pd.other$b_EM <- factor(pd.other$b_EM,levels=c(0.40,0.69,1.11),labels=c("moderate","strong",
                                                                        "very strong"))
pd.other$b_X <- factor(pd.other$b_X,levels=c(0.40,0.69,1.11),labels=c("moderate","strong",
                                                                       "very strong"))

### MAIN MANUSCRIPT PLOTS (THESE DO NOT INCLUDE STC-II) ###

### Nested loop plot for MSE
# pdf("Analysis/mse.pdf", width=18,height=15,pointsize=20)
tiff("Analysis/mse.tiff", res=800, width=18, height=15, units='in')
par(pty="m")
plot(pd.other$maic.mse,
     type="n",
     ylim=c(0,0.5), bty="n",
     xlab="Scenario",
     ylab="Mean square error (MSE)",
     las=1, xaxt="n",
     cex.axis=1.5, cex.lab=1.5) #1.15
lines(pd.other, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd.other, which="r", ymin.refline=0.375, ymax.refline=0.5, cex.ref=1.5) #0.9
lines(pd.other$maic.mse, col="red", lty=2,type="s", lwd=1.5)
lines(pd.other$stc.mse, col="green", lty=5, type="s", lwd=1.5)
lines(pd.other$bucher.mse, col="blue", lty=4, type="s", lwd=1.5)
legend(x=130, y=0.35,lwd=c(1.5,1.5,1.5),col=c("red","green","blue"), lty=c(2,5,4),
       cex=1.5,bty="n",c("MAIC", "STC","Bucher"))
dev.off()

### Nested loop plot for coverage
# pdf("Analysis/coverage.pdf", width=18,height=15,pointsize=20)
tiff("Analysis/coverage.tiff", res=800, width=8.64, height=7.2, units='in')
par(pty="m")
plot(pd.other$maic.cov*100, type="n", ylim=c(20, 100), bty="n", xlab="Scenario",
     ylab="Coverage of 95% confidence intervals (%)", las=1, xaxt="n",
     cex.axis=1, cex.lab=1)
abline(h=95, col="grey") # nominal 
lines(pd.other, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd.other, which="r", ymin.refline=20, ymax.refline=40, cex.ref=0.9)
lines(pd.other$maic.cov*100, col="red", lty=2,type="s", lwd=1)    
lines(pd.other$stc.cov*100, col="green", lty=5, type="s", lwd=1) 
lines(pd.other$bucher.cov*100, col="blue", lty=4, type="s", lwd=1) 
legend("left", lwd=c(1.5,1.5,1.5), col=c("red","green", "blue"), 
       lty=c(2,5,4), cex=1, bty="n",c("MAIC", "STC", "Bucher"))
dev.off()

### Nested loop plot for EmpSE 
# pdf("Analysis/EmpSE.pdf",width=18, height=15,pointsize=20)
tiff("Analysis/EmpSE.tiff", res=800, width=18, height=15, units='in')
par(pty="m")
plot(pd.other$maic.empse, type="n",ylim=c(0, 0.4), bty="n", xlab="Scenario",
     ylab="Empirical standard error (ESE)",las=1, xaxt="n",
     cex.axis=1.5, cex.lab=1.5)
lines(pd.other, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd.other, which="r", ymin.refline=0.3, ymax.refline=0.4, cex.ref=1.5)
lines(pd.other$maic.empse, col="red", lty=2,type="s",lwd=1.5)    
lines(pd.other$stc.empse, col="green", lty=5, type="s", lwd=1.5)   
lines(pd.other$bucher.empse, col="blue", lty=4, type="s", lwd=1.5) 
legend("bottom", lwd=c(1.5,1.5,1.5), col=c("red","green","blue"), lty=c(2,5,4), 
       cex=1.5, bty="n", c("MAIC", "STC", "Bucher"))
dev.off()

### Nested loop plot for variability ratio
# pdf("Analysis/variability_ratio.pdf", width=18,height=15, pointsize=20)
tiff("Analysis/variability_ratio.tiff", res=800, width=18, height=15, units='in')
par(pty="m")
plot(pd.other$maic.vr, type="n", ylim=c(0.8,2.2), bty="n", xlab="Scenario", ylab="Variability ratio",
     las=1, xaxt="n", cex.axis=1.5, cex.lab=1.5)
abline(h=1, col="grey") # unbiased variance
lines(pd.other, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd.other, which="r", ymin.refline=1.85, ymax.refline=2.2, cex.ref=1.5)
lines(pd.other$maic.vr, col="red", lty=2,type="s",lwd=1.5)    
lines(pd.other$stc.vr, col="green", lty=5, type="s", lwd=1.5)   
lines(pd.other$bucher.vr, col="blue", lty=4, type="s", lwd=1.5) 
legend(x=77, y=1.29, lwd=c(1.5,1.5,1.5), col=c("red","green","blue"),lty=c(2,5,4), 
       cex=1.5, bty="n", c("MAIC", "STC", "Bucher"))
dev.off()
 
### Nested loop plot for bias 
# pdf("Analysis/bias.pdf", width=18, height=15, pointsize=20)
tiff("Analysis/bias.tiff", res=800, width=18, height=15, units='in')
par(pty="m")
plot(pd.bias$maic.bias, type="n", ylim=c(-0.75,0.75), bty="n", xlab="Scenario", ylab="Bias",
     las=1, xaxt="n", cex.axis=1.5, cex.lab=1.5)
abline(h=0, col="grey") # no bias
lines(pd.bias, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE")) # add vertical lines
lines(pd.bias, which="r",ymin.refline=0.375, ymax.refline=0.75, cex.ref=1.5) # add reference lines
lines(pd.bias$maic.bias, col="red", lty=2,type="s", lwd=1.5) # performance measures 
lines(pd.bias$stc.bias, col="green", lty=5, type="s", lwd=1.5)
lines(pd.bias$bucher.bias, col="blue", lty=4, type="s", lwd=1.5)
legend("bottom",lwd=c(1.5,1.5,1.5),col=c("red","green", "blue"),
       lty=c(2,5,4),cex=1.5,bty="n", c("MAIC", "STC", "Bucher")) # legend
dev.off()

#### SUPPLEMENTARY MATERIAL PLOTS (INCLUDE STC-II) ####

### Nested loop plot for MSE
pdf("Analysis/Supplementary_Material/mse.pdf", width=18,height=15,pointsize=20)
par(pty="m")
plot(pd.other$maic.mse,
     type="n",
     ylim=c(0,0.6), bty="n",
     xlab="Scenario",
     ylab="Mean square error (MSE)",
     las=1, xaxt="n",
     cex.axis=1.15, cex.lab=1.15)
lines(pd.other, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd.other, which="r", ymin.refline=0.45, ymax.refline=0.6, cex.ref=0.9) 
lines(pd.other$maic.mse, col="red", lty=2,type="s", lwd=1.5)
lines(pd.other$stc.mse, col="green", lty=5, type="s", lwd=1.5)
lines(pd.other$stc.ii.mse, col="pink", lty=1, type="s", lwd=1.5)
lines(pd.other$bucher.mse, col="blue", lty=4, type="s", lwd=1.5)
legend(x=130, y=0.35, lwd=c(1.5,1.5,1.5,1.5),col=c("red","green","pink","blue"), lty=c(2,5,1,4),
       cex=1,bty="n",c("MAIC", "STC","STC-II","Bucher"))
dev.off()

### Nested loop plot for coverage
pdf("Analysis/Supplementary_Material/coverage.pdf", width=18,height=15,pointsize=20)
par(pty="m")
plot(pd.other$maic.cov*100, type="n", ylim=c(20, 100), bty="n", xlab="Scenario",
     ylab="Coverage of 95% confidence intervals (%)", las=1, xaxt="n",
     cex.axis=1.15, cex.lab=1.15)
abline(h=95, col="grey") # nominal 
lines(pd.other, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd.other, which="r", ymin.refline=20, ymax.refline=40, cex.ref=0.9)
lines(pd.other$maic.cov*100, col="red", lty=2,type="s", lwd=1.5)    
lines(pd.other$stc.cov*100, col="green", lty=5, type="s", lwd=1.5)   
lines(pd.other$stc.ii.cov*100, col="pink", lty=1, type="s", lwd=1.5) 
lines(pd.other$bucher.cov*100, col="blue", lty=4, type="s", lwd=1.5) 
legend("left",lwd=c(1.5,1.5,1.5,1.5),col=c("red","green","pink", "blue"), lty=c(2,5,1,4),
       cex=0.9, bty="n",c("MAIC", "STC", "STC-II", "Bucher"))
dev.off()

### Nested loop plot for EmpSE 
pdf("Analysis/Supplementary_Material/EmpSE.pdf",width=18, height=15,pointsize=20)
par(pty="m")
plot(pd.other$maic.empse, type="n",ylim=c(0, 0.4), bty="n", xlab="Scenario",
     ylab="Empirical standard error (ESE)",las=1, xaxt="n",
     cex.axis=1.15, cex.lab=1.15)
lines(pd.other, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd.other, which="r", ymin.refline=0.3, ymax.refline=0.4, cex.ref=0.9)
lines(pd.other$maic.empse, col="red", lty=2,type="s",lwd=1.5)    
lines(pd.other$stc.empse, col="green", lty=5, type="s", lwd=1.5)
lines(pd.other$stc.ii.empse, col="pink", lty=1, type="s", lwd=1.5) 
lines(pd.other$bucher.empse, col="blue", lty=4, type="s", lwd=1.5) 
legend("bottom", lwd=c(1.5,1.5,1.5,1.5), col=c("red","green","pink","blue"), lty=c(2,5,1,4), 
       cex=0.9, bty="n", c("MAIC", "STC", "STC-II", "Bucher"))
dev.off()

### Nested loop plot for variability ratio
pdf("Analysis/Supplementary_Material/variability_ratio.pdf", width=18,height=15, pointsize=20)
par(pty="m")
plot(pd.other$maic.vr, type="n", ylim=c(0.8,2.2), bty="n", xlab="Scenario", ylab="Variability ratio",
     las=1, xaxt="n", cex.axis=1.15, cex.lab=1.15)
abline(h=1, col="grey") # unbiased variance
lines(pd.other, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd.other, which="r", ymin.refline=1.85, ymax.refline=2.2, cex.ref=0.9)
lines(pd.other$maic.vr, col="red", lty=2,type="s",lwd=1.5)    
lines(pd.other$stc.vr, col="green", lty=5, type="s", lwd=1.5) 
lines(pd.other$stc.ii.vr, col="pink", lty=1, type="s", lwd=1.5) 
lines(pd.other$bucher.vr, col="blue", lty=4, type="s", lwd=1.5) 
legend(x=77, y=1.29, lwd=c(1.5,1.5,1.5,1.5), col=c("red","green","pink","blue"),lty=c(2,5,1,4), 
       cex=0.9, bty="n", c("MAIC", "STC","STC-II","Bucher"))
dev.off()

### Nested loop plot for bias 
pdf("Analysis/Supplementary_Material/bias.pdf", width=18, height=15, pointsize=20)
par(pty="m")
plot(pd.bias$maic.bias, type="n", ylim=c(-0.75,0.75), bty="n", xlab="Scenario", ylab="Bias",
     las=1, xaxt="n", cex.axis=1.15, cex.lab=1.15)
abline(h=0, col="grey") # no bias
lines(pd.bias, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE")) # add vertical lines
lines(pd.bias, which="r",ymin.refline=0.375, ymax.refline=0.75, cex.ref=0.9) # add reference lines
lines(pd.bias$maic.bias, col="red", lty=2,type="s", lwd=1.5) # performance measures 
lines(pd.bias$stc.bias, col="green", lty=5, type="s", lwd=1.5)
lines(pd.bias$stc.ii.bias, col="pink", lty=1, type="s", lwd=1.5)
lines(pd.bias$bucher.bias, col="blue", lty=4, type="s", lwd=1.5)
legend("bottom",lwd=c(1.5,1.5,1.5,1.5),col=c("red","green","pink","blue"),
       lty=c(2,5,1,4),cex=1,bty="n", c("MAIC", "STC","STC-II","Bucher")) # legend
dev.off()

### Nested loop plot for MAE
pdf("Analysis/Supplementary_Material/mae.pdf", width=18, height=15, pointsize=20)
par(pty="m")
plot(pd.other$maic.mae,type="n",ylim=c(0, 0.8),bty="n",xlab="Scenario", 
     ylab="Mean absolute error (MAE)",las=1, xaxt="n")
lines(pd.other, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd.other, which="r",ymin.refline=0.6, ymax.refline=0.8,cex.ref=0.9)
lines(pd.other$maic.mae, col="red", lty=2,type="s", lwd=1.5)    
lines(pd.other$stc.mae, col="green", lty=5, type="s", lwd=1.5)   
lines(pd.other$stc.ii.mae, col="pink", lty=1, type="s", lwd=1.5)
lines(pd.other$bucher.mae, col="blue", lty=4, type="s", lwd=1.5) 
legend("bottom",lwd=c(1.5,1.5,1.5,1.5),col=c("red","green","pink", "blue"),lty=c(2,5,1,4),
       cex=1,bty="n",c("MAIC", "STC", "STC-II", "Bucher"))
dev.off() 

### Nested loop plot of standardized biases
pdf("Analysis/Supplementary_Material/std_biases.pdf", width=18,height=15,pointsize=20)
par(pty="m")
plot(pd.other$maic.std.bias, ylim=c(-400,400), bty="n", xlab="Scenario", ylab="Standardized percentage bias",
     las=1, xaxt="n",cex.axis=1.15, cex.lab=1.15)
abline(h=0, col="grey") # no bias
lines(pd.other, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd.other, which="r", ymin.refline=220, ymax.refline=400, cex.ref=0.9)
lines(pd.other$maic.std.bias, col="red", lty=2,type="s", lwd=1.5)    
lines(pd.other$stc.std.bias, col="green", lty=5, type="s", lwd=1.5)   
lines(pd.other$stc.ii.std.bias, col="pink", lty=1, type="s", lwd=1.5)  
lines(pd.other$bucher.std.bias, col="blue", lty=4, type="s", lwd=1.5) 
legend(x=76,y=215,lwd=c(1.5,1.5,1.5,1.5), col=c("red","green","pink","blue"), lty=c(2,5,1,4), 
       cex=0.9, bty="n", c("MAIC", "STC", "STC-II", "Bucher"))
dev.off()

### Nested loop plot of confidence interval width
pdf("Analysis/Supplementary_Material/CI_width.pdf", width=18,height=15,pointsize=20)
par(pty="m")
plot(pd.other$maic.uci.mean-pd.other$maic.lci.mean, type="n", ylim=c(0, 2), bty="n", xlab="Scenario",
     ylab="95% confidence interval width", las=1, xaxt="n",cex.axis=1.15, cex.lab=1.15)
lines(pd.other, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd.other, which="r", ymin.refline=1.5, ymax.refline=2, cex.ref=0.9)
lines(pd.other$maic.uci.mean-pd.other$maic.lci.mean, col="red", lty=2,type="s", lwd=1.5)    
lines(pd.other$stc.uci.mean-pd.other$stc.lci.mean, col="green", lty=5, type="s", lwd=1.5) 
lines(pd.other$stc.ii.uci.mean-pd.other$stc.ii.lci.mean, col="pink", lty=1, type="s", lwd=1.5) 
lines(pd.other$bucher.uci.mean-pd.other$bucher.lci.mean, col="blue", lty=4, type="s", lwd=1.5) 
legend("bottom", lwd=c(1.5,1.5,1.5,1.5), col=c("red","green","pink","blue"), lty=c(2,5,1,4), 
       cex=0.9, bty="n", c("MAIC", "STC", "STC-II", "Bucher"))
dev.off()

### Nested loop plot of % replicates worse than Bucher
pdf("Analysis/Supplementary_Material/worse_than_bucher.pdf", width=18,height=15,pointsize=20)
par(pty="m")
plot(pd.other$maic.error.worse.than*100, type="n", ylim=c(0, 100), bty="n", xlab="Scenario",
     ylab="% of point estimates worse than Bucher", las=1, xaxt="n",cex.axis=1.15, cex.lab=1.15)
lines(pd.other, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd.other, which="r", ymin.refline=75, ymax.refline=100, cex.ref=0.9)
lines(pd.other$maic.error.worse.than*100, col="red", lty=2,type="s", lwd=1.5)    
lines(pd.other$stc.error.worse.than*100, col="green", lty=5, type="s", lwd=1.5)  
lines(pd.other$stc.ii.error.worse.than*100, col="pink", lty=1, type="s", lwd=1.5) 
legend(x=5, y=60, lwd=c(1.5,1.5,1.5), col=c("red","green","pink"), lty=c(2,5,1), 
       cex=0.9, bty="n", c("MAIC", "STC", "STC-II"))
dev.off()

### Nested loop plot of MAIC % reduction in effective sample size
pdf("Analysis/Supplementary_Material/maic_percentage_reduction_ess.pdf", width=18,height=15,pointsize=20)
par(pty="m")
plot((pd.other$N_AC-pd.other$maic.aess.mean)*100/pd.other$N_AC, type="n", ylim=c(0,100), 
     bty="n", xlab="Scenario", ylab="% reduction in effective sample size", las=1, xaxt="n",
     cex.axis=1.15, cex.lab=1.15)
lines(pd.other, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd.other, which="r", ymin.refline=75, ymax.refline=100, cex.ref=0.9)
lines((pd.other$N_AC-pd.other$maic.aess.mean)*100/pd.other$N_AC, col="red", lty=2,type="s", lwd=1.5)    
legend("bottom", lwd=c(1.5), col=c("red"), lty=c(2), cex=0.9, bty="n", c("MAIC"))
dev.off()

# Plot of bias converging over a simulation scenario (simulation scenario 1)
maic.rolling.bias.1 <- cumsum(maic.means.list[[1]])/(1:replicates) # moving averages
stc.rolling.bias.1 <- cumsum(stc.means.list[[1]])/(1:replicates)
stc.ii.rolling.bias.1 <- cumsum(stc.ii.means.list[[1]])/(1:replicates)
bucher.rolling.bias.1 <- cumsum(bucher.means.list[[1]])/(1:replicates)

pdf("Analysis/Supplementary_material/rolling_bias_convergence.pdf", width=18, height=15, pointsize=20)

plot(maic.rolling.bias.1, xlab="Simulation number", ylab="Rolling bias", col="red",
     bty="l",ylim=c(-0.2,0.4), lwd=1, lty=1, cex.axis=1.15, cex.lab=1.15)
abline(h=0, col="grey") # no bias
lines(stc.rolling.bias.1, col="green", lwd=1, lty=2)
lines(stc.ii.rolling.bias.1 , col="pink", lwd=1, lty=3)
lines(bucher.rolling.bias.1, col="blue", lwd=1, lty=4)
legend("topright", legend = c("MAIC", "STC", "STC-II", "Bucher"), 
       col = c("red", "green", "pink" , "blue"), bty = "n", 
       pt.cex = 2, cex = 1.2, text.col = "black", horiz = F , inset = c(0.1, 0.1),
       lwd=c(1,1,1,1), lty=c(1,2,3,4))
dev.off()
