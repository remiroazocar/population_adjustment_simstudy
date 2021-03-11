# Example R code implementing MAIC, STC and the Bucher method on a simulated example 

library("survival") # required for weighted and standard Cox regressions

# setwd("C:/Users/Antonio/Desktop/population_adjustment_simstudy/Example") 

### MAIC ###

rm(list=ls())
AC.IPD <- read.csv("AC_IPD.csv") # load AC patient-level data
BC.ALD <- read.csv("BC_ALD.csv") # load BC aggregate-level data

N <- nrow(AC.IPD) # number of subjects in AC
X.EM <- AC.IPD[,c("X1","X2")] # AC effect modifiers 
bar.X.EM.BC <- BC.ALD[c("mean_X1", "mean_X2")] # BC effect modifier means
K.EM <- ncol(X.EM) # number of effect modifiers 

# center the AC effect modifiers on the BC means
for (k in 1:K.EM) {
  X.EM[,k] <- X.EM[,k] - bar.X.EM.BC[,k]
}

# objective function to be minimized for weight estimation
Q <- function(alpha, X.EM) {
  return(sum(exp(X.EM %*% alpha)))
}

alpha <- rep(1,K.EM) # arbitrary starting point for the optimiser
# objective function minimized using BFGS
Q.min <- optim(fn=Q, X.EM=as.matrix(X.EM), par=alpha, method="BFGS")
hat.alpha <- Q.min$par # finite solution is the logistic regression parameters
log.hat.w <- rep(0, N)
for (k in 1:K.EM) {
  log.hat.w <- log.hat.w + hat.alpha[k]*X.EM[,k]
}
hat.w <- exp(log.hat.w) # estimated weights 
aess <- sum(hat.w)^2/sum(hat.w^2) # approximate effective sample size

# fit weighted Cox proportional hazards model using robust=TRUE for robust variance
outcome.fit <- coxph(Surv(time, status)~trt, robust=TRUE, weights=hat.w, data=AC.IPD)

# fitted treatment coefficient is relative effect for A vs. C
hat.Delta.AC <- summary(outcome.fit)$coef[1] 
hat.var.Delta.AC <- vcov(outcome.fit)[[1]] # estimated variance for A vs. C
hat.Delta.BC <- with(BC.ALD, logHR_B) # B vs. C
hat.var.Delta.BC <- with(BC.ALD, var_logHR_B)
hat.Delta.AB <- hat.Delta.AC - hat.Delta.BC # A vs. B
hat.var.Delta.AB <- hat.var.Delta.AC + hat.var.Delta.BC

### STC ###

rm(list=ls())
AC.IPD <- read.csv("AC_IPD.csv") # load AC patient-level data
BC.ALD <- read.csv("BC_ALD.csv") # load BC aggregate-level data

# fit regression of outcome on the baseline characteristics and treatment
# effect modifiers are centered at the mean BC values
# purely prognostic variables are included but not centered
outcome.fit <- coxph(Surv(time, status)~X3+X4+trt*I(X1-BC.ALD$mean_X1)+trt*I(X2-BC.ALD$mean_X2),
                     data=AC.IPD)

# estimated treatment coefficient is relative effect for A vs. C
hat.Delta.AC <- coef(outcome.fit)["trt"]
hat.var.Delta.AC <- vcov(outcome.fit)["trt", "trt"] # estimated variance for A vs. C
hat.Delta.BC <- with(BC.ALD, logHR_B) # B vs. C
hat.var.Delta.BC <- with(BC.ALD, var_logHR_B)
hat.Delta.AB <- hat.Delta.AC - hat.Delta.BC # A vs. B
hat.var.Delta.AB <- hat.var.Delta.AC + hat.var.Delta.BC

### Bucher method ###

rm(list=ls())
AC.IPD <- read.csv("AC_IPD.csv") # load AC patient-level data
BC.ALD <- read.csv("BC_ALD.csv") # load BC aggregate-level data

# simple regression of outcome on treatment
outcome.fit <- coxph(Surv(time, status)~trt, data=AC.IPD)

# fitted treatment coefficient is relative effect for A vs. C
hat.Delta.AC <- coef(outcome.fit)["trt"]
hat.var.Delta.AC <- vcov(outcome.fit)["trt", "trt"] # estimated variance for A vs. C
hat.Delta.BC <- with(BC.ALD, logHR_B) # B vs. C
hat.var.Delta.BC <- with(BC.ALD, var_logHR_B)  
hat.Delta.AB <- hat.Delta.AC - hat.Delta.BC # A vs. B
hat.var.Delta.AB <- hat.var.Delta.AC + hat.var.Delta.BC
