# This file contains user-defined MAIC functions, functions to evaluate 
# the performance measures of interest and functions to present simulation 
# results in a nested loop plot (Rücker, G., Schwarzer, G. Presenting simulation 
# results in a nested loop plot. BMC Med Res Methodol 14, 129 (2014) 
# doi:10.1186/1471-2288-14-129)

### MAIC functions
#
# function to estimate MAIC weights (assumes M>1)
maic <- function(A.X, B.summary) {
  M <- length(B.summary)
  for (i in 1:M) {
      A.X[,i] <- A.X[,i] - B.summary[,i] # only means are balanced
  }
  A.X <- as.matrix(A.X)
  # objective function for weight estimation
  objective.function <- function(beta.pars,X){
    return(sum(exp(X %*% beta.pars)))
  }
  # arbitrary starting point for the optimiser
  beta.start<-rep(1,M)
  # optimisation (BFGS)
  out.optim <- optim(fn=objective.function, X=A.X, par=beta.start, method = "BFGS")
  # betas
  beta.pars <- out.optim$par
  # weights
  temp.w<-rep(0,nrow(A.X))
  for (i in 1:M){
    temp.w<-temp.w+beta.pars[i]*A.X[,i]
  }
  w<-exp(temp.w)
  return(w)
}

# approx. effective sample size of weights
approx.ess <- function(w) {
  ess <- sum(w)^2/sum(w^2)
  return(ess)
}

### Functions to evaluate performance measures
#
# bias estimate
bias <- function(theta.hat, theta) {
  nsim <- length(theta.hat)
  est <- sum(theta.hat)/nsim - theta
  return(est)
}

# Monte Carlo SE of bias estimate
bias.mcse <- function(theta.hat) {
  nsim <- length(theta.hat)
  tmp <- sum((theta.hat - mean(theta.hat))^2)
  mcse <- sqrt(1/(nsim*(nsim-1))*tmp)
  return(mcse)
}

# coverage estimate
coverage <- function(theta.hat.low, theta.hat.upp, theta) {
  nsim <- length(theta.hat.low)
  est <- sum(ifelse(theta>=theta.hat.low & theta<=theta.hat.upp,1,0))/nsim
  return(est)
}

# Monte Carlo SE of coverage estimate
coverage.mcse <- function(coverage, nsim) {
  mcse <- sqrt((coverage*(1-coverage))/nsim)
  return(mcse)
}

# MSE estimate
mse <- function(theta.hat, theta) {
  nsim <- length(theta.hat)
  est <- sum((theta.hat-theta)^2)/nsim
  return(est)
}

# Monte Carlo SE of MSE estimate
mse.mcse <- function(theta.hat, theta) {
  nsim <- length(theta.hat)
  tmp <- (theta.hat-theta)^2
  mse.est <- sum(tmp)/nsim
  mcse <- sqrt(sum((tmp - mse.est)^2)/(nsim*(nsim-1)))
  return(mcse)
}

# MAE estimate
mae <- function(theta.hat, theta) {
  nsim <- length(theta.hat)
  est <- sum(abs(theta.hat-theta))/nsim
  return(est)
}

# Monte Carlo SE of any continuous estimate
mcse.estimate <- function(perf.measure) {
  nsim <- length(perf.measure)
  perf.measure.mean <- sum(perf.measure)/nsim
  mcse <- sqrt(sum((perf.measure-perf.measure.mean)^2)/(nsim*(nsim-1)))
  return(mcse)
}

# Empirical standard error 
empse <- function(theta.hat) {
  nsim <- length(theta.hat)
  tmp <- sum((theta.hat - mean(theta.hat))^2)
  est <- sqrt(tmp/(nsim-1))
  return(est)
}

# EmpSE MCSE
empse.mcse <- function(empse, nsim) {
  mcse <- empse/(sqrt(2*(nsim-1)))
  return(mcse)
} 

# Variability ratio
var.ratio <- function(theta.hat, std.err) {
  nsim <- length(theta.hat)
  num <- sum(std.err)/nsim
  denom <- sqrt(sum((theta.hat-mean(theta.hat))^2)/(nsim-1))
  est <- num/denom
  return(est)    
}

# Variability ratio MCSE
var.ratio.mcse <- function(avg.se, emp.se, var.avg.se, var.emp.se) {
  # approximation of ratio variance based on independence of avg. se and emp.se
  # see Wolter, K., 2007. Introduction to variance estimation. 
  mcse <- sqrt((1/emp.se^2)*var.avg.se + (((avg.se^2)/(emp.se^4))*var.emp.se))
  return(mcse)             
}

### functions to present simulation results in a nested loop plot
#
# These are by: Rücker, G., Schwarzer, G. Presenting simulation results in a nested loop plot. 
# BMC Med Res Methodol 14, 129 (2014) doi:10.1186/1471-2288-14-129
#
# function reorders simulation dataset x in order varnames
nestedloop <- function(x,
                       varnames, sign=rep(1, length(varnames)),
                       varlabels=NULL){
  ##
  if (!inherits(x, "data.frame"))
    stop("Argument 'x' must be a data.frame.")
  ##
  mo <- matrix(sign,
               nrow=dim(x)[[1]], ncol=length(varnames),
               byrow=TRUE)
  xo <- x[,varnames]
  ##
  ## Re-ordering:
  res <- x[do.call(order, mo*xo),]
  ##
  attr(res, "varnames")  <- varnames
  attr(res, "varlabels") <- varlabels
  attr(res, "sign") <- sign
  ##
  class(res) <- c("nestedloop", class(res))
  ##
  res
}

# plots lines in nested loop plot 
lines.nestedloop <- function(x,
                             varnames=attributes(x)$varnames,
                             varlabels=attributes(x)$varlabels,
                             which="v",
                             col=if (which=="r") "#999999" else "black",
                             ymin.refline, ymax.refline,
                             cex.ref=0.9,
                             log=TRUE,
                             ...){
  ##
  nvar <- length(varnames)
  ##
  if (length(col)==1)
    col <- rep(col, nvar)
  ##
  if (which=="v"){
    ##
    ## Vertical lines
    ##
    nlen <- rep(NA, nvar)
    ##
    for (i in 1:nvar)
      nlen[i] <- length(unique(x[,varnames[i]]))
    ##
    cnlen <- cumprod(nlen)
    ##
    for (i in (nvar-1):1)
      abline(v=cnlen[nvar]*(0:cnlen[i])/cnlen[i]+1,
             col=col[i])
  }
  else if (which=="r"){
    ##
    ## Reference lines
    ##
    if (is.null(varlabels))
      varlabels <- varnames
    ##
    labels.varnames <- rep("", nvar)
    ##
    for (i in 1:length(varnames)){
      if (is.factor(x[,varnames[i]]))
        varvals <- unique(x[,varnames[i]])
      else{
        varvals <- format(unique(x[,varnames[i]]))
        varvals <- sub("^[[:space:]]*(.*?)[[:space:]]*$",
                       "\\1",
                       varvals,
                       perl=TRUE)
      }
      ##
      labels.varnames[i] <- paste(varlabels[i],
                                  " (",
                                  paste(varvals,
                                        collapse=", "),
                                  ")", sep="")
    }
    ##
    if (log){
      ymax <- log(ymax.refline)
      ymin <- log(ymin.refline)
    }
    else{
      ymax <- ymax.refline
      ymin <- ymin.refline
    }
    ##
    distance <- (ymax-ymin)/nvar
    ##
    ypos <- ymax-0.2*distance-(1/nvar)*(0:(nvar-1))*(ymax-ymin)
    ypos.ref.max <- ypos-0.20*distance
    ypos.ref.min <- ypos-0.75*distance
    ##
    if (log){
      ypos <- exp(ypos)
      ypos.ref.max <- exp(ypos.ref.max)
      ypos.ref.min <- exp(ypos.ref.min)
    }
    ##
    for (i in 1:nvar){
      ##print(c(ypos.ref.max[i], ypos.ref.min[i]))
      ##
      ##abline(h=ypos[i], lwd=1, col="red")
      ##abline(h=ypos.ref.max[i], lwd=1, col="blue")
      ##abline(h=ypos.ref.min[i], lwd=1, col="green")
      text(1, ypos[i], labels.varnames[i], adj=0, cex=cex.ref)
      ##
      xvar <- x[,varnames[i]]
      if (is.factor(xvar))
        xvar <- as.numeric(xvar)
      xvar <- ypos.ref.min[i] +
        (xvar-min(xvar))/(max(xvar)-min(xvar))*
        (ypos.ref.max[i]-ypos.ref.min[i])
      lines(xvar, col=col[i], type="s", lwd=1)
    }
  }
  ##
  invisible(NULL)
}
