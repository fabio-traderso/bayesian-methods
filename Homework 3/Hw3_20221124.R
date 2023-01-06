
library(bayesm)
rm(list = ls())

rbprobitGibbsHW <- function(Data, Prior, Mcmc) {
  breg1 <- function(root, X, z, Abetabar) {
    #
    # Argument 'bounds' is a k x 2 matrix, specifying the lower bounds in the first,
    # upper bounds in the second column
    #
    beta <- betainit
    cov <- crossprod(root, root)
    betatilde <- cov %*% (crossprod(X, z) + Abetabar)
    
    for (i in seq_len(ncol(X))) {
      mu <- betatilde[i] - as.vector(beta[-i] - betatilde[-i]) %*% as.vector(cov[-i, i]) / cov[i, i]
      sigma <- chol(1 / cov[i, i])
      beta[i] <- mu + sigma%*% rnorm(1)
    }
    return(beta)
    
  }
  
  # Input validation
  pandterm <- function(message) {
    stop(message, call. = FALSE)
  }
  #
  # ----------------------------------------------------------------------
  #
  # check arguments
  #
  if (missing(Data)) {
    pandterm("Requires Data argument -- list of y and X")
  }
  if (is.null(Data$X)) {
    pandterm("Requires Data element X")
  }
  X <- Data$X
  if (is.null(Data$y)) {
    pandterm("Requires Data element y")
  }
  y <- Data$y
  nvar <- ncol(X)
  nobs <- length(y)
  #
  # check data for validity
  #
  if (length(y) != nrow(X)) {
    pandterm("y and X not of same row dim")
  }
  if (sum(unique(y) %in% c(0:1)) < length(unique(y))) {
    pandterm("Invalid y, must be 0,1")
  }
  #
  # check for Prior
  #
  if (missing(Prior)) {
    betabar <- c(rep(0, nvar))
    A <- .01 * diag(nvar)
  } else {
    if (is.null(Prior$betabar)) {
      betabar <- c(rep(0, nvar))
    } else {
      betabar <- Prior$betabar
    }
    if (is.null(Prior$A)) {
      A <- .01 * diag(nvar)
    } else {
      A <- Prior$A
    }
  }
  #
  # check dimensions of Priors
  #
  if (ncol(A) != nrow(A) || ncol(A) != nvar || nrow(A) != nvar) {
    pandterm(paste("bad dimensions for A", dim(A)))
  }
  if (length(betabar) != nvar) {
    pandterm(paste("betabar wrong length, length= ", length(betabar)))
  }
  #
  # check MCMC argument
  #
  if (missing(Mcmc)) {
    pandterm("requires Mcmc argument")
  } else {
    if (is.null(Mcmc$R)) {
      pandterm("requires Mcmc element R")
    } else {
      R <- Mcmc$R
    }
    if (is.null(Mcmc$keep)) {
      keep <- 1
    } else {
      keep <- Mcmc$keep
    }
  }
  #
  # print out problem
  #
  cat(" ", fill = TRUE)
  cat("Starting Gibbs Sampler for Binary Probit Model", fill = TRUE)
  cat("   with ", length(y), " observations", fill = TRUE)
  cat("Table of y Values", fill = TRUE)
  print(table(y))
  cat(" ", fill = TRUE)
  cat("Prior Parms: ", fill = TRUE)
  cat("betabar", fill = TRUE)
  print(betabar)
  cat("A", fill = TRUE)
  print(A)
  cat(" ", fill = TRUE)
  cat("MCMC parms: ", fill = TRUE)
  cat("R= ", R, " keep= ", keep, fill = TRUE)
  cat(" ", fill = TRUE)
  
  betadraw=matrix(double(floor(R/keep)*nvar),ncol=nvar)
  beta=c(rep(0,nvar))
  sigma=c(rep(1,nrow(X)))
  root=chol(chol2inv(chol((crossprod(X,X)+A))))
  Abetabar=crossprod(A,betabar)
  a=ifelse(y == 0,-100, 0)
  b=ifelse(y == 0, 0, 100)
  #
  #	start main iteration loop
  #
  itime=proc.time()[3]
  cat("MCMC Iteration (est time to end - min) ",fill=TRUE)
  
  
  for (rep in 1:R) 
  {
    # draw z given beta(i-1)
    mu=X%*%beta
    z=rtrun(mu,sigma,a,b)
    beta=breg1(root,X,z,Abetabar)
    #
    #       print time to completion and draw # every 100th draw
    #
    if(rep%%100 == 0)
    {ctime=proc.time()[3]
    timetoend=((ctime-itime)/rep)*(R-rep)
    cat(" ",rep," (",round(timetoend/60,1),")",fill=TRUE)
    }
    
    if(rep%%keep == 0) 
    {mkeep=rep/keep; betadraw[mkeep,]=beta}
  }
  ctime = proc.time()[3]
  cat('  Total Time Elapsed: ',round((ctime-itime)/60,2),'\n')
  attributes(betadraw)$class=c("bayesm.mat","mcmc")
  attributes(betadraw)$mcpar=c(1,R,keep)
  return(list(betadraw=betadraw))
}

simbprobit = function(X, beta) {
  y = ifelse((X%*%beta + rnorm(nrow(X)))<0, 0, 1)
  list(X=X, y=y, beta=beta)
}

nobs = 100
X = cbind(rep(1,nobs), runif(nobs), runif(nobs))
beta = c(0,1,-1)
nvar = ncol(X)
simout = simbprobit(X, beta)

R=50000
Data1 = list(X=simout$X, y=simout$y)
Mcmc1 = list(R=R, keep=1)
Prior= list(betabar=c(1,2,3), A=matrix(rep(0,9),ncol=3))

out = rbprobitGibbs(Data=Data1, Mcmc=Mcmc1)
summary(out$betadraw, tvalues=beta)
windows()
plot(out$betadraw, tvalues=beta)


# lower und upper truncation points for the three element paramter vector
a=c(-Inf,0,-Inf)
b=c(+Inf,+Inf,0)

