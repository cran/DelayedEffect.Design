# History: Jul 06 2017 Initial coding

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##                                                                     ##
## Power computation on N from 200:1000 for the following methods:     ##
## APPLE, SEPPLE, APPLE.mis, SEPPLE.mis, logrank_analytic, logrank_sim ##
##                                                                     ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

# Function to return lambda1, t1 and p
get_lambda_t1_p <- function(lambda1, t1, p) {

  p.flag   <- !is.null(p)
  t1.flag  <- !is.null(t1)
  lam.flag <- !is.null(lambda1)
  m        <- p.flag + t1.flag + lam.flag
  if (m < 2) stop("ERROR: at least two of (lambda1, t1, p) must be specified")
  if (p.flag) {
    if (p < 0) stop("ERROR: p < 0")
    if (t1.flag) {
     if (t1 <= 0) stop("ERROR: t1 <= 0")
     lambda1 <- -log(p)/t1
    } else {
      # lambda1 is specified
      if (lambda1 == 0) stop("ERROR: lambda1 = 0")
      t1 <- -log(p)/lambda1
    }
  } else {
    # lambda1 and t1 are specified
    p <- exp(-lambda1*t1)
  }

  list(lambda1=lambda1, t1=t1, p=p)

} # END: get_lambda_t1_p

# Function to compute the p-value using survdiff
getPval_0 <- function(X, evt, trt) {

  logrank <- try(survdiff(formula=Surv(X, evt) ~ trt))
  if (!("try-error" %in% class(logrank))) {
    pval <- 1 - pchisq(logrank$chisq, length(logrank$n) - 1)
  } else {
    pval <- NA
  }

  pval

} # END: getPval_0

# Function to compute the p-value using survdiff.fit 
#getPval_1 <- function(X, evt, trt) {
#
#  y   <- cbind(X, evt)
#  fit <- try(survival:::survdiff.fit(y, trt, rho=0))
#  if (!("try-error" %in% class(fit))) {
#    # This code below is from the survdiff function. We need it for
#    # the chi-squared test statistic
#    if (is.matrix(fit$observed)) {
#      otmp <- apply(fit$observed, 1, sum)
#      etmp <- apply(fit$expected, 1, sum)
#    } else {
#      otmp <- fit$observed
#      etmp <- fit$expected
#    }
#    df  <- (etmp > 0)
#    ndf <- sum(df)
#    if (ndf < 2) {
#      chi <- 0
#    } else {
#      temp2 <- ((otmp - etmp)[df])[-1]
#      vv  <- (fit$var[df, df])[-1, -1, drop = FALSE]
#      chi <- sum(solve(vv, temp2) * temp2)
#    }
#    pval <- 1 - pchisq(chi, df=ndf-1)
#  } else {
#    pval <- NA
#  }
#
#  pval
#
#} # END: getPval_1

# Function to compute the p-value by calling survdiff or survdiff.fit.
# Using surdiff.fit is more efficient
getPval <- function(which, X, evt, trt) {

  pval <- getPval_0(X, evt, trt)
  #if (which) {
  #  pval <- getPval_1(X, evt, trt)
  #} else {
  #  pval <- getPval_0(X, evt, trt) 
  #}

  pval

} # END: getPval

# Function to check some parms
checkParms <- function(which, t1, N, HR, tao, A, ap, alpha, nsim, beta) {

  if (N < 1) stop("ERROR with N")
  if (HR < 0) stop("ERROR with HR")
  if (A < 0) stop("ERROR with A")
  if (tao < A) stop("ERROR with tao")
  if (t1 > tao) stop("ERROR with t1")
  if ((ap <= 0) || (ap >= 1)) stop("ERROR with ap")
  if ((alpha <= 0) || (alpha >= 1)) stop("ERROR with tao")
  if (which == 1) {
    if (nsim < 1) stop("ERROR with nsim")
  } else if (which == 2) {
    if ((beta < 0) || (beta > 1)) stop("ERROR with beta")
  }

  NULL

} # END: checkParms

# Function to check if the survival.fit function can be called.
# If not, then we will use the survdiff function
#checkSurvDiffFit <- function() {
#  err  <- 0
#  # Below is a test example found in the survdiff documentation
#  time <- c(59, 115,  156, 421, 431, 448, 464, 475, 477, 563, 638, 744, 769, 770, 803, 855, 
#            1040, 1106, 1129, 1206, 1227,  268, 329, 353, 365, 377)
#  cens <- c(1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0)
#  trt  <- c(1, 1, 1, 2, 1, 1, 2, 2, 1, 2, 1, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 1, 1, 2, 2, 2)
#  fit  <- try(survival:::survdiff.fit(cbind(time, cens), trt, rep(1, length(time)), 0))
#  if ("try-error" %in% class(fit)) {
#    err <- 1
#  } else {
#    vece <- c(5.233531, 6.766469)
#    veco <- c(7, 5)
#    vec  <- c(2.936196, -2.936196, -2.936196, 2.936196)
#    matv <- matrix(vec, byrow=TRUE, nrow=2, ncol=2)
#    if (max(abs(fit$expected - vece)) > 1e-4) err <- 1
#    if (max(abs(fit$var - matv)) > 1e-4) err <- 1
#    if (max(abs(fit$observed - veco)) > 0) err <- 1
#  }
#  if (err) cat("Check survival package. The survdiff.fit function cannot be called directly\n")
#  err
#} # END: checkSurvDiffFit

pow.APPLE <- function(lambda1, t1, p, N, HR, tao, A, ap=0.5, alpha=0.05) 
{
  tmp     <- get_lambda_t1_p(lambda1, t1, p)
  lambda1 <- tmp$lambda1
  t1      <- tmp$t1
  p       <- tmp$p
  checkParms(0, t1, N, HR, tao, A, ap, alpha, NULL, NULL)

  lambda2 <- HR*lambda1  
  M1   <- exp((lambda2-lambda1)*t1)*(exp(-lambda2*t1)*A-(1/lambda2)*exp(-lambda2*tao)*(exp(lambda2*A)-1))
  M2   <- exp(-lambda1*t1)*A-(1/lambda1)*exp(-lambda1*tao)*(exp(lambda1*A)-1)
  tmp1 <- sqrt(ap*(1-ap))*log(1/HR)*sqrt(N*(M1+M2)/(2*A))
  tmp2 <- qnorm(1-alpha/2)   

  X    <- pnorm(tmp1 - tmp2) + pnorm(-tmp1 - tmp2)
  X

} # END: pow.APPLE

HR.APPLE <- function(lambda1, t1, p, N, tao, A, beta, ap=0.5, alpha=0.05) {

  tol     <- 1e-3
  tmp     <- get_lambda_t1_p(lambda1, t1, p)
  lambda1 <- tmp$lambda1
  t1      <- tmp$t1
  p       <- tmp$p
  checkParms(2, t1, N, 1, tao, A, ap, alpha, NULL, beta)

  M1 <- function(x) {
    tmp1 <- x*lambda1
    exp((tmp1 -lambda1)*t1)*(exp(-tmp1*t1)*A-(1/tmp1)*exp(-tmp1*tao)*(exp(tmp1*A)-1))
  }

  M2    <- exp(-lambda1*t1)*A-(1/lambda1)*exp(-lambda1*tao)*(exp(lambda1*A)-1)
  qtmp  <- qnorm(1-alpha/2)
  tmp2  <- 2*A
  tmp3  <- 1 - beta
  tmpap <- sqrt(ap*(1-ap))

  pow.HR <- function(x) {
    tmp <- tmpap*log(1/x)*sqrt(N*(M1(x)+M2)/tmp2)
    ret <- pnorm(tmp - qtmp) + pnorm(-tmp - qtmp) - tmp3
    ret
  }

  HR <- uniroot(pow.HR, lower=1e-6, upper=1, tol=tol)$root
  HR
}

## 1.2 Given power, HR to resolve N ##
N.APPLE <- function(lambda1, t1, p, HR, tao, A, beta, ap=0.5, alpha=0.05) 
{
  tol     <- 1e-3
  tmp     <- get_lambda_t1_p(lambda1, t1, p)
  lambda1 <- tmp$lambda1
  t1      <- tmp$t1
  p       <- tmp$p
  checkParms(2, t1, 100, HR, tao, A, ap, alpha, NULL, beta)

  lambda2 <- HR*lambda1  
  M1      <- exp((lambda2-lambda1)*t1)*(exp(-lambda2*t1)*A-(1/lambda2)*exp(-lambda2*tao)*(exp(lambda2*A)-1))
  M2      <- exp(-lambda1*t1)*A-(1/lambda1)*exp(-lambda1*tao)*(exp(lambda1*A)-1)
  qtmp    <- qnorm(1-alpha/2)
  tmp1    <- sqrt(ap*(1-ap))*log(1/HR)
  tmp2    <- 2*A
  tmp3    <- 1-beta
  tmp4    <- (M1 + M2)/tmp2

  pow.N <- function(x){
    tmp <- tmp1*sqrt(x*tmp4)
    ret <- pnorm(tmp - qtmp) + pnorm(-tmp - qtmp) - tmp3
  }
  N     <- uniroot(pow.N, lower=0, upper=1e7, tol=tol)$root

  N

} # END: N.APPLE

pow.SEPPLE <- function(lambda1, t1, p, N, HR, tao, A, ap=0.5, alpha=0.05, nsim=10000) 
{
  tmp     <- get_lambda_t1_p(lambda1, t1, p)
  lambda1 <- tmp$lambda1
  t1      <- tmp$t1
  p       <- tmp$p
  checkParms(1, t1, N, HR, tao, A, ap, alpha, nsim, NULL)

  Nupper   <- round(N+qnorm(0.99)*sqrt(N)) 
  a        <- N/A
  tvec     <- c(0, t1)
  ratevec1 <- c(lambda1, HR*lambda1)  
  ratevec2 <- c(lambda1, lambda1)
  p.val    <- rep(NA, nsim)
  #which    <- 1 - checkSurvDiffFit()
  which    <- 0

  for (i in 1:nsim)
  {
    t <- 0
    e <- rep(0, Nupper)
    ######## 1. Simulate a Poisson process with intensity A between 0 and A;
    tmp <- log(runif(Nupper))/a
    for (k in 1:Nupper)
    {
      t <- t - tmp[k]
      if (t > A) break 
      e[k] <- t
    }
    ef <- e[e!=0]
    ns <- length(ef)
    if (!ns) next   

    ######## 2. For each enrolled subject, randomize him/her to the treatment or control with 1:1 ratio
    Z <- rbinom(ns, 1, ap)

    #data=cbind(1:ns, ef, Z)
    #colnames(data)=c("id", "enroll", "trt")
    ######## 3. Simulate the time to event from enrollment from a piecewise exponential distribution
    #data.trt=data[data[,3]==1,]
    #data.ctr=data[data[,3]==0,]
    tmp    <- Z == 1
    ef.trt <- ef[tmp]
    Z.trt  <- Z[tmp]
    tmp    <- !tmp
    ef.ctr <- ef[tmp]
    Z.ctr  <- Z[tmp]
    n1     <- length(Z.trt)
    n0     <- length(Z.ctr)
    if ( !n1 || !n0 ) next

    Tt <- rpexp(n1, rate=ratevec1, t=tvec)
    Tc <- rpexp(n0, rate=ratevec2, t=tvec)

    #set.seed(2*N+1+3*N*(i-1))
    #Tt = rpexp(length(data.trt[,3]), rate=c(lambda1, HR*lambda1), t=c(0, t1))
    #data.trt1=cbind(data.trt, Tt)
    #set.seed(2*N+length(Tt)+1+3*N*(i-1))
    #Tc=rpexp(length(data.ctr[,3]), rate=c(lambda1, lambda1), t=c(0, t1))
    #data.ctr1=cbind(data.ctr, Tc)
    #data1=rbind(data.trt1, data.ctr1)
    #colnames(data1)=c("id", "enroll", "trt", "T")
    ######## 4. Calculate the subset of subjects who experience events from t1 to tao during the Phase II 
    #data1.sub=data1[data1[,4]>t1,]
    #data2.sub=cbind(data1.sub[,4]-t1, tao-t1-data1.sub[,2])
    #data3.sub=cbind(data1.sub, data2.sub, apply(data2.sub, 1, min), (data2.sub[,1]<=data2.sub[,2]))
    #colnames(data3.sub)=c("id", "enroll", "trt", "T", "Ttrunc", "tao-t1-enroll", "X", "evt")
    #data3subf=data.frame(data3.sub)
    #logrank=survdiff(formula = Surv(data3subf$X, data3subf$evt) ~ data3subf$trt, data=data3subf)
 
    T        <- c(Tt, Tc)
    eff      <- c(ef.trt, ef.ctr)
    Z        <- c(Z.trt, Z.ctr)
    tmp      <- T > t1 
    data2.1  <- T[tmp] - t1
    data2.2  <- tao - t1 - eff[tmp]
    trt      <- Z[tmp]
    X        <- pmin(data2.1, data2.2)
    evt      <- as.numeric(data2.1 <= data2.2)
    p.val[i] <- getPval(which, X, evt, trt)
  }
  p.val <- p.val[is.finite(p.val)]
  m     <- length(p.val)
  if (!m) {
    stop("ERROR: p-value could not be computed")
  } else if (m < nsim) {
    warning(paste("WARNING: only ", m, " simulations used in computing the p-value", sep=""))
  }
  X <- mean(p.val <= alpha)

  X

} #END: pow.SEPPLE

pow.sim.logrk <- function(lambda1, t1, p, N, HR, tao, A, ap=0.5, alpha=0.05, nsim=10000) 
{
  tmp     <- get_lambda_t1_p(lambda1, t1, p)
  lambda1 <- tmp$lambda1
  t1      <- tmp$t1
  p       <- tmp$p
  checkParms(1, t1, N, HR, tao, A, ap, alpha, nsim, NULL)

  Nupper   <- round(N+qnorm(0.99)*sqrt(N)) 
  a        <- N/A
  tvec     <- c(0, t1)
  ratevec1 <- c(lambda1, HR*lambda1)  
  ratevec2 <- c(lambda1, lambda1)
  p.val    <- rep(NA, nsim)
  #which    <- 1 - checkSurvDiffFit()
  which    <- 0  

  for (i in 1:nsim)
  {
    t <- 0
    e <- rep(0, Nupper)
    ######## 1. Simulate a Poisson process with intensity A between 0 and A;
    tmp <- log(runif(Nupper))/a
    for (k in 1:Nupper)
    {
      t <- t - tmp[k]
      if (t > A) break 
      e[k] <- t
    }
    ef <- e[e!=0]
    ns <- length(ef)
    if (!ns) next
      
    ######## 2. For each enrolled subject, randomize him/her to the treatment or control with 1:1 ratio
    Z <- rbinom(ns, 1, ap)

    #data=cbind(1:ns, ef, Z)
    #colnames(data)=c("id", "enroll", "trt")
    ######## 3. Simulate the time to event from enrollment from a piecewise exponential distribution
    #data.trt=data[data[,3]==1,]
    #data.ctr=data[data[,3]==0,]

    tmp    <- Z == 1
    ef.trt <- ef[tmp]
    Z.trt  <- Z[tmp]
    tmp    <- !tmp
    ef.ctr <- ef[tmp]
    Z.ctr  <- Z[tmp] 
    n1     <- length(Z.trt)
    n0     <- length(Z.ctr)
    if ( !n1 || !n0 ) next

    Tt <- rpexp(n1, rate=ratevec1, t=tvec)
    Tc <- rpexp(n0, rate=ratevec2, t=tvec)

    #Tt = rpexp(length(data.trt[,3]), rate=c(lambda1, HR*lambda1), t=c(0, t1))
    #data.trt1=cbind(data.trt, Tt)
    #Tc=rpexp(length(data.ctr[,3]), rate=c(lambda1, lambda1), t=c(0, t1))
    #data.ctr1=cbind(data.ctr, Tc)
    #data1=rbind(data.trt1, data.ctr1)
    #colnames(data1)=c("id", "enroll", "trt", "T")
    ######## 6. Apply log-rank test to the whole set
    #data2=cbind(data1[,4], tao-data1[,2])
    #data3=cbind(data1, data2[,2], apply(data2, 1, min), (data2[,1]<=data2[,2]))
    #colnames(data3)=c("id", "enroll", "trt", "T", "tao-enroll", "X", "evt")
    #data3f=data.frame(data3)

    T        <- c(Tt, Tc)
    eff      <- c(ef.trt, ef.ctr)
    Z        <- c(Z.trt, Z.ctr)
    data2.1  <- T
    data2.2  <- tao - eff
    trt      <- Z
    X        <- pmin(data2.1, data2.2)
    evt      <- as.numeric(data2.1 <= data2.2)
    p.val[i] <- getPval(which, X, evt, trt)
  }

  p.val <- p.val[is.finite(p.val)]
  m     <- length(p.val)
  if (!m) {
    stop("ERROR: p-value could not be computed")
  } else if (m < nsim) {
    warning(paste("WARNING: only ", m, " simulations used in computing the p-value", sep=""))
  }
  X <- mean(p.val <= alpha)

  X  

} #END: pow.sim.logrk

