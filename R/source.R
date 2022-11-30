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

# Function to check some parms
checkParms0 <- function(which, t1, N, HR, tao, A, ap, alpha, nsim, beta) {

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

# Function to get rpexp values
my_rpexp0 <- function(n, rate, t1, use.c.code=1) {

  ret <- double(n)

  #if (use.c.code) {
    tmp <- .C("myrpexp0", as.integer(n), as.numeric(rate), as.numeric(t1), ret=ret, PACKAGE="DelayedEffect.Design")
    ret <- tmp$ret
  #} else {
  #  ret <- rpexp(n, rate=rate, t=c(0, t1))
  #}

  ret

} # END: my_rpexp


pow.APPLE <- function(lambda1, t1, p, N, HR, tao, A, ap=0.5, alpha=0.05) 
{
  tmp     <- get_lambda_t1_p(lambda1, t1, p)
  lambda1 <- tmp$lambda1
  t1      <- tmp$t1
  p       <- tmp$p
  checkParms0(0, t1, N, HR, tao, A, ap, alpha, NULL, NULL)

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
  checkParms0(2, t1, N, 1, tao, A, ap, alpha, NULL, beta)

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
  checkParms0(2, t1, 100, HR, tao, A, ap, alpha, NULL, beta)

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
  checkParms0(1, t1, N, HR, tao, A, ap, alpha, nsim, NULL)

  Nupper   <- round(N+qnorm(0.99)*sqrt(N)) 
  a        <- N/A
  tvec     <- c(0, t1)
  ratevec1 <- c(lambda1, HR*lambda1)  
  ratevec2 <- c(lambda1, lambda1)
  p.val    <- rep(NA, nsim)
  
  for (i in 1:nsim)
  {
    ef <- get_ef(Nupper, A, a)
    ns <- length(ef)
    if (!ns) next   

    ######## 2. For each enrolled subject, randomize him/her to the treatment or control with 1:1 ratio
    Z      <- rbinom(ns, 1, ap)
    tmp    <- Z == 1
    ef.trt <- ef[tmp]
    Z.trt  <- Z[tmp]
    tmp    <- !tmp
    ef.ctr <- ef[tmp]
    Z.ctr  <- Z[tmp]
    n1     <- length(Z.trt)
    n0     <- length(Z.ctr)
    if ( !n1 || !n0 ) next

    Tt <- my_rpexp0(n1, ratevec1, t1)
    Tc <- my_rpexp0(n0, ratevec2, t1)
    #Tt <- rpexp(n1, rate=ratevec1, t=tvec)
    #Tc <- rpexp(n0, rate=ratevec2, t=tvec)

    T        <- c(Tt, Tc)
    eff      <- c(ef.trt, ef.ctr)
    Z        <- c(Z.trt, Z.ctr)
    tmp      <- T > t1 
    data2.1  <- T[tmp] - t1
    data2.2  <- tao - t1 - eff[tmp]
    trt      <- Z[tmp]
    X        <- pmin(data2.1, data2.2)
    evt      <- as.numeric(data2.1 <= data2.2)
    tmp      <- try(get_logrankP(X, evt, trt), silent=TRUE)
    if (!("try-error" %in% class(tmp))) p.val[i] <- tmp
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
  checkParms0(1, t1, N, HR, tao, A, ap, alpha, nsim, NULL)

  Nupper   <- round(N+qnorm(0.99)*sqrt(N)) 
  a        <- N/A
  tvec     <- c(0, t1)
  ratevec1 <- c(lambda1, HR*lambda1)  
  ratevec2 <- c(lambda1, lambda1)
  p.val    <- rep(NA, nsim)

  for (i in 1:nsim)
  {
    ef <- get_ef(Nupper, A, a)
    ns <- length(ef)
    if (!ns) next
      
    ######## 2. For each enrolled subject, randomize him/her to the treatment or control with 1:1 ratio
    Z      <- rbinom(ns, 1, ap)
    tmp    <- Z == 1
    ef.trt <- ef[tmp]
    Z.trt  <- Z[tmp]
    tmp    <- !tmp
    ef.ctr <- ef[tmp]
    Z.ctr  <- Z[tmp] 
    n1     <- length(Z.trt)
    n0     <- length(Z.ctr)
    if ( !n1 || !n0 ) next

    Tt <- my_rpexp0(n1, ratevec1, t1)
    Tc <- my_rpexp0(n0, ratevec2, t1)
    #Tt <- rpexp(n1, rate=ratevec1, t=tvec)
    #Tc <- rpexp(n0, rate=ratevec2, t=tvec)

    T        <- c(Tt, Tc)
    eff      <- c(ef.trt, ef.ctr)
    Z        <- c(Z.trt, Z.ctr)
    data2.1  <- T
    data2.2  <- tao - eff
    trt      <- Z
    X        <- pmin(data2.1, data2.2)
    evt      <- as.numeric(data2.1 <= data2.2)
    tmp      <- try(get_logrankP(X, evt, Z), silent=TRUE)
    if (!("try-error" %in% class(tmp))) p.val[i] <- tmp
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

