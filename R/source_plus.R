
# Function to check some parms
checkParms <- function(N, lambda1, HR, tl, tu, tao, A, ap, alpha, nsim, beta, dist, shape1, shape2) {

  if (N < 1) stop("ERROR with N")
  if (HR < 0) stop("ERROR with HR")
  if (A < 0) stop("ERROR with A")
  if (tao < A) stop("ERROR with tao")
  if (tl > tu) stop("ERROR with tl/tu")
  if ((ap <= 0) || (ap >= 1)) stop("ERROR with ap")
  if ((alpha <= 0) || (alpha >= 1)) stop("ERROR with tao")
  if (nsim < 1) stop("ERROR with nsim")
  if ((beta < 0) || (beta > 1)) stop("ERROR with beta")
  if (!(dist %in% c("uniform", "beta", "gamma"))) stop("ERROR: dist must be uniform, beta or gamma")
  if (!(dist == "uniform")) {
    if ((is.null(shape1)) || (shape1 <= 0)) stop("ERROR with shape1")
    if ((is.null(shape2)) || (shape2 <= 0)) stop("ERROR with shape2")
  }

  NULL

} # END: checkParms

# Function to get rpexp values
my_rpexp <- function(n, rate, t.star, use.c.code=1) {

  ret <- double(n)

  #if (use.c.code) {
    tmp <- .C("myrpexp", as.integer(n), as.numeric(rate), as.numeric(t.star), ret=ret, PACKAGE="DelayedEffect.Design")
    ret <- tmp$ret
  #} else {
  #  for (j in 1:n) ret[j] <- rpexp(1, rate=rate, t=c(0, t.star[j]))
  #}

  ret

} # END: my_rpexp

# Function to return t.star
get_t.star <- function(n.trt, tl, tu, dchar, shape1, shape2) {

  if (dchar == "u") {
    t.star <- tl+(tu-tl)*runif(n.trt)
  } else if (dchar == "b") {
    t.star <- tl+(tu-tl)*rbeta(n.trt, shape1, shape2)
  } else {
    t.star <- rgamma(n.trt, shape1, shape2)
  }
  tmp <- t.star < 1e-12
  if (any(tmp)) t.star[tmp] <- 1e-12

  t.star

} # END: get_t.star

get_ef <- function(Nupper, A, a, use.c.code=1) {

  vec <- log(runif(Nupper))/a
  e   <- double(Nupper)

  #if (use.c.code) {
    tmp <- .C("getef", as.integer(Nupper), as.numeric(vec), as.numeric(A), ret=e, PACKAGE="DelayedEffect.Design")
    e   <- tmp$ret
  #} else {
  #  t   <- 0
  #  for (k in 1:Nupper) {
  #    t <- t - vec[k]
  #    if (t>A) break 
  #    e[k] <- t
  #  }
  #}
  ef <- e[e!=0]

  ef

} # END: get_ef

# Function to call surdiff 
call_survdiff <- function(X, evt, trt, use.c.code=1) {

  #if (use.c.code) {
    rho    <- 0
    n      <- length(trt)
    ngroup <- 2
    trt    <- match(trt, unique(trt))

    # Define the strat vector for the c function
    nstrat   <- 1
    strat    <- rep(0, n)
    strat[n] <- 1
    ord      <- order(X, -evt)
   
    xx <- .C("logrnk", as.integer(n), as.integer(ngroup), as.integer(nstrat), 
              as.double(rho), as.double(X[ord]), as.integer(evt[ord]), 
              as.integer(trt[ord]), as.integer(strat), 
              observed=double(ngroup), expected=double(ngroup), 
              var.e=double(ngroup*ngroup), double(ngroup), double(n), PACKAGE="DelayedEffect.Design")
  
    ret <- list(expected=xx$expected, observed=xx$observed, 
                var=matrix(xx$var.e, ngroup, ngroup))
    
  #} else {
  #  ret <- survdiff(formula=Surv(X, evt) ~ trt)
  #}

  ret

} # END call_survdiff

# Function to compute logrank p-value 
get_logrankP <- function(X, evt, trt, use.c.code=1) {

    fit <- call_survdiff(X, evt, trt, use.c.code=use.c.code)

    #if (use.c.code) {
      # This code below is from the survdiff function. We need it for
      # the chi-squared test statistic
      if (is.matrix(fit$observed)) {
        otmp <- apply(fit$observed, 1, sum)
        etmp <- apply(fit$expected, 1, sum)
      } else {
        otmp <- fit$observed
        etmp <- fit$expected
      }
      df  <- (etmp > 0)
      ndf <- sum(df)
      if (ndf < 2) {
        chi <- 0
      } else {
        temp2 <- ((otmp - etmp)[df])[-1]
        vv  <- (fit$var[df, df])[-1, -1, drop = FALSE]
        chi <- sum(solve(vv, temp2) * temp2)
      }
      pval <- 1 - pchisq(chi, df=ndf-1)
    #} else {
    #  pval <- 1 - pchisq(fit$chisq, length(fit$n) - 1)
    #}
  
  pval

} # END: get_logrankP

## 1. Power of APPLE+ ##
pow.APPLE.plus <- function(lambda1, tl, tu, N, HR, tao, A, ap=0.5, alpha=0.05) 
{
  checkParms(N, lambda1, HR, tl, tu, tao, A, ap, 0.5, 1000, 0.5, "uniform", 0.5, 0.5)

  E <- exp(-lambda1*tl)*(2/(((HR*lambda1)^2)*((tu-tl)^3)*(lambda1-HR*lambda1))-2*HR*lambda1/((lambda1^3)*((tu-tl)^3)*(lambda1-HR*lambda1))+2*(1-3/(lambda1*(tu-tl)))/((lambda1*(tu-tl))^2))+exp(-lambda1*tu)*(-1/((tu-tl)*(lambda1-HR*lambda1))+HR*lambda1*(1+2/(lambda1*(tu-tl))+2/((lambda1*(tu-tl))^2))/(lambda1*(tu-tl)*(lambda1-HR*lambda1))+(1+4/(lambda1*(tu-tl))+6/((lambda1*(tu-tl))^2))/(lambda1*(tu-tl)))+exp(-lambda1*tl)*exp(-HR*lambda1*(tu-tl))*(-2/(((tu-tl)^2)*(lambda1-HR*lambda1)*(HR*lambda1))-2/(((HR*lambda1)^2)*((tu-tl)^3)*(lambda1-HR*lambda1)))
  D <- 2*(exp(-lambda1*tl)-exp(-lambda1*tu))/((lambda1*(tu-tl))^2)-2*exp(-lambda1*tu)/(lambda1*(tu-tl))
  C <- (HR*lambda1)*(exp(-(lambda1-HR*lambda1)*tl)-exp(-(lambda1-HR*lambda1)*tu))/((tu-tl)*(lambda1-HR*lambda1))
  d_bar <- N*(E+D)/2-N*exp(-lambda1*tao)*(exp(lambda1*A)-1)/(2*lambda1*A)-N*C*(exp(-HR*lambda1*(tao-A))-exp(-HR*lambda1*tao))/(2*((HR*lambda1)^2)*A)

  tmp1 <- sqrt(ap*(1-ap))*log(1/HR)*sqrt(d_bar)
  tmp2 <- qnorm(1-alpha/2)

  pow  <- pnorm( tmp1 - tmp2) + pnorm(-tmp1 - tmp2)
  pow
}

## 1.1 Given power, N to resolve HR ##
HR.APPLE.plus <- function(lambda1, tl, tu, N, tao, A, beta, ap=0.5, alpha=0.05)
{
  checkParms(N, lambda1, 1, tl, tu, tao, A, ap, alpha, 1000, beta, "uniform", 0.5, 0.5)

  E <- function(x) {exp(-lambda1*tl)*(2/(((x*lambda1)^2)*((tu-tl)^3)*(lambda1-x*lambda1))-2*x*lambda1/((lambda1^3)*((tu-tl)^3)*(lambda1-x*lambda1))+2*(1-3/(lambda1*(tu-tl)))/((lambda1*(tu-tl))^2))+exp(-lambda1*tu)*(-1/((tu-tl)*(lambda1-x*lambda1))+x*lambda1*(1+2/(lambda1*(tu-tl))+2/((lambda1*(tu-tl))^2))/(lambda1*(tu-tl)*(lambda1-x*lambda1))+(1+4/(lambda1*(tu-tl))+6/((lambda1*(tu-tl))^2))/(lambda1*(tu-tl)))+exp(-lambda1*tl)*exp(-x*lambda1*(tu-tl))*(-2/(((tu-tl)^2)*(lambda1-x*lambda1)*(x*lambda1))-2/(((x*lambda1)^2)*((tu-tl)^3)*(lambda1-x*lambda1)))}
  D <- 2*(exp(-lambda1*tl)-exp(-lambda1*tu))/((lambda1*(tu-tl))^2)-2*exp(-lambda1*tu)/(lambda1*(tu-tl))
  C <- function(x) {(x*lambda1)*(exp(-(lambda1-x*lambda1)*tl)-exp(-(lambda1-x*lambda1)*tu))/((tu-tl)*(lambda1-x*lambda1))}
  d_bar <- function(x) {N*(E(x)+D)/2-N*exp(-lambda1*tao)*(exp(lambda1*A)-1)/(2*lambda1*A)-N*C(x)*(exp(-x*lambda1*(tao-A))-exp(-x*lambda1*tao))/(2*((x*lambda1)^2)*A)}

  tmp0 <- sqrt(ap*(1-ap))
  tmp2 <- qnorm(1-alpha/2)
  tmp3 <- 1 - beta


  pow.HR <- function(x) {
    tmp1 <- tmp0*log(1/x)*sqrt(d_bar(x))

    pnorm( tmp1 - tmp2) + pnorm(-tmp1 - tmp2)  - tmp3
  }

  HR=uniroot(pow.HR, lower=0.0001, upper=0.9999, tol = 0.01)$root
  HR
}

## 1.2 Given power, HR to resolve N ##
N.APPLE.plus <- function(lambda1, tl, tu, HR, tao, A, beta, ap=0.5, alpha=0.05) 
{
  checkParms(1000, lambda1, HR, tl, tu, tao, A, ap, alpha, 1000, 0.5, "uniform", 0.5, 0.5)

  E <- exp(-lambda1*tl)*(2/(((HR*lambda1)^2)*((tu-tl)^3)*(lambda1-HR*lambda1))-2*HR*lambda1/((lambda1^3)*((tu-tl)^3)*(lambda1-HR*lambda1))+2*(1-3/(lambda1*(tu-tl)))/((lambda1*(tu-tl))^2))+exp(-lambda1*tu)*(-1/((tu-tl)*(lambda1-HR*lambda1))+HR*lambda1*(1+2/(lambda1*(tu-tl))+2/((lambda1*(tu-tl))^2))/(lambda1*(tu-tl)*(lambda1-HR*lambda1))+(1+4/(lambda1*(tu-tl))+6/((lambda1*(tu-tl))^2))/(lambda1*(tu-tl)))+exp(-lambda1*tl)*exp(-HR*lambda1*(tu-tl))*(-2/(((tu-tl)^2)*(lambda1-HR*lambda1)*(HR*lambda1))-2/(((HR*lambda1)^2)*((tu-tl)^3)*(lambda1-HR*lambda1)))
  D <- 2*(exp(-lambda1*tl)-exp(-lambda1*tu))/((lambda1*(tu-tl))^2)-2*exp(-lambda1*tu)/(lambda1*(tu-tl))
  C <- (HR*lambda1)*(exp(-(lambda1-HR*lambda1)*tl)-exp(-(lambda1-HR*lambda1)*tu))/((tu-tl)*(lambda1-HR*lambda1))
  d_bar <- function(x) {x*(E+D)/2-x*exp(-lambda1*tao)*(exp(lambda1*A)-1)/(2*lambda1*A)-x*C*(exp(-HR*lambda1*(tao-A))-exp(-HR*lambda1*tao))/(2*((HR*lambda1)^2)*A)}

  tmp0 <- sqrt(ap*(1-ap))*log(1/HR)
  tmp2 <- qnorm(1-alpha/2)
  tmp3 <- 1 - beta

  pow.N <- function(x) {
    tmp1 <- tmp0*sqrt(d_bar(x))

    pnorm( tmp1 - tmp2) + pnorm(-tmp1 - tmp2) - tmp3
  }

  N=uniroot(pow.N, lower=0, upper=10000000, tol = 0.01)$root
  N
}

## 2. Power of SEPPLE+ ##
pow.SEPPLE.plus <- function(lambda1, tl, tu, N, HR, tao, A, dist="uniform",
                   shape1=NULL, shape2=NULL, ap=0.5, alpha=0.05, nsim=10000) 
{
  checkParms(N, lambda1, HR, tl, tu, tao, A, ap, alpha, nsim, 0.5, dist, shape1, shape2)

  Nupper <- round(N+qnorm(0.99)*sqrt(N)) 
  dchar  <- tolower(substr(dist, 1, 1))
  rate1  <- c(lambda1, HR*lambda1)
  p.val  <- rep(NA, nsim)
  a      <- N/A

  for (i in 1:nsim)
  {
    ef        <- get_ef(Nupper, A, a)
    ns        <- length(ef)
    if (!ns) next
    Z         <- rbinom(ns, 1, ap)
    tmp       <- as.logical(Z)
    Z.trt     <- Z[tmp]
    ef.trt    <- ef[tmp]
    tmp       <- !tmp
    Z.ctr     <- Z[tmp]
    ef.ctr    <- ef[tmp]
    n.trt     <- length(Z.trt)
    n.ctr     <- length(Z.ctr)
    t.star    <- get_t.star(n.trt, tl, tu, dchar, shape1, shape2)
    Tt        <- my_rpexp(n.trt, rate1, t.star)
    Tc        <- rexp(n.ctr, rate=lambda1)
    T         <- c(Tt, Tc)
    eff       <- c(ef.trt, ef.ctr)
    Z         <- c(Z.trt, Z.ctr)
    data2.1   <- T
    data2.2   <- tao - eff
    data3.X   <- pmin(data2.1, data2.2)
    data3.evt <- as.numeric(data2.1 <= data2.2)
    data3.o   <- order(data3.X)
    Z         <- Z[data3.o]
    T         <- T[data3.o]
    data3.evt <- data3.evt[data3.o]
    D1        <- as.numeric((data3.evt == 1) & (Z == 1))
    N         <- length(Z)
    Y1        <- rep(n.trt, N)
    tmp       <- Z == 0
    for (j in 2:N) {
      if (tmp[j-1]) {
        Y1[j] <- Y1[j-1]
      } else {
        Y1[j] <- Y1[j-1] - 1
      }
    }
    w       <- rep(1, N)
    tmp1    <- T < tl
    tmp2    <- (T >= tl) & (T <= tu)
    tmp     <- !tmp1 & tmp2
    w[tmp1] <- 0 
    w[tmp]  <- (T[tmp]-tl)/(tu-tl)
    tmp     <- data3.evt == 1
    Y       <- N:1
    vec     <- w*(D1 - Y1/Y)
    UL      <- sum(vec[tmp])
    vec     <- w*w*Y1*(Y - Y1)/(Y*Y)
    VL      <- sum(vec[tmp])
     
    p.val[i] = pchisq(UL^2/VL, 1, lower.tail=FALSE)
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
}

GPW.logrank <- function(data, obs.time, time.to.event, event.status, trt.group, tl, tu) {
      
 data <- check_dataAndVars(data, obs.time, time.to.event, event.status, trt.group)
 check_tl_tu(tl, tu) 
 ret <- GPW.logrank.main(data[, obs.time, drop=TRUE], data[, time.to.event, drop=TRUE], 
                         data[, event.status, drop=TRUE], data[, trt.group, drop=TRUE], 
                         tl, tu)
 ret
}

GPW.logrank.main <- function(X, T, evt, Z, tl, tu) {

  # X     observational time
  # T     time to event
  # evt   event status
  # Z     treatment group
  
  # Order by event time
  tmp   <- order(X)
  T     <- T[tmp]
  evt   <- evt[tmp]
  Z     <- Z[tmp]
  evt1  <- evt == 1
  Z1    <- Z == 1
  n.trt <- sum(Z1)

  D1    <- as.numeric(evt1 & Z1)
  N     <- length(Z)
  Y1    <- rep(n.trt, N)
  tmp   <- Z == 0
  for (j in 2:N) {
    if (tmp[j-1]) {
      Y1[j] <- Y1[j-1]
    } else {
      Y1[j] <- Y1[j-1] - 1
    }
  }
  w       <- rep(1, N)
  tmp1    <- T < tl
  tmp2    <- (T >= tl) & (T <= tu)
  tmp     <- !tmp1 & tmp2
  w[tmp1] <- 0 
  w[tmp]  <- (T[tmp]-tl)/(tu-tl)
  tmp     <- evt1
  Y       <- N:1
  vec     <- w*(D1 - Y1/Y)
  UL      <- sum(vec[tmp])
  vec     <- w*w*Y1*(Y - Y1)/(Y*Y)
  VL      <- sum(vec[tmp])
     
  p.val   <- pchisq(UL^2/VL, 1, lower.tail=FALSE)
  p.val
}

## 3. Power of SEPPLE on date simulated under the random delayed effect scenario ##
pow.SEPPLE.random.DE <- function(lambda1, tl, tu, N, HR, tao, A, t.fixed, 
           dist="uniform", shape1=NULL, shape2=NULL, ap=0.5, alpha=0.05, nsim=10000) 
{
  checkParms(N, lambda1, HR, tl, tu, tao, A, ap, alpha, nsim, 0.5, dist, shape1, shape2)

  dchar  <- tolower(substr(dist, 1, 1))
  rate1  <- c(lambda1, HR*lambda1)
  p.val  <- rep(NA, nsim)
  Nupper <- round(N+qnorm(0.99)*sqrt(N)) 
  a      <- N/A  

  for (i in 1:nsim)
  {
    ef        <- get_ef(Nupper, A, a)
    ns        <- length(ef) 
    if (!ns) next
    Z         <- rbinom(ns, 1, ap)
    tmp       <- as.logical(Z)
    Z.trt     <- Z[tmp]
    ef.trt    <- ef[tmp]
    tmp       <- !tmp
    Z.ctr     <- Z[tmp]
    ef.ctr    <- ef[tmp]
    n.trt     <- length(Z.trt)
    n.ctr     <- length(Z.ctr)
    if ((!n.trt) || (!n.ctr)) next
    t.star    <- get_t.star(n.trt, tl, tu, dchar, shape1, shape2)
    Tt        <- my_rpexp(n.trt, rate1, t.star)
    Tc        <- rexp(n.ctr, rate=lambda1)
    T         <- c(Tt, Tc)
    eff       <- c(ef.trt, ef.ctr)
    Z         <- c(Z.trt, Z.ctr)
    t.star    <- c(t.star, rep(0, n.ctr))
    data2.1   <- T - t.fixed
    data2.2   <- tao - t.fixed - eff
    data3.X   <- pmin(data2.1, data2.2)
    data3.evt <- as.numeric(data2.1 <= data2.2)
    tmp       <- T > t.fixed
    pval      <- try(get_logrankP(data3.X[tmp], data3.evt[tmp], Z[tmp]), silent=TRUE)
    if (!("try-error" %in% class(pval))) p.val[i] <- pval
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
}

## 4. Power of regular log-rank test on date simulated under the random delayed effect scenario ##
pow.sim.logrk.random.DE <- function(lambda1, tl, tu, N, HR, tao, A, dist="uniform",
                             shape1=NULL, shape2=NULL, ap=0.5, alpha=0.05, nsim=10000) 
{
  checkParms(N, lambda1, HR, tl, tu, tao, A, ap, alpha, nsim, 0.5, dist, shape1, shape2)

  dchar  <- tolower(substr(dist, 1, 1))
  rate1  <- c(lambda1, HR*lambda1)
  p.val  <- rep(NA, nsim)
  Nupper <- round(N+qnorm(0.99)*sqrt(N)) 
  a      <- N/A  

  for (i in 1:nsim)
  {
    ef        <- get_ef(Nupper, A, a)
    ns        <- length(ef)
    if (!ns) next
    Z         <- rbinom(ns, 1, ap)
    tmp       <- as.logical(Z)
    Z.trt     <- Z[tmp]
    ef.trt    <- ef[tmp]
    tmp       <- !tmp
    Z.ctr     <- Z[tmp]
    ef.ctr    <- ef[tmp]
    n.trt     <- length(Z.trt)
    n.ctr     <- length(Z.ctr)
    if ((!n.trt) || (!n.ctr)) next
    t.star    <- get_t.star(n.trt, tl, tu, dchar, shape1, shape2)
    Tt        <- my_rpexp(n.trt, rate1, t.star)
    Tc        <- rexp(n.ctr, rate=lambda1)
    T         <- c(Tt, Tc)
    eff       <- c(ef.trt, ef.ctr)
    Z         <- c(Z.trt, Z.ctr)
    t.star    <- c(t.star, rep(0, n.ctr))
    data2.1   <- T
    data2.2   <- tao - eff
    data3.X   <- pmin(data2.1, data2.2)
    data3.evt <- as.numeric(data2.1 <= data2.2)
    tmp       <- try(get_logrankP(data3.X, data3.evt, Z), silent=TRUE)
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
}
