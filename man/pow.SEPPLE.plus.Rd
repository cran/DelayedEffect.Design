
\name{pow.SEPPLE.plus}
\alias{pow.SEPPLE.plus}
\title{SEPPLE+ power computation}
\description{
  Perform the power calculation using the numeric SEPPLE+ method based on the 
 generalized piecewise
weighted log-rank test when the treatment time-lag effect is present and the lag duration 
varies heterogeneously from individual to individual or from study to study, within a
certain domain and following a specific pattern.                      
}
\usage{pow.SEPPLE.plus(lambda1, tl, tu, N, HR, tao, A, dist="uniform",
              shape1=NULL, shape2=NULL, ap=0.5, alpha=0.05, nsim=10000) 
}
\arguments{
  \item{lambda1}{Baseline hazard}
  \item{tl}{Lower bound of delayed duration domain}
  \item{tu}{Upper bound of delayed duration domain}  
  \item{N}{Sample size}
  \item{HR}{Post-delay hazard ratio after tu, defined as the post-delay hazard rate of the treatment group compared to that of the control group}
  \item{tao}{Total study duration}
  \item{A}{Total enrollment duration}
  \item{dist}{One of "uniform", "beta" or "gamma", for the lag distribution}
  \item{shape1}{NULL or a positive parameter value for the \code{beta} or
                \code{gamma} distribution.}
  \item{shape2}{NULL or a positive parameter value for the \code{beta} or
                \code{gamma} distribution.}
  \item{ap}{Experimental-control allocation ratio. The default is 0.5.}
  \item{alpha}{Type I error rate (two-sided). The default is 0.05.}
  \item{nsim}{Number of simulations. The default is 10000.}
}
\details{
 SEPPLE+ is an acronym for: \cr
 Simulation-based Empirical Power calculation method based on generalized Piecewise weighted Log-rank
   tEst with random treatment time-lag effect. 
 See the reference for details of this method. 
}
\value{ The power
}
\references{
Xu, Z., Park, Y., Zhen, B. & Zhu, B. (2017). 
 Achieving optimal power of logrank test with random treatment time-lag effect. 
 Biometrika. Under review. 

Xu, Z., Zhen, B., Park, Y., & Zhu, B. (2017). Designing therapeutic cancer vaccine trials with
delayed treatment effect. Statistics in medicine, 36(4), 592-605.
}
\author{Zhenzhen Xu <Zhenzhen.Xu@fda.hhs.gov> , Boguang Zhen<Boguang.Zhen@fda.hhs.gov>,
Yongsoek Park <yongpark@pitt.edu> and Bin Zhu <bin.zhu@nih.gov>}
\examples{
  lambda1 <- 0.001982
  tl      <- 30
  tu      <- 30*11
  N       <- 200
  HR      <- 0.55
  tao     <- 365*3
  A       <- 365
  shape1  <- 5
  shape2  <- 5
  pow.SEPPLE.plus(lambda1, tl, tu, N, HR, tao, A, dist="beta", 
                  shape1=shape1, shape2=shape2, nsim=1000)
}
\keyword{ power } 
\seealso{\code{\link{pow.SEPPLE.random.DE}}, \code{\link{pow.sim.logrk.random.DE}}}
