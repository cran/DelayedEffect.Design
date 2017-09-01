
\name{pow.sim.logrk}
\alias{pow.sim.logrk}
\title{Simulated log-rank power computation}
\description{
  Perform the power calculation using a simulation-based method based on the regular log-rank
test when the treatment time-lag effect is present and the lag duration is homogeneous across the 
 individual subject                          
}
\usage{pow.sim.logrk(lambda1, t1, p, N, HR, tao, A, ap=0.5, alpha=0.05, nsim=10000) 
}
\arguments{
  \item{lambda1}{Baseline hazard or NULL (see details)}
  \item{t1}{Delayed duration or NULL (see details)}
  \item{p}{Proportion of subjects who survive beyond the delayed period or NULL (see details)}
  \item{N}{Sample size}
  \item{HR}{Post-delay hazard ratio}
  \item{tao}{Total study duration}
  \item{A}{Total enrollment duration}
  \item{ap}{Experimental-control allocation ratio. The default is 0.5.}
  \item{alpha}{Type I error rate (two-sided). The default is 0.05.}
  \item{nsim}{Number of simulations. The default is 10000.}
}
\details{
 Out of the three input parameters \code{lambda1}, \code{t1} and \code{p},
 only two need to be specified, the remaining one will be computed
 internally from the formula \code{lambda1 = -log(p)/t1}.
  If all three are not NULL, then 
 \code{lambda1} will be set to \code{-log(p)/t1} regardless of the user input value. 
}
\value{ The power
}
\references{
Xu, Z., Zhen, B., Park, Y., & Zhu, B. (2017). Designing therapeutic cancer vaccine trials with
delayed treatment effect. Statistics in medicine, 36(4), 592-605.
}
\author{Zhenzhen Xu <Zhenzhen.Xu@fda.hhs.gov>, Boguang Zhen<Boguang.Zhen@fda.hhs.gov>,
Yongsoek Park <yongpark@pitt.edu> and Bin Zhu <bin.zhu@nih.gov>
}
\examples{
  lambda1 <- NULL
  t1      <- 183
  p       <- 0.7
  N       <- 200
  HR      <- 0.55
  tao     <- 365*3
  A       <- 365
  pow.sim.logrk(lambda1, t1, p, N, HR, tao, A, nsim=1000)
}
\keyword{ power } 
\seealso{\code{\link{pow.APPLE}}, \code{\link{pow.SEPPLE}}}
