
\name{pow.APPLE}
\alias{pow.APPLE}
\title{APPLE power computation}
\description{
   Perform the power calculation using the close-form APPLE method based on the piecewise
weighted log-rank test when the treatment time-lag effect is present and the lag duration 
is homogeneous across the individual subject                           
}
\usage{ pow.APPLE(lambda1, t1, p, N, HR, tao, A, ap=0.5, alpha=0.05)
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
}
\details{
 APPLE is an acronym for: \cr
 Analytic Power calculation method based on Piecewise weighted Log-rank tEst.
 See the reference for details of this method. 

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
\author{
Zhenzhen Xu <Zhenzhen.Xu@fda.hhs.gov>, Boguang Zhen<Boguang.Zhen@fda.hhs.gov>,
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
  pow.APPLE(lambda1, t1, p, N, HR, tao, A)
}
\keyword{ power } 
\seealso{\code{\link{N.APPLE}}, \code{\link{HR.APPLE}}, \code{\link{pow.SEPPLE}}, \code{\link{pow.sim.logrk}}}
