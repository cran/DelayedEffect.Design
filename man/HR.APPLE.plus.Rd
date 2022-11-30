
\name{HR.APPLE.plus}
\alias{HR.APPLE.plus}
\title{APPLE+ hazard ratio computation}
\description{
   Perform the post-delay hazard ratio calculation given power and sample size using the close-form
APPLE+ method based on the generalized piecewise weighted log-rank test when the treatment time-lag effect is
present and the lag duration varies heterogeneously from individual to individual
or from study to study, within a certain domain and following a specific pattern.                     
}
\usage{ HR.APPLE.plus(lambda1, tl, tu, N, tao, A, beta, ap=0.5, alpha=0.05) 
}
\arguments{
  \item{lambda1}{Baseline hazard}
  \item{tl}{Lower bound of delayed duration domain}
  \item{tu}{Upper bound of delayed duration domain}  
  \item{N}{Sample size}
  \item{tao}{Total study duration}
  \item{A}{Total enrollment duration}
  \item{beta}{Type II error rate; Power=1-beta}
  \item{ap}{Experimental-control allocation ratio. The default is 0.5.}
  \item{alpha}{Type I error rate (two-sided). The default is 0.05.}
}
\details{
 APPLE+ is an acronym for: \cr
 Analytic Power calculation method based on generalized Piecewise weighted Log-rank tEst
with random treatment time-lag effect.
 See the reference for details of this method.
}
\value{ The hazard ratio
}
\references{
Xu, Z., Park, Y., Zhen, B. & Zhu, B. (2017). 
 Achieving optimal power of logrank test with random treatment time-lag effect. 
 Biometrika. Under review. 

Xu, Z., Zhen, B., Park, Y., & Zhu, B. (2017). Designing therapeutic cancer vaccine trials with
delayed treatment effect. Statistics in medicine, 36(4), 592-605.
}
\author{
Zhenzhen Xu <Zhenzhen.Xu@fda.hhs.gov>, Boguang Zhen<Boguang.Zhen@fda.hhs.gov>,
Yongsoek Park <yongpark@pitt.edu> and Bin Zhu <bin.zhu@nih.gov>
}
\examples{
  lambda1 <- 0.001982
  tl      <- 30
  tu      <- 30*11
  N       <- 200
  tao     <- 365*3
  A       <- 365
  beta    <- 0.2
  HR.APPLE.plus(lambda1, tl, tu, N, tao, A, beta)
}
\keyword{ power } 
\seealso{\code{\link{pow.APPLE.plus}}, \code{\link{N.APPLE.plus}}}
