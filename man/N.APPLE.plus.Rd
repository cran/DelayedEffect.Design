
\name{N.APPLE.plus}
\alias{N.APPLE.plus}
\title{APPLE+ sample size computation}
\description{
   Perform the sample size calculation given the power and post-delay hazard ratio using the close-form
APPLE+ method based on the generalized piecewise weighted log-rank test when the treatment time-lag
effect is present and the lag duration varies heterogeneously from individual to individual
or from study to study, within a certain domain and following a specific pattern.                       
}
\usage{ N.APPLE.plus(lambda1, tl, tu, HR, tao, A, beta, ap=0.5, alpha=0.05)
}
\arguments{
  \item{lambda1}{Baseline hazard}
  \item{tl}{Lower bound of delayed duration domain}
  \item{tu}{Upper bound of delayed duration domain}  
  \item{HR}{Post-delay hazard ratio after \code{tu}, defined as the post-delay hazard rate of the treatment group compared to that of the control group}
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
\value{ The sample size
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
  HR      <- 1.3
  tao     <- 365*3
  A       <- 365
  beta    <- 0.2
  N.APPLE.plus(lambda1, tl, tu, HR, tao, A, beta)
}
\keyword{ power } 
\seealso{\code{\link{pow.APPLE.plus}}, \code{\link{HR.APPLE.plus}}}
