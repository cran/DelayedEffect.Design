
\name{DelayedEffect.Design}
\alias{DelayedEffect.Design}
\docType{package}
\title{
 Sample size and power calculations using APPLE and SEPPLE
}
\description{
 An R package for sample size and power calculation when the treatment time-lag effect is
present and the lag duration is homogeneous across the individual subject using the APPLE and
SEPPLE methods based on the piecewise weighted log-rank test. For comparison, this package
also performs the power calculation based on the regular log-rank test which ignores 
the existence of lag effect.
}
\details{
 The two new methods in this package for performing the sample size and power calculations are: \cr
1. Analytic Power calculation method based on Piecewise weighted Log-rank tEst (APPLE), \cr
2. Simulation-based Empirical Power calculation method based on Piecewise weighted Log-rank
   tEst (SEPPLE). \cr
  See the reference for details of these methods and the piecewise weighted log-rank test.
  The functions for computing power corresponding to the above methods are \code{\link{pow.APPLE}} and
  \code{\link{pow.SEPPLE}}. These can be compared to \code{\link{pow.sim.logrk}},
  which computes the power from a simulation-based algorithm using
  the regular log-rank test which ignores the existence of lag effect.

 This package also includes the function \code{\link{N.APPLE}} 
 to back calculate the sample size given the power and hazard ratio, 
 and the function \code{\link{HR.APPLE}} 
 to back calculate the hazard ratio given
the power and sample size, respectively, using the close-form APPLE method.
}
\references{
  Xu, Z., Zhen, B., Park, Y., & Zhu, B. (2017). Designing therapeutic cancer vaccine trials with
delayed treatment effect. Statistics in medicine, 36(4), 592-605.
}
\author{Zhenzhen Xu <Zhenzhen.Xu@fda.hhs.gov> , Boguang Zhen<Boguang.Zhen@fda.hhs.gov>,
Yongsoek Park <yongpark@pitt.edu> and Bin Zhu <bin.zhu@nih.gov>
}
\keyword{package}