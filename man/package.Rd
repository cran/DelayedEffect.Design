
\name{DelayedEffect.Design}
\alias{DelayedEffect.Design}
\docType{package}
\title{
 Sample size and power calculations using the APPLE, SEPPLE, APPLE+ and SEPPLE+ methods
}
\description{
 An R package for sample size and power calculation when the treatment time-lag 
effect is present. The package incorporates two specific lag assumptions: \cr
1. the lag duration is homogeneous across the individual subject; \cr
2. the lag duration varies heterogeneously from individual to individual within a 
   certain domain and following a specific pattern.
}
\details{
The four new methods in this package for performing the sample size and power calculations are:\cr
1. Analytic Power calculation method based on Piecewise weighted Log-rank tEst (APPLE), \cr
2. Simulation-based Empirical Power calculation method based on Piecewise weighted Log-rank tEst (SEPPLE), \cr
3.  Analytic Power calculation method based on generalized Piecewise weighted Log-rank tEst
with random treatment time-lag effect (APPLE+), \cr
4. Simulation-based Empirical Power calculation method based on generalized Piecewise weighted 
Log-rank tEst with random treatment time-lag effect (SEPPLE+). \cr
See the reference for details of these methods.
 Specifically, APPLE and SEPPLE assume that the lag duration is homogeneous across the individual subject,
 whereas APPLE and SEPPLE assume that the lag duration varies heterogeneously from individual to individual 
or from study to study within a certain domain and following a specific pattern.  
The functions for computing power corresponding to the above methods are pow.APPLE, pow.SEPPLE, pow.APPLE.plus,
 pow.SEPPLE.plus and pow.SEPPLE.random.DE. These can be compared to pow.sim.logrk and pow.sim.logrk.rankdom.DE 
which compute the power from a simulation-based algorithm using the regular log-rank test which ignores the 
existence of lag effects. 
The package also includes the function N.APPLE, N.APPLE.plus to back calculate the sample size given the power 
and hazard ratio, and the functions HR.APPLE and HR.APPLE.plus to back calculate the hazard ratio given the 
power and sample size, respectively, using the close-from APPLE and APPLE+ methods.
}
\references{
 Xu, Z., Park, Y., Zhen, B. & Zhu, B. (2017). 
 Achieving optimal power of logrank test with random treatment time-lag effect. 
 Biometrika. Under review. 

 Xu, Z., Zhen, B., Park, Y., & Zhu, B. (2017). Designing therapeutic cancer vaccine trials with
delayed treatment effect. Statistics in medicine, 36(4), 592-605.
}
\author{Zhenzhen Xu <Zhenzhen.Xu@fda.hhs.gov> , Boguang Zhen<Boguang.Zhen@fda.hhs.gov>,
Yongsoek Park <yongpark@pitt.edu> and Bin Zhu <bin.zhu@nih.gov>
}
\keyword{package}