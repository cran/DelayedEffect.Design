
\name{GPW.logrank}
\alias{GPW.logrank}
\title{Generalized Piecewise Weighted Logrank Test}
\description{
  Compute the p-value based on the 
 generalized piecewise weighted log-rank test when the treatment 
 time-lag effect is present and the lag duration 
varies heterogeneously from individual to individual or from study to study, within a
certain domain and following a uniform pattern.                      
}
\usage{GPW.logrank(data, obs.time, time.to.event, event.status, trt.group, tl, tu) 
}
\arguments{
  \item{data}{Data frame}
  \item{obs.time}{Column name in \code{data} for the observational time.}
  \item{time.to.event}{Column name in \code{data} for the event time.}
  \item{event.status}{Column name in \code{data} for the event status, where
                      0 denotes being censored, and 1 denotes events.}
  \item{trt.group}{Column name in \code{data} for the treatment group, where
                      0 denotes controls, and 1 denotes treated subjects.}
  \item{tl}{Lower bound of delayed duration domain}
  \item{tu}{Upper bound of delayed duration domain}  
}
%\details{}
\value{ The p-value of the test.
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
  data(data, package="DelayedEffect.Design")
  GPW.logrank(data, "X", "time_to_event", "event_status", "Z", 30, 30*11) 
}
\keyword{ test } 
\seealso{\code{\link{pow.SEPPLE.plus}} }
