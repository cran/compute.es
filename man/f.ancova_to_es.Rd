\name{a.fes}
\alias{a.fes}
\title{ANCOVA F-statistic to Effect Size
}
\description{
Converts an ANCOVA \eqn{F} to an effect size of  \eqn{d} (mean difference), \eqn{g} (unbiased estimate of \eqn{d}), \eqn{r} (correlation coefficient), \eqn{z'} (Fisher's \eqn{z}), and log odds ratio. The variances of these estimates are also computed.
}
\usage{
a.fes(f, n.1, n.2, R, q)
}
\arguments{
  \item{f}{\eqn{F} value from ANCOVA.
}
  \item{n.1}{Treatment group sample size.
}
  \item{n.2}{Comparison group sample size.
}
  \item{R}{Covariate outcome correlation or multiple correlation.
}
  \item{q}{number of covariates.
}
}
\value{
\item{d}{Standardized mean difference (\eqn{d}).}
\item{var.d }{Variance of \eqn{d}.}
\item{g }{Unbiased estimate of \eqn{d}.}
\item{var.g }{Variance of \eqn{g}.}
\item{r }{Correlation coefficient.}
\item{var.r }{Variance of \eqn{r}.}
\item{log.odds}{Log odds ratio.}
\item{var.log.odds}{Variance of log odds ratio.}
\item{n}{Total sample size.}
}
\note{
\bold{Detailed information regarding output values of}

   (1) \emph{Cohen's \eqn{d}, Hedges' \eqn{g} (unbiased estimate of \eqn{d}) and variance}
   
   (2) \emph{Correlation coefficient (\eqn{r}), Fisher's \eqn{z}, and variance}
   
   (3) \emph{Log odds and variance}

is provided below:   
   
\bold{Cohen's d, Hedges' g and Variance of g}:

In this particular formula Cohen's \eqn{d} is calculated from the ANCOVA \eqn{F} with independent groups

\deqn{d=%
       \sqrt{\frac{F(n_{1}+n_{2})}%
         {n_{1}n_{2}}}%
         \sqrt{1-R^2}}{%
      d=%
       (F(n_(1)+n_(2)))/%
          (n_(1)n_(2))%
          sqrt(1-R^2)}

The variance of \eqn{d} is derived from

\deqn{v_{d}=%
        \frac{(n_{1}+n_{2})(1-R^2)}%
        {n_{1}n_{2}}+%
        \frac{d^2}%
        {2(n_{1}+n_{2})}}{%
      v_(d)=%
        ((n_(1)+n_(2))(1-R^2))/%
        (n_(1)n_(2))+%
        (d^2)/%
        (2(n_(1)+n_(2)))}  


The effect size estimate \eqn{d} has a small upward bias (overestimates the population parameter effect size) which can be removed using a correction formula to derive the unbiased estimate of Hedges' \eqn{g}. The correction factor, \eqn{j}, is defined as

\deqn{J=%
       1-%
         \frac{3}%
         {4df-1}}{%
      J=%
       1-%
         (3)/%
          (4df-1)}

where \eqn{df}= degrees of freedom, which is \eqn{n_{1}+n_{2}-2}{n_(1)+n_(2)-2} for two independent groups. Then, to calculate \eqn{g}

\deqn{g=%
       Jd}{%
      g=%
       Jd }
and the variance of \eqn{g}

\deqn{v_{g}=%
       J^2v_{d}}{%
       v_(g)=%
       J^2v_(d)}        

      
\bold{Correlation Coefficient r, Fisher's z, and Variances}:

In this particular formula \eqn{r} is calculated as follows

\deqn{r=%
      \frac{d}%
      {\sqrt{d^2+a}}}{%
       r=%
      (d)/%
      (sqrt(d^2+a))}
      
where \eqn{a} corrects for inbalance in \eqn{n_{1}}{n_(1)} & \eqn{n_{2}}{n_(2)} and is defined as

\deqn{a=%
      \frac{(n_{1}+n_{2})^2}%
      {n_{1}n_{2}}}{%
       a=%
      ((n_(1)+n_(2))^2)/%
      (n_(1)n_(2))}

The variance of \eqn{r} is then defined as

\deqn{v_{r}=%
      \frac{a^2v_{d}}%
      {(d^2+a)^3}}{%
      v_(r)=%
      (a^2v_(d))/%
      ((d^2+a)^3)}

Often researchers are interested in transforming \eqn{r} to \eqn{z'} (Fisher's \eqn{z}) because \eqn{r} is not normally distributed, particularly at large values of \eqn{r}. Therefore, converting to \eqn{z'} will help to normally distribute the estimate. Converting from \eqn{r} to \eqn{z'} is defined as

\deqn{z=%
      .5^*log(\frac{1+r}%
      {1-r})}{%
      z=%
      .5^*log((1+r)/%
      (1-r)}  

and the variance of \eqn{z}

\deqn{v_{z}=%
      \frac{1}%
      {n-3}}{%
      v_(z)=%
      (1)/%
      (n-3)}

where \eqn{n} is the total sample size for groups 1 and 2.                         


\bold{Log Odds Ratio & Variance of Log Odds}:

In this particular formula, log odds is calculated as follows

\deqn{\log(o)=%
      \frac{\pi d}%
      {\sqrt{3}}}{%
      log(o)=%
      (pi d)/%
      (sqrt(3)}

where \eqn{pi} = 3.1459. The variance of log odds is defined as

\deqn{v_{log(o)}=%
      \frac{\pi^2v_{d}}%
      {3}}{%
      v_(log(o))=%
      (pi^2v_(d))/%
      (3)}
   

}
\author{ AC Del Re 

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\references{Borenstein (2009). Effect sizes for continuous data. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 279-293). New York: Russell Sage Foundation.
}
\seealso{
\code{\link{fes}}
}
\keyword{ arith }