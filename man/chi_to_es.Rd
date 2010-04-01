\name{chies}
\alias{chies}
\title{Chi-Squared Statistic to Effect Size
}
\description{Converting Chi-squared (\eqn{\chi^2}) statistic with 1 degree of freedom to to an effect size of  \eqn{d} (mean difference), \eqn{g} (unbiased estimate of \eqn{d}), \eqn{r} (correlation coefficient), \eqn{z} (Fisher's \eqn{z}), and log odds ratio. The variances of these estimates are also computed.
}
\usage{
chies(chi.sq, n)
}
\arguments{
  \item{chi.sq}{Chi squared statistic from primary study.
}
  \item{n}{Sample size in primary study.
}
}
\details{
The chi-squared statistic (\eqn{\chi^2}) is defined as

\deqn{
  \chi^2=%
  \sum{\frac{(o-e)^2}%
       {e}}}{%
   chi^2=%
     sum(((o-e)^2)/%
     (e)}
where \eqn{o} is the observed value and \eqn{e} is the expected value. NOTE: This function requires the \eqn{\chi^2} value to have been derived with 1 degree of freedom (indicating 2 independent groups are used in the calculation).      
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

In this particular formula Cohen's \eqn{d} is calculated after \eqn{r} is computed and then derived from it

\deqn{d=%
       \frac{2r}%
         {\sqrt{1-r^2}}}{%
      d=%
       (2r)/%
          (sqrt(1-r^2))}

The variance of \eqn{d} is derived from

\deqn{v_{d}=%
        \frac{4v}%
        {(1-r^2)^3}}{%
      v_(d)=%
        (4v)/%
        ((1-r^2)^3)}             


The effect size estimate \eqn{d} has a small upward bias (overestimates the population parameter effect size) which can be removed using a correction formula to derive the unbiased estimate of Hedges' \eqn{g}. The correction factor, \eqn{j}, is defined as

\deqn{J=%
       1-%
         \frac{3}%
         {4df-1}}{%
      J=%
       1-%
         (3)/%
          (4df-1)}

where \eqn{df}= degrees of freedom, which is equal to 1 since the \eqn{\chi^2} degree of freedom = 1. Then, to calculate \eqn{g}

\deqn{g=%
       Jd}{%
      g=%
       Jd}
and the variance of \eqn{g}

\deqn{v_{g}=%
       J^2v_{d}}{%
       v_(g)=%
       J^2v_(d)}        

\bold{Correlation Coefficient r, Fisher's z, and Variances}:

In this particular formula \eqn{r} is calculated as follows
  
\deqn{r=%
       \sqrt{\frac{\chi^2}%
       {n}}}{%
       r=%
      sqrt(chi^2/%
      (n))}
      
where \eqn{\chi^2}{chi^2} is the chi-squared value with 1 degree of freedom and \eqn{n} is the total sample size. 

The variance of \eqn{r} is then defined as

\deqn{v_{r}=%
      \frac{(1-r^2)^2}%
      {n-1}}{%
      v_(r)=%
      ((1-r^2)^2)/%
      (n-1)}

Often researchers are interested in transforming \eqn{r} to \eqn{z'} (Fisher's \eqn{z}) because \eqn{r} is not normally distributed, particularly at large values of \eqn{r}. Therefore, converting to \eqn{z'} will help to normally distribute the estimate. Converting from \eqn{r} to \eqn{z'} is defined as

\deqn{z=%
      .5^*log(\frac{1+r}%
      {1-r})}{%
      z=%
      .5^*log((1+r)/%
      (1-r))}  

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
      (sqrt(3))}

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
\keyword{ arith }