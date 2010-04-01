\name{des}
\alias{des}
\title{Mean Difference (d) to Effect size
}
\description{Converts \eqn{d} (mean difference) to an effect size of \eqn{g} (unbiased estimate of \eqn{d}), \eqn{r} (correlation coefficient), \eqn{z} (Fisher's \eqn{z}), and log odds ratio. The variances of these estimates are also computed.
}
\usage{
des(d, n.1, n.2)
}
\arguments{
  \item{d}{Mean difference statistic (\eqn{d}).
 }
  \item{n.1}{Sample size of group one.
}
  \item{n.2}{Sample size of group one.
}
}
\details{\bold{Information regarding input (d):}

In a study comparing means from independent groups, the population standardized mean difference is defined as

\deqn{\delta=%
        \frac{\mu_{2}-\mu_{1}}%
          {\sigma}}{%
      delta=%
        (mu_(2)-mu_(1))/%
         (sigma)}
         
 where \eqn{\mu_{2}}{(mu_(2))} is the population mean of the second group, \eqn{\mu_{1}}{(mu_(1))} is the population mean of the first group, and \eqn{\sigma}{(sigma)} is the population standard deviation (assuming \eqn{\sigma_{2}}{(sigma_(2))} = \eqn{\sigma_{1}}{(sigma_(1))}).
 
The estimate of  \eqn{\delta}{delta} from independent groups is defined as

\deqn{d=%
        \frac{\bar Y_{2}-\bar Y_{1}}%
          {S_{within}}}{%
      d=%
        (Y_(2) bar-Y_(1) bar)/%
         (S_(within))}  
where \eqn{\bar Y_{2}}{Y_(2) bar} and \eqn{\bar Y_{1}}{Y_(1) bar} are the sample means in each group and \eqn{S_{within}}{S_(within)} is the standard deviation pooled across both groups and is defined as
 
\deqn{S_{within}=%
        \sqrt{\frac{(n_{1}-1)S^2_{1}+(n_{2}-1)S^2_{2}}%
          {n_{1}+n_{2}-2}}}{%
      S_(within)=%
        ((sqrt((n_(1)-1)S^2_(1)+(n_(2)-1)S^2_(2)/%
         (n_(1)+n_(2)-2)}
               
where \eqn{n_{1}}{n_(1)} and \eqn{n_{2}}{n_(2)} are the sample sizes of group 1 and 2 respectively and \eqn{S^2_{1}}{S^2_(1)} and \eqn{S^2_{2}}{S^2_(2)} are the standard deviations of each group. The variance of \eqn{d}{d} is then defined as

\deqn{v_{d}=%
       \frac{n_{1}+n_{2}}%
       {n_{1}n_{2}}+%
       \frac{d^2}%
       {2(n_{1}n_{2})}}{%
       v_(d)=%
       (n_(1)+n_(2))/%
       (n_(1)n_(2))+%
       (d^2)/%
       (2(n_(1)n_(2)))} 
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

   (1) \emph{Hedges' \eqn{g} (unbiased estimate of \eqn{d}) and variance}
   
   (2) \emph{Correlation coefficient (\eqn{r}), Fisher's \eqn{z}, and variance}
   
   (3) \emph{Log odds and variance}

is provided below:    
   
\bold{Hedges' g and Variance of g}:


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
\keyword{ arith }

