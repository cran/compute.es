\name{propes}
\alias{propes}
\title{Proportions to Effect Size
}
\description{Converts proportions (typically seen in studies reporting odds ratios) to an effect size of  \eqn{d} (mean difference), \eqn{g} (unbiased estimate of \eqn{d}), \eqn{r} (correlation coefficient), \eqn{z'} (Fisher's \eqn{z}), and log odds ratio. The variances of these estimates are also computed.
}
\usage{
propes(p1, p2, n.ab, n.cd)
}
\arguments{
  \item{p1}{Proportion one.
}
 \item{p2}{Proportion two.
}
  \item{n.ab}{Total sample size for group A and B.
}
  \item{n.cd}{Total sample size for group C and D.
}
}
\details{
This formula will first compute an odds ratio and then transform to log odds and its variance. Then, Cohen's \eqn{d} will be calculated and this value will then be used to compute the remaining effect size estimates. The odds ratio is derived as follows
  
\deqn{
  or=%
    \frac{p_{1}(1-p_{2})}%
    {p_{2}(1-p_{1})}}{%
    or=%
    (p_(1)(1-p_(2)))/%
    (p_(2)(1-p_(1)))}


The conversion to a log odds and its variance is defined as
\deqn{ln(o)=%
       log(or)}

\deqn{v_{ln(o)}=%
       \frac{1}%
       {n_{AB}p_{1}(1-p_{1})}+%
        \frac{1}%
        {n_{CD}p_{2}(1-p_{2})}}{%
        v_(ln(o))=%
       (1)/%
       (n_(AB)p_(1)(1-p_(1)))+%
        (1)/%
        (n_(CD)p_(2)(1-p_(2)))}

where \eqn{n_{AB}}{n_(AB)} is the sum of group A and B sample size, \eqn{n_{CD}}{n_(CD)} is the sum of group C and D sample size, \eqn{p_{1}}{p_(1)} is the proportion for group 1 and \eqn{p_{2}}{p_(2)} is the proportion for group 2.       
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

In this particular formula Cohen's \eqn{d} is calculated from the log odds as follows

\deqn{d=%
       \frac{ln(o)\sqrt{3}}%
         {\pi}}{%
      d=%
       (ln(o)sqrt(3))/%
          (pi)}

The variance of \eqn{d} is derived from

\deqn{v_{d}=%
        \frac{3v_{ln(o)}}%
        {\pi^2}}{%
      v_(d)=%
        (3v_(ln(o)))/%
        (pi^2)}  

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



}
\author{ AC Del Re 

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\references{Borenstein (2009). Effect sizes for continuous data. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 279-293). New York: Russell Sage Foundation.
}
\seealso{
\code{\link{failes}},
\code{\link{lores}}
}
\keyword{ arith }