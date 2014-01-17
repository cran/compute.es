\name{res}
\alias{res}
\title{Correlation coefficient (r) to Effect Size
}
\description{Converts correlation (r) to an effect size of  \eqn{d} (mean difference), \eqn{g} (unbiased estimate of \eqn{d}), \eqn{r} (correlation coefficient), \\eqn{z'} (Fisher's \eqn{z}), and log odds ratio. The variances, confidence intervals and p-values of these estimates are also computed, along with NNT (number needed to treat), U3 (Cohen's \eqn{U_(3)} overlapping proportions of distributions), CLES (Common Language Effect Size) and Cliff's Delta.
}
\usage{
res(r, var.r = NULL, n, level=95, dig=2, id=NULL, data=NULL)
}
\arguments{
  \item{r}{Correlation coefficient.
 }
  \item{var.r}{Variance of r. If value is not reported then leave it blank and variances will be computed based on sample size. Otherwise, enter this value (e.g., \code{r_to_es(.27, var.r = .02, 30)}.
}
  \item{n}{Total sample size.
}
\item{level}{Confidence level. Default is \code{95\%}.
}
\item{dig}{Number of digits to display. Default is \code{2} digits.
}
\item{id}{Study identifier. Default is \code{NULL}, assuming a scalar is used as input. If input is a vector dataset (i.e., \code{data.frame}, with multiple values to be computed), enter the name of the study identifier here.
}
\item{data}{name of \code{data.frame}. Default is \code{NULL}, assuming a scalar is used as input. If input is a vector dataset (i.e., \code{data.frame}, with multiple values to be computed), enter the name of the \code{data.frame} here.
}
}
\value{
\item{d}{Standardized mean difference (\eqn{d}).}
\item{var.d }{Variance of \eqn{d}.}
\item{l.d }{lower confidence limits for \eqn{d}.}
\item{u.d }{upper confidence limits for \eqn{d}.}
\item{U3.d }{Cohen's \eqn{U_(3)}, for \eqn{d}.}
\item{cl.d }{ Common Language Effect Size for \eqn{d}.}
\item{cliffs.d }{Cliff's Delta for \eqn{d}.}
                 \item{p.d }{p-value for \eqn{d}.}
                 \item{g }{Unbiased estimate of \eqn{d}.}
                 \item{var.g }{Variance of \eqn{g}.}
                 \item{l.g }{lower confidence limits for \eqn{g}.}
                 \item{u.g }{upper confidence limits for \eqn{g}.}
                 \item{U3.g }{Cohen's \eqn{U_(3)}, for \eqn{g}.}
\item{cl.g }{ Common Language Effect Size for \eqn{g}.}
\item{p.g }{p-value for \eqn{g}.}
\item{r }{Correlation coefficient.}
\item{var.r }{Variance of \eqn{r}.}
\item{l.r }{lower confidence limits for \eqn{r}.}
\item{u.r }{upper confidence limits for \eqn{r}.}
\item{p.r }{p-value for \eqn{r}.}
\item{z }{Fisher's z (\eqn{z'}).}
\item{var.z }{Variance of \eqn{z'}.}
\item{l.z }{lower confidence limits for \eqn{z'}.}
\item{u.z }{upper confidence limits for \eqn{z'}.}
\item{p.z}{p-value for \eqn{z'}.}
\item{OR}{Odds ratio.}
\item{l.or }{lower confidence limits for \eqn{OR}.}
\item{u.or }{upper confidence limits for \eqn{OR}.}
\item{p.or}{p-value for \eqn{OR}.}
\item{lOR}{Log odds ratio.}
\item{var.lor}{Variance of log odds ratio.}
\item{l.lor }{lower confidence limits for \eqn{lOR}.}
\item{u.lor }{upper confidence limits for \eqn{lOR}.}
\item{p.lor}{p-value for \eqn{lOR}.}
\item{N.total}{Total sample size.}
\item{NNT}{Number needed to treat.}
}

\note{
  
  \bold{Detailed information regarding output values of:}
  
  (1) \emph{Cohen's \eqn{d}, Hedges' \eqn{g} (unbiased estimate of \eqn{d}) and variance}
  
  (2) \emph{Correlation coefficient (\eqn{r}), Fisher's \eqn{z'}, and variance}

(3) \emph{Log odds and variance}

is provided below (followed by general information about NNT, U3, Common Language Effect Size, and Cliff's Delta):
  
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



\bold{General information about NNT, U3, Common Language Effect Size, and Cliff's Delta:}

\emph{Number needed to treat (NNT).} NNT is interpreted as the number of participants that would need to be treated in one group (e.g., intervention group) in order to have one additional positive outcome over that of the outcome of a randomly selected participant in the other group (e.g., control group). In the \code{compute.es} package, NNT is calculated directly from d (Furukawa & Leucht, 2011), assuming relative normality of distribution and equal variances across groups, as follows:
  
  \deqn{NNT=%
        \frac{1}%
{ 2 * \Phi(\frac{d}%
             {\sqrt{2}})-1}
  }{
    NNT=%
    1/(2*Phi(d/sqrt(2))-1)       
  }

\emph{U3.} Cohen (1988) proposed a method for characterizing effect sizes by expressing them in terms of (normal) distribution overlap, called U3. This statistic describes the percentage of scores in one group that are exceeded by the mean score in another group. If the population means are equal then half of the scores in the treatment group exceed half the scores in the comparison group, and U3 = 50\%. As the population mean difference increases, U3 approaches 100\% (Valentine & Cooper, 2003). 

\emph{Common Language Effect Size (CLES).} CLES (McGraw & Wong, 1992) expresses the probability that a randomly selected score from one population will be greater than a randomly sampled score from another population. CLES is computed as the percentage of the normal curve that falls between negative infinity and the effect size (Valentine & Cooper, 2003).  

\emph{Cliff's Delta/success rate difference.} Cliff's delta (or success rate difference; Furukawa & Leucht (2011)) is a robust alternative to Cohen's d, when data are either non-normal or ordinal (with truncated/reduced variance). Cliff's Delta is a non-parametric procedure that provides the probability that individual observations in one group are likely to be greater than the observations in another group. It is the probability that a randomly selected participant of one population has a better outcome than a randomly selected participant of the second population (minus the reverse probability). Cliff's Delta of negative 1 or positive 1 indicates no overlap between the two groups, whereas a value of 0 indicates complete overlap and equal group distributions. 
      
      \deqn{\delta=%
      2 * \Phi(\frac{d}%
             {\sqrt{2}})-1
      }{
        Cliff's Delta=%
        2*Phi(d/sqrt(2))-1       
      }
        
}
        
        \author{ AC Del Re
        
        Much appreciation to Dr. Jeffrey C. Valentine for his contributions in implementing \eqn{U3} and \eqn{CLES} procedures and related documentation.
        
        Maintainer: AC Del Re \email{acdelre@gmail.com}
        }
        \references{Borenstein (2009). Effect sizes for continuous data. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 279-293). New York: Russell Sage Foundation.
        
        Cohen, J. (1988). \emph{Statistical power for the behavioral sciences (2nd ed.)}. Hillsdale, NJ: Erlbaum.
        
        Furukawa, T. A., & Leucht, S. (2011). How to obtain NNT from Cohen's d: comparison of two methods. \emph{PloS one, 6}(4), e19070.

McGraw, K. O. & Wong, S. P. (1992). A common language effect size statistic. \emph{Psychological Bulletin, 111,} 361-365.

Valentine, J. C. & Cooper, H. (2003). \emph{Effect size substantive interpretation guidelines: 
Issues in the interpretation of effect sizes.} Washington, DC: What Works Clearinghouse.

}
\examples{ 
# CALCULATE SEVERAL EFFECT SIZES BASED ON CORRELATION STATISTIC: 

res(.3, n=30)
}
\keyword{ arith }