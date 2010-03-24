\name{chi_to_es}
\alias{chi_to_es}
\title{Chi-Squared Statistic to Effect Size
}
\description{Converting Chi-squared statistic with 1 df to to an effect size of  d (mean difference), g (unbiased estimate of d), r (correlation coefficient), z (Fisher's z), and log odds ratio. The variances of these estimates are also computed.
}
\usage{
chi_to_es(chi.sq, n)
}
\arguments{
  \item{chi.sq}{Chi squared statistic from primary study.
}
  \item{n}{Sample size in primary study.
}
}
\value{
\item{d}{Standardized mean difference (d).}
\item{var.d }{Variance of d.}
\item{g }{Unbiased estimate of d.}
\item{var.g }{Variance of g.}
\item{r }{Correlation coefficient.}
\item{var.r }{Variance of r.}
\item{log_odds}{Log odds ratio.}
\item{var.log_odds}{Variance of log odds ratio.}
\item{n}{Total sample size.}
}
\author{ AC Del Re 

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\references{Borenstein (2009). Effect sizes for continuous data. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 279-293). New York: Russell Sage Foundation.
}
\keyword{ arith }