\name{d_to_es}
\alias{d_to_es}
\title{Mean Difference (d) to Effect size
}
\description{Converts d (mean difference) to an effect size of g (unbiased estimate of d), r (correlation coefficient), z (Fisher's z), and log odds ratio. The variances of these estimates are also computed.
}
\usage{
d_to_es(d, n.1, n.2)
}
\arguments{
  \item{d}{Mean difference statistic (d).
 }
  \item{n.1}{Sample size of group one.
}
  \item{n.2}{Sample size of group one.
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

