\name{fail_to_es}
\alias{fail_to_es}
\title{Failure groups to Effect Size
}
\description{Converts number of 'failure' groups to an effect size of  d (mean difference), g (unbiased estimate of d), r (correlation coefficient), z (Fisher's z), and log odds ratio. The variances of these estimates are also computed.
}
\usage{
fail_to_es(B, D, n.1, n.0)
}
\arguments{
  \item{B}{Treatment failure.
}
 \item{D}{Non-treatment failure.
}
  \item{n.1}{Treatment sample size.
}
  \item{n.0}{Control/comparison sample size.
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
\seealso{
\code{\link{lor_to_es}},
\code{\link{prop_to_es}}
}
\keyword{ arith }