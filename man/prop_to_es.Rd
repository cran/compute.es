\name{prop_to_es}
\alias{prop_to_es}
\title{Proportion to Effect Size
}
\description{Converts proportions (typically seen in studies reporting odds ratio statistics) to an effect size of  d (mean difference), g (unbiased estimate of d), r (correlation coefficient), z (Fisher's z), and log odds ratio. The variances of these estimates are also computed.
}
\usage{
prop_to_es(p1, p2, n.ab, n.cd)
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
\code{\link{fail_to_es}},
\code{\link{lor_to_es}}
}
\keyword{ arith }