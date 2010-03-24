\name{lor_to_es}
\alias{lor_to_es}
\title{Log Odds Ratio to Standardized Mean Difference (d)
}
\description{Converts a log odds ratio to an effect size of  d (mean difference), g (unbiased estimate of d), r (correlation coefficient), z (Fisher's z), and original log odds value. The variances of these estimates are also computed.
}
\usage{
lor_to_es(lor, var.lor, n.1, n.2)
}
\arguments{
  \item{lor}{Log odds ratio reported in the primary study.
}
  \item{var.lor}{Variance of the log odds ratio.
}
 \item{n.1}{Sample size of treatment group.
}
  \item{n.2}{Sample size of comparison group.
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
\code{\link{prop_to_es}},
\code{\link{fail_to_es}}
}
\keyword{ arith }