\name{p_to_es}
\alias{p_to_es}
\title{p-value to Effect Size
}
\description{One or two tailed p-value to an effect size of  d (mean difference), g (unbiased estimate of d), r (correlation coefficient), z (Fisher's z), and log odds ratio. The variances of these estimates are also computed.
}
\usage{
p_to_es(p, n.1, n.2, tail = "two")
}
\arguments{
  \item{p}{p-value. 
}
  \item{n.1}{Sample size of treatment group.
}
  \item{n.2}{Sample size of comparison group.
}
 \item{tail}{One or two-tailed p-value. Default is \code{two}.
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
\code{\link{p.ancova_to_es}}
}
\keyword{ arith }