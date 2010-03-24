\name{p.ancova_to_es}
\alias{p.ancova_to_es}
\title{One or Two-tailed p-value from ANCOVA to Effect Size
}
\description{Converts a one or two-tailed p-value from ANCOVA to an effect size of  d (mean difference), g (unbiased estimate of d), r (correlation coefficient), z (Fisher's z), and log odds ratio. The variances of these estimates are also computed.
}
\usage{
p.ancova_to_es(p, n.1, n.2, R, q, tail = "two")
}
\arguments{
  \item{p}{One-tailed p-value reported in primary study.
}
   \item{n.1}{Treatment group sample size.
}
  \item{n.2}{Comparison group sample size.
}
  \item{R}{Covariate outcome correlation or multiple correlation.
}
  \item{q}{number of covariates.
}
\item{tail}{One or two-tailed p-value. Default is \code{two}.
}
}
\value{
\item{d}{Standardized mean difference (d).}
\item{var.d}{Variance of d.}
\item{g}{Unbiased estimate of d.}
\item{var.g}{Variance of g.}
\item{r}{Correlation coefficient.}
\item{var.r}{Variance of r.}
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
\code{\link{p_to_es}}
}
\keyword{ arith }