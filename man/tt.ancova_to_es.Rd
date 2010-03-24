\name{tt.ancova_to_es}
\alias{tt.ancova_to_es}
\title{t-test Value from ANCOVA to Effect Size
}
\description{Converts a t-test value from ANCOVA to an effect size of  d (mean difference), g (unbiased estimate of d), r (correlation coefficient), z (Fisher's z), and log odds ratio. The variances of these estimates are also computed.
}
\usage{
tt.ancova_to_es(t, n.1, n.2, R, q)
}
\arguments{
  \item{t}{t-test value reported in primary study.
}
   \item{n.1}{Treatment group sample size.
}
  \item{n.2}{Comparison group sample size.
}
  \item{R}{Covariate outcome correlation or multiple correlation.
}
  \item{q}{number of covariates.
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
\code{\link{tt_to_es}}
}
\keyword{ arith }