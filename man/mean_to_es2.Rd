\name{mean_to_es2}
\alias{mean_to_es2}
\title{Means with Pooled SD to Effect Size
}
\description{Converts raw mean scores with pooled standard deviation to an effect size of  d (mean difference), g (unbiased estimate of d), r (correlation coefficient), z (Fisher's z), and log odds ratio. The variances of these estimates are also computed.
}
\usage{
mean_to_es2(m.1, m.2, s.pooled, n.1, n.2)
}
\arguments{
  \item{m.1}{Mean of group one.
}
  \item{m.2}{Mean of group two.
}
  \item{s.pooled}{Pooled standard deviation.
}
  \item{n.1}{Sample size of group one. 
}
  \item{n.2}{Sample size of group two. 
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
\code{\link{mean_to_es}},
\code{\link{mean_anc_to_es}},
\code{\link{mean_anc_to_es2}}
}
\keyword{ arith }