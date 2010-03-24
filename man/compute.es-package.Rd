\name{compute.es-package}
\alias{compute.es-package}
\docType{package}
\title{Compute Effect Sizes
}
\description{ This package contains several functions that will convert a variety of statistics, such as a t-test or p-value and sample size, to effect sizes of d (mean difference), g (unbiased estimate of d), r (correlation coefficient), z (Fisher's z), and log odds ratio. The variances of these estimates are also computed. These conversion functions are a particularly useful resource during the preliminary stages of a meta-analytic project, when coding for relevant effect size data (and the appropriate statistic is not reported). This package uses recommended coversion formulas as described in The Handbook of Research Synthesis and Meta-Analysis (Cooper, Hedges, & Valentine, 2009). 
}
\details{
\tabular{ll}{
Package: \tab compute.es\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2010-03-24\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}

}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\references{ 
Borenstein (2009). Effect sizes for continuous data. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 279-293). New York: Russell Sage Foundation.

Cooper, H., Hedges, L.V., & Valentine, J.C. (2009). \emph{The handbook of research synthesis and meta-analysis} (2nd edition). New York: Russell Sage Foundation.
}
\keyword{ package }
\examples{ 

# 1. Computations to Calculate Effect Sizes:
 
# For example, suppose the primary study reported a t-test 
# value for differences between 2 groups. Then, running:

tt_to_es(1.74, 30, 31)  

# where the reported t-value = 1.74, treatment sample 
# size = 30, and the control/comparison sample size = 31 will
# output effect sizes of d, g, r, and log odds ratio. 
# The variances will also be computed for these estimates.
# Note: If only the total sample size is reported simply split 
# the number in half for entry into the function. 

}
\seealso{
\bold{MAc} package:  \url{http://CRAN.R-project.org/package=MAc};
\bold{RcmdrPlugin.MAc} package:  \url{http://CRAN.R-project.org/package=RcmdrPlugin.MAc};
\bold{MAd} package:  \url{http://CRAN.R-project.org/package=MAd};
\bold{RcmdrPlugin.MAd} package:  \url{http://CRAN.R-project.org/package=RcmdrPlugin.MAd}
}