\name{compute.es-package}
\alias{compute.es-package}
\docType{package}
\title{Compute Effect Sizes in R
}
\description{This package provides a comprehensive set of tools/functions to easily derive and/or convert statistics generated from one's study (or from those reported in a published study) to all of the common effect size estimates, along with their variances, confidence intervals, and p-values. Several additional statistics are generated, including NNT (number needed to treat), U3 (Cohen's U3 distribution overlap statistic), CLES (Common Language Effect Size) and Cliff's Delta (success rate difference). The \code{compute.es} package's functions will convert a variety of statistics, such as means and standard deviations, t-test or p-value and sample size, to estimates of:

\tabular{ll}{
 (1) Cohen's \eqn{d} (mean difference) \cr
 (2) Hedges' \eqn{g} (unbiased estimate of \eqn{d})\cr
 (3) \eqn{r} (correlation coefficient)\cr
 (4) \eqn{z'} (Fisher's \eqn{z})\cr
 (5) log odds ratio \cr
 (6) the variances, confidence intervals and p-values of the above estimates \cr
 (7) Other statistics: NNT, U3, CLES, Cliff's Delta \cr
 }

The functions in this package can compute the effect sizes from a single study or from multiple studies simultaneously.  The \code{compute.es} package uses recommended conversion formulas as described in \emph{The Handbook of Research Synthesis and Meta-Analysis} (Cooper, Hedges, & Valentine, 2009). 
}
\details{
\tabular{ll}{
Package: \tab compute.es\cr
Type: \tab Package\cr
Version: \tab 0.2-4\cr
Date: \tab 2014-09-16\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}

\bold{Structure of Functions}

The function names for this package are designed for quick processing, such that the first part of the function corresponds to the input method (statistical information reported in the study) and the remaining part corresponds to the output values, which are the effect size estimates ('es' at the end of each function). For example, the function \code{des()} has the input of a Cohen's \eqn{d} and will output various effect size ('es') estimates.

The other function inputs and names are as follows:

\tabular{ll}{
ANCOVA F-test:  \tab \code{a.fes()} \cr  
ANCOVA means: \tab \code{a.mes()} \cr
ANCOVA means (pooled \eqn{sd}): \tab \code{a.mes2()} \cr
ANCOVA p-value: \tab \code{a.pes()}  \cr
ANCOVA t-test: \tab \code{a.tes()}  \cr
Chi-squared (1 \eqn{df}): \tab \code{chies()} \cr
Correlation: \tab \code{res()} \cr
d-statistic: \tab \code{des()}  \cr
Failure group (binary): \tab \code{failes()} \cr
F-test: \tab \code{fes()} \cr
Log odds ratio: \tab \code{lores()} \cr
Means: \tab \code{mes()}  \cr
Means (pooled \eqn{sd}): \tab \code{mes2()} \cr
Proportions (binary): \tab \code{propes()} \cr
p-value: \tab \code{pes()} \cr
t-test: \tab \code{tes()} \cr
}
}
\author{AC Del Re with contributions from Jeffrey C. Valentine

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\references{ 
Borenstein (2009). Effect sizes for continuous data. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 279-293). New York: Russell Sage Foundation.

Cooper, H., Hedges, L.V., & Valentine, J.C. (2009). \emph{The handbook of research synthesis and meta-analysis} (2nd edition). New York: Russell Sage Foundation.

Furukawa, T. A., & Leucht, S. (2011). How to obtain NNT from Cohen's d: comparison of two methods. \emph{PloS one, 6}(4), e19070.

McGraw, K. O. & Wong, S. P. (1992). A common language effect size statistic. \emph{Psychological Bulletin, 111,} 361-365.

Valentine, J. C. & Cooper, H. (2003). \emph{Effect size substantive interpretation guidelines: 
Issues in the interpretation of effect sizes.} Washington, DC: What Works Clearinghouse.

}
\keyword{ package }
\examples{ 

## 1. Computations to Calculate Effect Sizes:
 
# For example, suppose the primary study reported a t-test 
# value for differences between 2 groups. Then, running:

tes(t=1.74, n.1=30, n.2=31)

# Or, more simply:

tes(1.74, 30, 31)  

# where the reported t-value = 1.74, treatment sample 
# size = 30, and the control/comparison sample size = 31 will
# output effect sizes of d, g, r, z, OR, and log odds ratio. 
# The variances, confidence intervals, p-values and other 
# statistics will also be computed.
# Note: If only the total sample size is reported simply split 
# the number in half for entry into the function. 

# Now suppose one has a dataset (i.e., data.frame in R-speak)
# with several t-values to be converted into effect sizes:

# First, we will generate sample data:

dat <- data.frame(id=1:5,t=rnorm(5, 2, .5), 
                  n.t=round(rnorm(5, 25),0), 
                  n.c=round(rnorm(5, 25),0))

# Running the fuction as follows will generate a new 
# data.frame with several effect size estimates

tes(t=t, n.1=n.t, n.2=n.c, level=95, dig=2, id=id, data=dat)

}
\seealso{
For information and user-friendly R packages to conduct a meta-analysis see:

Menu-Driven Meta-Analysis (Graphical User Interface): 

\bold{RcmdrPlugin.MA package}:  \url{https://CRAN.R-project.org/package=RcmdrPlugin.MA}

Meta-Analysis with Correlations: 

\bold{MAc package}:  \url{https://CRAN.R-project.org/package=MAc}

Meta-Analysis with Mean Differences: 

\bold{MAd package}:  \url{https://CRAN.R-project.org/package=MAd}
}
