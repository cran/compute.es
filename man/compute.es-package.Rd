\name{compute.es-package}
\alias{compute.es-package}
\docType{package}
\title{Compute Effect Sizes
}
\description{ At the outset of a meta-analytic endeavor, particularly when the literature search and study retrieval has been completed, the focus lies on deriving effect size estimates from included studies. These studies will report various statistical information to compute the estimates. For example, some studies will report means and standard deviations or mean differences, others will report correlations, and still others will report odds ratios. Can all these different estimates be combined for use in a meta-analysis? Yes and no. The statistical information reported in each study will need to be converted to a common index, whether that be the standardized mean difference, correlation, or binary estimate--it really depends on the type of meta-analysis being conducted. 

The goal in creating this package was to provide a comprehensive set of tools/functions to easily derive and/or convert statistics reported in the primary studies to all of the common effect size estimates and their variances. The \code{compute.es} package's functions will convert a variety of statistics, such as a t-test or p-value and sample size, to estimates of:

\tabular{ll}{
 (1) Cohen's \eqn{d} (mean difference) \cr
 (2) Hedges' \eqn{g} (unbiased estimate of \eqn{d})\cr
 (3) \eqn{r} (correlation coefficient)\cr
 (4) \eqn{z'} (Fisher's \eqn{z})\cr
 (5) log odds ratio \cr
 (6) the variances of these estimates \cr
 }

This package also contains functions to calculate \eqn{d} and \eqn{g} or \eqn{r} and \eqn{z'} and their variances for a full data set. The \code{compute.es} package uses recommended conversion formulas as described in \emph{The Handbook of Research Synthesis and Meta-Analysis} (Cooper, Hedges, & Valentine, 2009). 
}
\details{
\tabular{ll}{
Package: \tab compute.es\cr
Type: \tab Package\cr
Version: \tab 0.2\cr
Date: \tab 2010-04-01\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}

\bold{Structure of Functions}

The function names for this package are designed for quick processing, such that the first part of the function corresponds to the input method (statistical information reported in the study) and the remaining part corresponds to the output values, which are the effect size estimates ('es' at the end of each function). For example, the function

\code{des()}

has the input of a Cohen's \eqn{d} and will output various effect size ('es') estimates.

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
\author{AC Del Re 

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

tes(1.74, 30, 31)  

# where the reported t-value = 1.74, treatment sample 
# size = 30, and the control/comparison sample size = 31 will
# output effect sizes of d, g, r, z, and log odds ratio. 
# The variances will also be computed for these estimates.
# Note: If only the total sample size is reported simply split 
# the number in half for entry into the function. 

}
\seealso{
For information and user-friendly R packages to conduct a meta-analysis see:

Meta-Analysis with Correlations: 

\bold{MAc package}:  \url{http://CRAN.R-project.org/package=MAc}

Menu-Driven Meta-Analysis with Correlations:

\bold{RcmdrPlugin.MAc package}:  \url{http://CRAN.R-project.org/package=RcmdrPlugin.MAc}

Meta-Analysis with Mean Differences: 

\bold{MAd package}:  \url{http://CRAN.R-project.org/package=MAd}

Menu-Driven Meta-Analysis with Mean Differences: 

\bold{RcmdrPlugin.MAd package}:  \url{http://CRAN.R-project.org/package=RcmdrPlugin.MAd}
}