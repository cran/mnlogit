\name{Hypothesis Testing}
\alias{lrtest}
\alias{waldtest}
\alias{scoretest}
\alias{hmftest}
\alias{lrtest.mnlogit}
\alias{waldtest.mnlogit}
\alias{scoretest.mnlogit}
\alias{hmftest.formula}
\alias{hmftest.mnlogit}

\title{Hypothesis testing for multinomial logit models}

\description{
Three hypothesis tests applicable to any MLE (Likelihood ratio test, Wald test, Rao score test) and the Hausman-McFadden test for IIA (independence of irrelevant alternatives) are provided. 
}

\usage{
lrtest(object, ...)
waldtest(object, ...)
scoretest(object, ...)
\method{lrtest}{mnlogit}(object, ...)
\method{scoretest}{mnlogit}(object, ...)
\method{waldtest}{mnlogit}(object, ...)

hmftest(x, ...)
\method{hmftest}{formula}(x, alt.subset, ...)
\method{hmftest}{mnlogit}(x, z, ...)
}
\arguments{
    \item{object}{An fitted model which is an object of class \code{mnlogit}.}
    \item{...}{For \code{lrtest} and \code{waldtest} a fitted \code{mnlogit} object or a \code{formula} object maybe given. However for \code{scoretest} ONLY an fitted \code{mnlogit} is accepted. For \code{hmftest} either a subset of alternatives or an \code{mnlogit} object estimated using a subset of alternatives must be given.}
    \item{x}{A fitted object of class \code{mnlogit} or a \code{formula} object.}
    \item{z}{An object of class \code{mnlogit} or a subset of alternatives for the \code{hmftest.mnlogit} method. Must be the same model as \code{x} estimated on a subset of alternatives. }
    \item{alt.subset}{A subset of alternatives to do testing on.}
}   

\value{
    An object of class \code{htest}, with elements:
    \item{statistic}{The value of the test statistic.}
    \item{parameter}{the value of the underlying test distribution's parameter. In this case, the numbe r of degrees of freedom of chi-squared distribution.}
    \item{data.name}{The data frame used.}
    \item{p.value}{Probability for accepting the null hypothesis.}
    \item{method}{The name of the hypothesis test.}
    \item{alternative}{Alternative hypothesis.}
}

\references{
Code for the last two tests (Score and IIA test) is gratefully adapted from the CRAN package \pkg{mlogit}, while the first 2 tests are performed by functions in the CRAN package \pkg{lmtest}. 

  Croissant, Yves. 
    \emph{Estimation of multinomial logit models in R: The mlogit Packages.}
    \url{http://cran.r-project.org/web/packages/mlogit/index.html}

   Achim Zeileis, Torsten Hothorn (2002)
   \emph{Diagnostic Checking in Regression Relationships},
   R News 2(3), 7-10. \url{http://CRAN.R-project.org/doc/Rnews/}.
}

\author{Asad Hasan, Wang Zhiyu, Alireza S. Mahani}

\examples{
  library(mnlogit)
  data(Fish, package = "mnlogit")

  # Unconstrained model
  fm <- formula(mode ~ price | income | catch)
  fit <- mnlogit(fm, Fish) 
  # Constrained model - intercep dropped
  fm.c <- formula(mode ~ price | income - 1 | catch)
  fit.c <- mnlogit(fm.c, Fish) 

  ## MLE hypothesis tests
  lrtest(fit, fit.c)
  waldtest(fit, fit.c)
  scoretest(fit, fit.c)

  ## IIA test
  alt.subset <- c("beach", "boat", "charter")
  hmftest(fit, alt.subset)
}
