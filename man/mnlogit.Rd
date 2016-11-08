\name{mnlogit}
\alias{mnlogit}
\alias{fitted.mnlogit}
\alias{residuals.mnlogit}
\alias{df.residual.mnlogit}
\alias{terms.mnlogit}
\alias{update.mnlogit}
\alias{print.mnlogit}
\alias{vcov.mnlogit}
\alias{logLik.mnlogit}
\alias{summary.mnlogit}
\alias{print.summary.mnlogit}
\alias{index}
\alias{index.mnlogit}
\alias{predict.mnlogit}
\alias{coef.mnlogit}
\alias{print.est.stats}
\alias{print.model.size}

\title{Fast estimation of multinomial logit models}

\description{
Time and memory efficient estimation of multinomial logit models using maximum likelihood method. Targeted at large scale multiclass classification problems in econometrics and machine learning. 
Numerical optimization is performed by the Newton-Raphson method using an optimized, parallel C++ library to achieve fast computation of Hessian matrices.
The user interface closely related to the CRAN package \pkg{mlogit}.
}

\usage{
mnlogit(formula, data, choiceVar=NULL, maxiter = 50, ftol = 1e-6,
        gtol = 1e-6, weights = NULL, ncores = 1, na.rm = TRUE, 
        print.level=0, linDepTol = 1e-6, start=NULL, alt.subset=NULL, ...)
\method{fitted}{mnlogit}(object, outcome=TRUE, ...)
\method{residuals}{mnlogit}(object, outcome=TRUE, ...)
\method{df.residual}{mnlogit}(object, ...)
\method{terms}{mnlogit}(x, ...)
\method{update}{mnlogit}(object, new, ...)
\method{print}{mnlogit}(x, digits = max(3, getOption("digits") - 2),
                        width = getOption("width"), 
                        what = c("obj", "eststat", "modsize"), ...)
\method{vcov}{mnlogit}(object, ...)
\method{logLik}{mnlogit}(object, ...)
\method{summary}{mnlogit}(object, ...)
\method{print.summary}{mnlogit}(x, digits = max(3, getOption("digits") - 2),
                       width = getOption("width"), ... )
\method{index}{mnlogit}(object, ...)
\method{predict}{mnlogit}(object, newdata = NULL, probability = TRUE, 
                          returnData=FALSE, choiceVar=NULL, ...)
\method{coef}{mnlogit}(object, order=FALSE, as.list = FALSE, ...)
}
\arguments{
    \item{formula}{\code{formula} object or string specifying the model to be estimated (see Note).}
    \item{data, newdata}{A \code{data.frame} object with data organized in the 'long' format (see Note).This can also be a \code{mlogit.data} class object. \code{newdata} is used in the \code{predict} method.}
    \item{choiceVar}{A string naming the column in 'data' which has the list of choices. \emph{Note:} This argument is not used if \code{data} or \code{newdata} is a \code{mlogit.data} object.}
    \item{maxiter}{An integer indicating maximum number of Newton's iterations. If \code{maxiter <= 0}, then only Hessian, gradient and the loglikelihood are calculated at initial point.}
    \item{ftol}{A real number indicating tolerance on the difference of two subsequent loglikelihood values.}
    \item{gtol}{A real number indicating tolerance on norm of the gradient.}
    \item{weights}{Optional vector of (positive) frequency weights, one for each observation.}
    \item{ncores}{An integer indicating number of processors allowed for Hessian calculations.}
    \item{na.rm}{a logical variable which indicates whether rows of the data frame containing NAs will be removed.}
    \item{print.level}{An integer which controls the amount of information to be printed during execution.}
    \item{linDepTol}{Tolerance for detecting linear dependence between columns in input data. Dependent columns are removed from the estimation.}
    \item{start}{Named vector of coefficients to use as initial guess. Use naming convention as given by \code{names(coeffit())}, where \code{fit} is a \code{mnlogit} class object.}
    \item{alt.subset}{Subset of alternatives to perform estimation on.}
    \item{...}{Currently unused.}
    \item{object, x}{An object of class \code{mnlogit}.}
    \item{outcome}{ a boolean which indicates, for the \code{fitted} and the \code{residuals} methods whether a matrix (for each choice, one value for each alternative) or a vector (for each choice, only a value for the alternative chosen) should be returned.}
    \item{new}{An \code{formula} for the \code{update} method. It must obey all rules specified for the \code{formula} argument.}
    \item{digits}{Number of digits to print.}
    \item{width}{The width of printing.}
    \item{what}{Specifies what to print. Default option is 'obj' is the print function for mnlogit objects. Option 'eststat' prints etimation stats and option 'mdsize' prints model size information.}
    \item{probability}{If TRUE predict output the probability matrix, otherwise the chocice with the highest probability for each observation is returned.}
    \item{returnData}{If \code{TRUE} a data attribute is added to the returned object.}
    \item{order}{If \code{TRUE} coefficients are ordered by variable name.}
    \item{as.list}{Returns estimated model coefficients grouped by variable type.}
}   

\value{
    An object of class \code{mnlogit}, with elements:
    \item{coefficients}{the named vector of coefficients.}
    \item{logLik}{the value of the log-likelihood function at exit.}
    \item{gradient}{the gradient of the log-likelihood function at exit.}
    \item{hessian}{the Hessian of the log-likelihood function at exit.}
    \item{est.stat}{Newton Raphson stats.}
    \item{fitted.values}{Estimated probabilities of the alternative selected in each observation.}
    \item{probabilities}{the probability matrix: \code{(i,j)} entry denotes the probability of the \code{jth} alternative being chosen in the \code{ith} observation.}
    \item{residuals}{The residual. Has attribute \code{outcome} which is the probability of not choosing the selected alternative.}
    \item{df}{The number of estimated coefficients in the model.}
    \item{AIC}{The AIC value of the fitted model.}
    \item{choices}{The vector of alternatives's names.}
    \item{model.size}{Information about number of parameters in model.}
    \item{ordered.coeff}{Vector of coefficients ordered by variable name.}
    \item{model}{The \code{data.frame} used in model estimation.}
    \item{freq}{The relative frequency of each choice in input data.}
    \item{formula}{The \code{formula} specifying the  model.}
    \item{call}{The \code{mnlogit} function call that user made, }
}

\note{
1. The data must be in the 'long' format. This means that for each observation there must be as many rows as there are alternatives (which should be grouped together).

2. The formula should  be specified  in the format: responseVar ~ choice specific variables with generic coefficients | individual specific variables | choice specific variables with choice specific coefficients. These are the 3 available variable types.

3. Any type of variables may be omitted. To omit use "1" as a placeholder.

4. An alternative specific intercept is included by default in the estimation. To omit it, use a '-1' or '0' anywhere in the formula.
}

\references{
  Asad Hasan, Zhiyu Wang, Alireza S. Mahani (2016).Fast Estimation of Multinomial Logit Models: R Package mnlogit. Journal of Statistical Software, 75(3), 1-24. doi:10.18637/jss.v075.i03
 
  Croissant, Yves. 
    Estimation of multinomial logit models in R: The mlogit Packages.
    \url{https://cran.r-project.org/package=mlogit}

  Train, K. (2004). Discrete Choice Methods with Simulation,
  Cambridge University Press.
}

\author{Asad Hasan, Wang Zhiyu, Alireza S. Mahani}
\keyword{mnlogit, logistic, classification, multinomial, mlogit, parallel}

\examples{

  library(mnlogit)
  data(Fish, package = "mnlogit")
  fm <- formula(mode ~ price | income | catch)

  fit <- mnlogit(fm, Fish, ncores = 2)

 \dontrun{
   fit <- mnlogit(fm, Fish, choiceVar="alt", ncores = 2) # same effect as previous
   summary(fit)
   print(fit)
   predict(fit)
   print(fit, what = "eststat")
   print(fit, what = "modsize")
  
   # Formula examples (see also Note)
   fm <- formula(mode ~ 1 | income)    # Only type-2 with intercept
   fm <- formula(mode ~ price - 1)     # Only type-1, no intercept
   fm <- formula(mode ~ 1 | 1 | catch) # Only type-3, including intercept 
  }
}
