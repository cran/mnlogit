\name{mnlogit}
\alias{mnlogit}
\alias{print.mnlogit}
\alias{summary.mnlogit}
\alias{predict.mnlogit}
\alias{coef.mnlogit}
\alias{print.est.stats}
\alias{print.model.size}
\alias{print.summary.mnlogit}

\title{Efficient estimation of multinomial logit model}

\description{
Time and memory efficient estimation of multinomial logit models using maximum likelihood method and targeted towards large scale multiclass classification problems in econometrics and machine learning. 
Numerical optimization performed by Newton-Raphson method using an optimized, parallel C++ library to achieve fast computation of Hessian matrices.
The user interface closely related to the CRAN package \pkg{mlogit}.
}

\usage{
mnlogit(formula, data, choiceVar, maxiter = 25, ftol = 1e-6,
        gtol = 1e-6, ncores = 1, na.rm = TRUE, print.level = 0, 
         linDepTol = 1e-6, ...)
\method{predict}{mnlogit}(object, newdata = NULL, probability = FALSE, ...)
\method{coef}{mnlogit}(object, as.list = FALSE, ...)
}
\arguments{
    \item{formula}{formula object or string specifying the model to be estimated (see Note).}
    \item{data}{A data.frame object with data organized in the 'long' format (see Note).}
    \item{choiceVar}{A string naming the column in 'data' which has the list of choices.}
    \item{maxiter}{An integer indicating maximum number of Newton's iterations,}
    \item{gtol}{A real number indicating tolerance on 2-norm of the gradient.}
    \item{ftol}{A real number indicating tolerance on the difference of two subsequent loglikelihood values.}
    \item{ncores}{An integer indicating number of processors allowed for Hessian calculations.}
    \item{na.rm}{a logical variable which indicates whether rows of the data frame containing NAs will be removed.}
    \item{print.level}{An integer which controls the amount of information to be printed during execution.}
    \item{linDepTol}{Tolerance for detecting linear dependence between columns in input data. Dependent columns are removed from the estimation.}
    \item{...}{Used in mnlogit for prediction - not meant for external usage.}
    \item{object}{A fitted mnlogit object.}
    \item{newdata}{A data.frame object to used for prediction.}
    \item{probability}{If TRUE predict output the probability matrix, otherwise the chocice with the highest probability for each observation is returned.}
    \item{as.list}{Returns estimated model coefficients grouped by variable type.}
}   

\value{
    An object of class mnlogit, with elements:
    \item{coeff}{the named vector of coefficients.}
    \item{probabilities}{the probability matrix: (i,j) entry denotes the probability of the jth choice being choosen in the ith observation.}
    \item{residuals}{the named vector of residuals which is the probability of not choosing the alternative which was chosen.}
    \item{logLik}{the value of the log-likelihood function at exit.}
    \item{df}{the number of parameters in the model.}
    \item{gradient}{the gradient of the log-likelihood function at exit.}
    \item{hessian}{the Hessian of the log-likelihood function at exit.}
    \item{AIC}{the AIC value of the fitted model.}
    \item{formula}{the formula specifying the  model.}
    \item{data}{the data.frame used in model estimation.}
    \item{choices}{the vector of alternatives.}
    \item{freq}{the frequencies of alternatives.}
    \item{model.size}{Information about number of parameters in model.}
    \item{est.stat}{Newton Raphson stats.}
    \item{freq}{the frequency of each choice in input data.}
    \item{call}{the mnlogit function call that user made, }
}

\note{
1. The data must be in the 'long' format. This means that for each observation there must be as many rows as there are alternatives (which should be grouped together).

2. The formula should  be specified  in the format: responseVar ~ choice specific variables with generic coefficients | individual specific variables | choice specific variables with choice specific coefficients. These are the 3 available variable types.

3. Any type of variables may be omitted. To omit use "1" as a placeholder.

4. An alternative specific intercept is included by default in the estimation. To omit it, use a '-1' or '0' anywhere in the formula.
}

\references{
  Croissant, Yves. 
    \emph{Estimation of multinomial logit models in R: The mlogit Packages.}
    \url{http://cran.r-project.org/web/packages/mlogit/index.html}

  Train, K. (2004) \emph{Discrete Choice Methods with Simulation},
  Cambridge University Press.
}

\author{Wang Zhiyu, Asad Hasan}
\keyword{mnlogit, logistic, classification, multinomial, mlogit, parallel}

\examples{
  library(mnlogit)

  data(Fish, package = "mnlogit")

  fm <- formula(mode ~ price | income | catch)
  result <- mnlogit(fm, Fish, "alt", ncores = 2)
  predict(result)

 \dontrun{

  print(result)
  print(result$est.stats)
  print(result$model.size)
  summary(result)
  
  # Formula examples (see also Note)
  fm <- formula(mode ~ 1 | income) # Only type-2 with intercept
  fm <- formula(mode ~ price - 1) # Only type-1, no intercept
  fm <- formula(mode ~ 1 | 1 | catch) # Only type-3, including intercept 

  }
}
