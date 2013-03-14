###############################################
#                                             #
#	Methods for mnlogit objects           # 
#                                             #
###############################################

# Extract the coefficients of the mnlogit object
coef.mnlogit <- function(object, as.list = FALSE, ...)
{
    if (!as.list) return(object$ordered.coeff)

    vec <- object$coeff
    p = object$model.size$p
    f = object$model.size$f
    d = object$model.size$d
    K = object$model.size$K
    np = object$model.size$nparams
     
    coeffList <- list(genericCoeff = NULL, indvCoeff = NULL, chSpCoeff = NULL)
    coeffList$indvCoeff <- if (p) vec[1:(K-1)*p] else NULL
    coeffList$chSpCoeff <- if (f) vec[((K-1)*p+1):((K-1)*p+K*f)] else NULL
    coeffList$genericCoeff <- if (d) vec[((K-1)*p+K*f+1):np] else NULL
    return(coeffList)
}


residuals.mnlogit <- function(object,...){
    object$residuals
}

# Extract the covariance matrix of the mnlogit object
vcov.mnlogit <- function(object, ...)
{
    result <- solve(object$hessian, tol = 1e-25)
    return(result)
}

# Extract the logliklihood info of the mnlogit object
logLik.mnlogit <- function(object, ...){
    object$logLik
}

# Print estimation statistics 
print.est.stats <- function(x, ...)
{
    cat("Maximum likelihood estimation using Newton-Raphson iterations.")
    cat(paste0("\n  Number of iterations: ", x$niters))
    cat(paste0("\n  Number of linesearch iterations: ", x$LSniters))
    cat(paste0("\nAt termination: "))
    cat(paste0("\n  Gradient norm = ", round(x$gradNorm, 8)))
    cat(paste0("\n  Diff between last 2 loglik values = ", round(x$funcDiff, 8)))
    cat(paste0("\n  Stopping reason: ", x$stopCond))
    cat(paste0("\nTotal estimation time (sec): ", round(x$totalMins*60, 3)))
    cat(paste0("\nTime for Hessian calculations (sec): ",
               round(x$hessMins*60, 3), " using ", x$ncores, " processors.\n"))
}

# Print model parameters
print.model.size <- function(x, ...)
{
    cat(paste0("Number of observations in training data = ", x$N))
    cat(paste0("\nNumber of alternatives = ", x$K))
    if (x$intercept) cat("\nIntercept turned: ON.")
    else cat("\nIntercept turned: OFF.")
    cat(paste0("\nNumber of parameters in model = ", x$nparams))
    cat(paste0("\n  # individual specific variables = ", x$p))
    cat(paste0("\n  # choice specific coeff variables = ", x$f))
    cat(paste0("\n  # generic coeff variables = ", x$d, "\n"))
}

# Print function
print.mnlogit <- function(x, digits = max(3, getOption("digits") - 2),
                          width = getOption("width"), ...)
{
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    if (length(coef(x))){
        cat("Coefficients:\n")
        print.default(format(coef(x), digits = digits), print.gap = 2 ,
                      quote = FALSE) 
    }
    else cat("No coefficients\n")
    cat("\n")
    invisible(x)
}

# Summary function
summary.mnlogit <- function(object,...)
{
    b <- coef(object)
    std.err2 <- sqrt(diag(vcov(object)))
    std.err <- b
    std.err[names(std.err2)] <- std.err2

    z <- b / std.err
    p <- 2 * ( 1 - pnorm(abs(z)))
    CoefTable <- cbind(b, std.err, z, p)
    colnames(CoefTable) <- c("Estimate", "Std.Error", "t-value", "Pr(>|t|)")
    object$CoefTable <- CoefTable
    class(object) <- c("summary.mnlogit", "mnlogit")
    return(object)
}

print.summary.mnlogit <- function(x, digits = max(3, getOption("digits") -2),
                                  width = getOption("width"), ...)
{
    cat("\nCall:\n")
    print(x$call)
    cat("\n")
    cat("Frequencies of alternatives in input data:\n")
    print(prop.table(x$freq), digits = digits)
    cat("\n")
    print(x$model.size)
    cat("\n")
    print(x$est.stat)
    cat("\nCoefficients : \n")
    printCoefmat(x$CoefTable, digits = digits)
    cat("\n")
    cat(paste0("Log-Likelihood: ", signif(x$logLik, digits), ", df = ",
               x$model.size$nparams, "\n"))
    cat(paste("AIC: ", x$AIC, "\n"))
    cat("\n")
}
