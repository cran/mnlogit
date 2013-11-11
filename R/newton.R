#############################################################################
#     Optimizer implementing the Newton-Raphson Method (with linesearch)
#
# Args:
#   reponse - vector giving the response (for training data) 
#   X - Design Matrix for Individual Specific Variable
#   Y - Design for choice attributes with individual variation & ch sp coeff
#   Z - Design for choice attributes and generic coeff
#   K - number of  choices
#   Dimension specification:
#    N - number of observations
#    p - # individual specific variables
#    f - # choice sp variables choice coeff
#    d - # choice specific variables gen coeff
#    dim(X) - N * p
#    dim(Y) - N*K * f (f - cols)
#    dim(Z) - N*(K-1) * d (rows corresponding to base choice are absent)
#    NOTE: Data for each choice is contiguous in matrices Y & Z.
#   Criteria of terminating Newton's Iterative process
#     maxiter - maximum number of Newton-Raphson iterations to run 
#     ftol    - function tolerance. 
#               Difference of two consecutive function evaluation
#     gtol    - gradient norm tolerance.
#   ncores  - number of processors allowed to use
#   print.level - increase from 0 to progressively print more running info
#   coeff.names - names of coefficients, vector of length = nparams 
#   predict - if NOT NULL, then it must vecotr of coeff of length nparam
#             and is used to compute predicted probability matrix
#
# Output: 
#  A list with various entries (see end of function) 
#
# Notes:
#   A simple linesearch procedure to forbid the Newton step from increasing
#   loglikelihood is also implemented.
#############################################################################
newtonRaphson <- function(response, X, Y, Z, K, maxiter, gtol, ftol, 
                          ncores, print.level, coeff.names, predict = NULL)
{
    initTime <- proc.time()[3]
    # Determine problem parameters
    N <- length(response)/(K - 1)       # num of observations
    p <- ifelse(is.null(X), 0, ncol(X)) # num of ind sp variables
    f <- ifelse(is.null(Y), 0, ncol(Y)) # num of ind sp variables
    d <- ifelse(is.null(Z), 0, ncol(Z)) # num of ind sp variables
    nparams <- (K - 1)*p + K*f + d # total number of coeffs
    size <- structure(list(
              "N" = N,
              "K" = K, 
              "p" = p, 
              "f" = f, 
              "d" = d, 
              "nparams" = nparams),
              class = "model.size")
    if (!is.null(predict)) {
        if (length(predict) != size$nparams)
            stop("Arg predict of incorrect length or not a vector.")
        return(likelihood(response, X, Y, Z, size, predict, ncores,
                          hess=FALSE, predict = TRUE))
    }

    # Initialize 'guess' for NR iteration
    coeffVec <- rep(0, size$nparams)
    funcTime <- gradTime <- hessTime <- lineSearch <- solveTime <- 0.0
    # Compute log-likelihood, gradient, Hessian at initial guess
    fEval <- likelihood(response, X, Y, Z, size, coeffVec, ncores)
    funcTime <- funcTime + attr(fEval, "funcTime")
    gradTime <- gradTime + attr(fEval, "gradTime")
    hessTime <- hessTime + attr(fEval, "hessTime")
    loglik <- fEval[[1]]
    gradient <- attr(fEval, "gradient")
    hessian <- attr(fEval, "hessian")
    gradNorm <- as.numeric(sqrt(crossprod(gradient, gradient)))
    lineSearchIters <- 0
    stop.code <- "null"

    # Newton-Raphson Iterations
    for (iternum in 1:maxiter) {
        oldLogLik <- loglik
        oldCoeffVec <- coeffVec

        if (print.level) {
          cat("====================================================")
          cat(paste0("\nNewton Raphson iteration # ", iternum))
          cat(paste0("\n  loglikelihood = ", loglik))
          cat(paste0("\n  gradient 2-norm = ", gradNorm))
          cat("\n  Approx Hessian condition number = ")
          cat(paste0(round(1.0/rcond(hessian), 2), "\n"))
        }        
        if (print.level > 1) {
            names(gradient) <- names(coeffVec) <- coeff.names
            coefTable <- rbind(coeffVec, gradient)
            rownames(coefTable) <- c("coef", "grad")
            print(coefTable, digits = 3)
            if (print.level > 2) {
                cat("\nPrinting the Hessian matrix.\n")
                colnames(hessian) <- rownames(hessian) <- coeff.names 
                print(hessian, digits = 3)
            }
        } 
        
        # Find NR update vector by solving a linear system
        t0 <- proc.time()[3] 
        dir <- -1*as.vector(solve(hessian, gradient, tol = 1e-24)) 
        solveTime <- solveTime + proc.time()[3] - t0 
    
        # Do the linesearch trying first to take the full Newton step
        t1 <- proc.time()[3] 
        alpha <- 1.0
        niter <- 0
        newloglik <- NULL
        while(alpha >= 1e-4) {
            niter <- niter + 1
            coeffVec <- oldCoeffVec + alpha*dir
            newloglik <- likelihood(response, X, Y, Z, size, coeffVec, ncores,
                                    hess = FALSE)
            funcTime <- funcTime + attr(newloglik, "funcTime")
            gradTime <- gradTime + attr(newloglik, "gradTime")
            hessTime <- hessTime + attr(newloglik, "hessTime")
            if (newloglik[[1]] < oldLogLik) break
            else alpha <- alpha/2.0
        } 
        t2 <- proc.time()[3]
        lineSearchIters <- lineSearchIters + niter
        lineSearch <- lineSearch + t2 - t1

        # Compute log-likelihood, gradient, Hessian
        fEval <- likelihood(response, X, Y, Z, size, coeffVec, ncores,
                            loglikObj = newloglik)
        funcTime <- funcTime + attr(fEval, "funcTime")
        gradTime <- gradTime + attr(fEval, "gradTime")
        hessTime <- hessTime + attr(fEval, "hessTime")
        loglik <- fEval[[1]]
        gradient <- attr(fEval, "gradient")
        hessian <- attr(fEval, "hessian")
        gradNorm <- as.numeric(sqrt(crossprod(gradient, gradient)))
         
        # First Success Criteria: function values converge
        loglikDiff <- abs(loglik - oldLogLik)    
        if (loglikDiff < ftol) {
          stop.code <- paste0("Succesive loglik difference < ftol (", ftol, ").")
          break
        }
        # Second Success Criteria: gradient reduces below gtol 
        if (gradNorm < gtol) {
          stop.code <- paste0("Gradient 2-norm < gtol (", gtol, ").")
          break
        }
    }  # end Newton-Raphson loop
    if (stop.code != "null" && iternum == maxiter)
          stop.code <- paste("Number of iterations:", iternum, "== maxiters.")
    if (print.level) {
        cat("====================================================")
        cat(paste("\nTermination reason:", stop.code, "\n"))
    }
     
    # Newton Raphson Stats
    stats <- structure(list(
               funcDiff = loglikDiff,
               gradNorm = gradNorm,
               niters = iternum,
               LSniters = lineSearchIters,
               stopCond = stop.code,
               totalMins = (proc.time()[3] - initTime)/60,
               hessMins = hessTime/60,
               ncores = ncores
            ), class = "est.stats")

    # Residual - probability of NOT making the choice that was really made
    responseMat <- matrix(response, nrow=size$N, ncol=(size$K-1))
    baseResp <- rep(1, size$N) - rowSums(responseMat)
    probMat <- attr(fEval, "probMat")
    baseProbVec <- 1 - rowSums(probMat)
    # Pch - prob of choice that was actually made
    Pch <- matrix(c(as.vector(ifelse(baseResp > 0, baseProbVec, NA)),
                    as.vector(ifelse(responseMat > 0, probMat, NA))),
                    nrow = size$K, ncol = size$N, byrow=TRUE)
    Pch <- Pch[!is.na(Pch)]
    Pch <- 1 - Pch   # the residual (oen for each observation) 

    result <- list(coeff = coeffVec, loglikelihood = loglik, grad  = gradient,
                   hessMat = hessian, probability = probMat, residual = Pch,
                   model.size = size, est.stats = stats)
    return (result)
}
