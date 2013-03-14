###############################################################################
#                        Package: mnlogit                                     #
#                                                                             #
# Multinomial logit, maximum likelihood estimation by Newton-Raphson method   #
#                                                                             #
# Implementors: Wang Zhiyu, Asad Hasan                                        # 
#               Scientific Computing Group, Sentrana Inc.                     #
###############################################################################

###############################################################################
#                 Main user function                                           
# Args:  									
#   formula     - same format as mlogit. See help(formula)
#   data        - input data (as a data.frame object) in "long" format
#   choiceVar   - the data column containing alternative names
#   maxiter     - maximum number of Newton-Raphson iterations to run 
#   ftol        - function tolerance. 
#                 Difference of two consecutive function evaluation
#                 Criteria of terminating Newton's Iterative process
#   gtol        - gradient norm tolerance.
#   weights     - an optional vector of positive frequency weights.
#   ncores      - number of processors allowed to use
#   na.rm       - if FALSE then stop(), else remove rows of data with NA
#   print.level - increase from 0 to progressively print more runing info
#   linDepTol   - tolerance with which linear dep among cols is detected
#   ...         - currently unused 
#        
# Output: 
#   mnlogit object
#        
# Note:    
#   'ftol', 'gtol' & 'maxiter' specify Newton-Raphson termination criteria 
###############################################################################
mnlogit <- function(formula, data, choiceVar, maxiter = 50, ftol = 1e-6,
                    gtol = 1e-6, weights = NULL, ncores = 1, na.rm = TRUE,
                    print.level = 0, linDepTol = 1e-6, ...)
{
    startTime <- proc.time()[3]
    initcall <- match.call()    # Store original function call
    predict <- NULL             # Disused 

    # Basic parameter checking
    if (!is.data.frame(data))
        stop("'data' must be a data.frame in 'long' format")
    if (ncores < 1) {
        ncores <- 1
        warning("Setting ncores equal to: 1")
    }

    # Extract various types of variables from formula
    formula <- parseFormula(formula) 
    response <- attr(formula, "response")     # response variable
    interceptOn <- attr(formula, "Intercept") # check if intercept is in model
    csvChVar <- attr(formula, "csvChCoeff") 
    indspVar <- attr(formula, "indSpVar") 
    csvGenVar <- attr(formula, "csvGenCoeff") 
    covariates <- c(csvChVar, indspVar, csvGenVar)    
    varNames <- attr(formula, "varNames")    

    if (is.null(covariates) && !interceptOn) 
        stop("Error! Predictor variable(s) must be specified")
    if (is.null(response)) 
        stop("Error! Alternative variable must be specified")
     
    # Determine relevant parameters
    choice.set <- unique(data[[choiceVar]])
    K <- length(choice.set) # number of choices
    if (nrow(data) %% K)
        stop("Mismatch between number of rows in data and number of choices.")
    N <- nrow(data)/K       # number of individuals

    # Check if weights is OK 
    if (!is.null(weights) && length(weights) != N)
      stop("Length of 'weights' arg must match number of observations in data.")
    if (!is.null(weights) && !all(weights > 0))
      stop("All entries in 'weights' must be strictly positive.")      

    # Work with only the columns appearing in formula
    data <- data[c(varNames, choiceVar)]

    # Handle NA; Find out row numbers with atleast one NA
    na.rows <- c()
    for (col in 1:ncol(data))
        na.rows <- union(na.rows, which(is.na(data[[col]])))
    Ndropped <- 0
    if (length(na.rows) > 0) {
        if (!na.rm)
            stop("NA present in input data.frame with na.rm = FALSE.")
        # Mark rows with NA for deletion
        keepRows <- rep(TRUE, nrow(data))
        keepRows[na.rows] <- FALSE 
        # Starting with 1st, mark rows for deletion in groups of K
        for (i in 1:N) {
            if (!all(keepRows[((i-1)*K + 1):(i*K)]))
                keepRows[((i-1)*K + 1):(i*K)] <- FALSE
        }
        data <- data[keepRows, , drop=FALSE]
        # Drop weights corresponding to dropped rows 
        if (!is.null(weights)) {
            weights <- weights[keepRows[seq(1, N * K, K)]]
        }
        N <- nrow(data)/K
        Ndropped <- (length(keepRows) - sum(keepRows))/K
    }
    if (print.level && Ndropped > 0) 
      cat(paste("Num of dropped observations (due to NA)  =", Ndropped, "\n"))
      
    # Rearrange the input data.frame object
    # Sort according to choices: data for an atlernative should be contiguous
    data <- data[order(data[[choiceVar]]), ]

    # Obtain response vector as a vector of 0,1
    respVec <- as.numeric(data[[attr(formula, "response")]])
    min.respVec <- min(respVec)
    spread <- max(respVec) - min.respVec
    if (spread != 1) {
        stop(paste("Response variable", attr(formula, "response"), 
                   "must be a factor with exactly 2 levels."))
    }
    respVec <- respVec - min.respVec
    freq.choices <- colSums(matrix(respVec, nrow = N, ncol = K))/N
    loFreq <- min(freq.choices)
    loChoice <- choice.set[which(loFreq == freq.choices)]
    names(freq.choices) <- choice.set
    if (loFreq < 1e-7) {
        cat("Frequencies of alternatives in input data:\n")
        print(prop.table(freq.choices), digits = 4)
        stop(paste("Frequency, in response, of choice:", loChoice, "< 1e-7."))
    }

    # Form design matrices 
    formDesignMat <- function(varVec = NULL, includeIntercept = TRUE)
    {
        if (is.null(varVec) && !includeIntercept) return(NULL) 
        fm <- paste(attr(formula, "response"), "~")
        if (!is.null(varVec))
            fm <- paste(fm, paste(varVec, collapse = "+"))
        if (!includeIntercept) fm <- paste(fm, "-1 ")
        else fm <- paste(fm, "+1 ")
        modMat <- model.matrix(as.formula(fm), data)
    } 
    X <- formDesignMat(varVec = attr(formula, "indSpVar"), 
                       includeIntercept = attr(formula, "Intercept"))
    X <- if (!is.null(X)) X[1:N, , drop=FALSE]   # Matrix of ind sp vars
    Y <- formDesignMat(varVec = attr(formula, "csvChCoeff"), 
                       includeIntercept = FALSE)
    Z <- formDesignMat(varVec = attr(formula, "csvGenCoeff"), 
                       includeIntercept = FALSE)

    # Detect bad columns (these are linearly dependent on other columns)
    badColsList <- list("indSpVar"=NULL, "csvChCoeff"=NULL, "csvGenCoeff"=NULL)
    getNullSpaceCols <- function(mat, tol = 1e-7)
    {
        if (is.null(mat)) return(NULL)
        if (ncol(mat)==1) return(NULL)
        qrdecomp <- qr(mat, tol = tol)
        rank <- qrdecomp$rank
        if (rank == ncol(mat)) return(NULL)
        nullSpCols <- qrdecomp$pivot[(rank + 1):ncol(mat)]
        return(nullSpCols)
    }
    badColsList$indSpVar <- getNullSpaceCols(X, tol = linDepTol)    
    for (i in 1:K) {
        init <- (i-1)*N + 1
        fin <- i*N
        badColsList$csvChCoeff <- union(badColsList$csvChCoeff,
            getNullSpaceCols(Y[init:fin, , drop=FALSE], tol = linDepTol))
    }
    badColsList$csvGenCoeff <- getNullSpaceCols(Z, tol = linDepTol)

    # Get names of variables to be dropped from estimation
    badVarsList <- list()
    badVarsList$indSpVar <- colnames(X[, badColsList$indSpVar, drop=FALSE])
    badVarsList$csvChCoeff <- colnames(Y[, badColsList$csvChCoeff, drop=FALSE])
    badVarsList$csvGenCoeff <- colnames(Z[,badColsList$csvGenCoeff,drop=FALSE])
    badCoeffNames <- makeCoeffNames(badVarsList, choice.set)

    # Eliminate linearly dependent columns
    if (!is.null(X) && is.null(predict))
      X <- X[ , setdiff(1:ncol(X), badColsList$indSpVar), drop=FALSE]
    if (!is.null(Y) && is.null(predict))
      Y <- Y[ , setdiff(1:ncol(Y), badColsList$csvChCoeff), drop=FALSE]
    if (!is.null(Z) && is.null(predict))
      Z <- Z[ , setdiff(1:ncol(Z), badColsList$csvGenCoeff), drop=FALSE]
 
    # Get names of variables 
    varNamesList <- list()
    varNamesList$indSpVar <- colnames(X)
    varNamesList$csvChCoeff <- colnames(Y)
    varNamesList$csvGenCoeff <- colnames(Z)
    coeffNames <- makeCoeffNames(varNamesList, choice.set)
    
    # Do the subtraction: Z_ik - Zi0 (for Generic coefficients data)
    ### NOTE: Base choice (with respect to normalization) is fixed here
    ###       Base choice is the FIRST alternative
    baseChoiceName <- choice.set[1]
    if(!is.null(Z)) { 
        for (ch_k in 2:K) {
            Z[((ch_k - 1)*N + 1):(ch_k*N), ] <-
              Z[((ch_k-1)*N+1):(ch_k*N), , drop=FALSE] - Z[1:N, , drop=FALSE]
        }
    }
    # Drop rows for base alternative
    Z <- Z[(N + 1):(K*N), , drop=FALSE]
    respVec <- respVec[(N + 1):(K*N)]
 
    t1 <- proc.time()[3]    # Time at end of pre-processing

    # Predict probability matrix and return to caller
    if (!is.null(predict)) {
        if (!is.null(badCoeffNames)) 
            stop("Collinear columns in data. Prediction cancelled.")
        predict <- eval.parent(predict)
        probmat <- newtonRaphson(respVec, X, Y, Z, K, maxiter, gtol, ftol,
                                 ncores, 0, NULL, predict)
        colnames(probmat) <- choice.set
        return(probmat)
    }

    gc()  # Invoke garbage collector at end of pre-processing 
    if (print.level > 1) {
      cat(paste0("Base alternative is: ", baseChoiceName))
      cat(paste0("\nPreprocessing data for estimation took ", 
                  round(t1 - startTime, 3), " sec.\n"))
    } 

    # Solve MLE using Newton-Raphson
    result <- newtonRaphson(respVec, X, Y, Z, K, maxiter, gtol, ftol, ncores,
                            print.level, coeffNames, weights)    
    # Post-processing
    colnames(result$hessMat) <- coeffNames 
    rownames(result$hessMat) <- coeffNames
    names(result$grad)   <- coeffNames
    # Reorder coeff so that those for a choice are together 
    od <- reordering(varNamesList, choice.set)
    coeffNames <- makeCoeffNames(varNamesList, choice.set)
    coefficients <- c(result$coeff, if (is.null(badCoeffNames)) NULL
                                    else rep(NA, length(badCoeffNames)))
    names(coefficients) <- c(coeffNames,
        badCoeffNames[reordering(badVarsList, choice.set)])
    reordered_coeff <- c(result$coeff[od], if (is.null(badCoeffNames)) NULL
                                           else rep(NA, length(badCoeffNames)))
    names(reordered_coeff) <- c(coeffNames[od],
        badCoeffNames[reordering(badVarsList, choice.set)])
    colnames(result$probability) <- choice.set[-1]
    AIC <- 2*(result$model.size$nparams - log(abs(result$loglikelihood)))
    result$model.size$intercept <- interceptOn

    fit <- structure(list(
              coeff = coefficients,
              probabilities = result$probability,
              residuals = result$residual,
              logLik = -result$loglikelihood,
              df = result$model.size$nparams,
              gradient = -result$grad,
              hessian = result$hessMat,
              AIC = AIC,
              formula = formula,
              data = data,
              choices = choice.set,
              freq = freq.choices,
              model.size = result$model.size,
              est.stats = result$est.stats,
              ordered.coeff = reordered_coeff,
              call = initcall
           ), class = "mnlogit")

    if (print.level)
      cat(paste0("\nTotal time spent in mnlogit = ",
                 round(proc.time()[3] - startTime, 3), " seconds.\n"))
    return(fit)
}

# Makes names of model coefficients
makeCoeffNames <- function (varNames, choices)
{
    if (length(varNames) == 0) return(NULL)
    choices <- as.vector(choices)
    coeffName <- c(outer(varNames$indSpVar, choices[-1], paste, sep=":"),
                   outer(varNames$csvChCoeff, choices, paste, sep=":"),
                   varNames$csvGenCoeff)
}

# Generate a re-ordering of coeff names (group choices together)
reordering <- function(varList, choices)
{
    if (length(varList) == 0) return(NULL)
    K <- length(as.vector(choices))
    p <- length(as.vector(varList$indSpVar))
    f <- length(as.vector(varList$csvChCoeff))
    d <- length(as.vector(varList$csvGenCoeff))
    orig <-  c(if (p > 0) rep(1:p, K-1) else NULL,
               if (f > 0) rep((p+1):(p+f), K) else NULL,
               if (d > 0) (p+f+1):(p+f+d) else NULL)
    order(orig) 
}
