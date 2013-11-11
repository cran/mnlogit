##############################################################################
#  Formula Parsing Routine  
#  
#  A formula for mnlogit is of type:
#
#  response ~   choice specific variables with generic coefficients
#             | individual specific variables with choice specific coefficients
#             | choice specific variables with individual variation and choice
#               specific coefficiients                    
#                           
###############################################################################
# Argument: 
#   f - formula object or anything coercable to a formula object
#
# Output: 
#  Input formula object with attributes:
#    varNames: name of all variables in formula (inlcuding response)
#    response: name of response variable
#    Intercept: logical variable indicating whether an additional 
#              'intercept' parameter (which is choice sp) is estimated.
#    csvGenCoeff: vector of choice specific variables with generic coeff
#    indSpVar: vector of indivicual specific variables 
#    csvChCoeff: vector of choice specific variables with choice sp coeff
#
# Note: 
#    1. Presence of '-1' or '0' indicates intercept is turned of 
#    2. To NOT include any variable sometype, use a '1' as a placeholder
###############################################################################
parseFormula <- function(f)
{     
    # Coerce to formula type
    f <- as.formula(f)
    call <- f
    attr(call, "varNames") <- all.vars(f) # get all variable names

    # Split the formula into LHS & RHS, separator: ~
    f <- as.character(f) 
    if (f[1] != "~" || length(f) < 3)
        stop("Not a valid formula")
    response <- f[2]

    # Start parsing the formula string
    args <- unlist(strsplit(f[3], split= "|", fixed=TRUE))
    args <- gsub(" ","",args) # Delete all SiNGLE SPACES

    interceptON  <- TRUE
    if (args[1] != "") { 
        vars <- getTerms(paste(response, "~", args[1])) 
        attr(call, "csvGenCoeff")  <- attr(vars, "covariates")
        interceptON <- (interceptON && attr(vars, "intercept"))
    } 
    if (length(args) > 1 && args[2] != "") { 
        vars <- getTerms(paste(response, "~", args[2])) 
        attr(call, "indSpVar")  <- attr(vars, "covariates")
        interceptON <- (interceptON && attr(vars, "intercept"))
    } 
    if (length(args) > 2 && args[3] != "") { 
        vars <- getTerms(paste(response, "~", args[3])) 
        attr(call, "csvChCoeff")  <- attr(vars, "covariates")
        interceptON <- (interceptON && attr(vars, "intercept"))
    } 
    attr(call, "Intercept") <- interceptON 
    attr(call, "response")  <- response 
    return(call)
}

getTerms <- function(formulaStr)
{
    fm <- as.formula(formulaStr)
    Terms <- terms(fm)
    attr(fm, "covariates") <- if (length(attr(Terms, "term.labels")))
                                  attr(Terms, "term.labels") else NULL
    attr(fm, "intercept") <- attr(Terms, "intercept")
    return(fm)
}
