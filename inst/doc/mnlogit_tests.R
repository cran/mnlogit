############################################################################
### Sample code for results of Table 2 of the paper
############################################################################

library(mnlogit)
library(mlogit)
source("simChoiceModel.R") #  Has makeModel() to generate simulated data

# Type 'X' problems
numChoices <- 10
data <- makeModel('X', K=numChoices) # generate data 
# Default args set: p = 50 variables, N = K * p * 20 observations

# Make formula 
vars <- paste("X", 1:50, sep="", collapse=" + ")
fm <- formula(paste("response ~ 1|", vars, " - 1 | 1"))

# Run mnlogit
system.time(fit.mnlogit <- mnlogit(fm, data, "choices"))  # runs on 1 proc

# Run mlogit
mdat <- mlogit.data(data[order(data$indivID), ], "response", shape="long",
    alt.var="choices")
system.time(fit.mlogit <- mlogit(fm, mdat))   # Newton-Raphson
system.time(fit.mlogit <- mlogit(fm, mdat, method='bfgs'))

# Type 'Y, 'Z' & 'YZ' problems
data <- makeModel('Y', K=numChoices)  # generate data 
# Default args set: p = 50 variables, N = K * p * 20 observations

# Make formula for type 'Y' problems 
vars <- paste("X", 1:50, sep="", collapse=" + ")
fm <- formula(paste("response ~ 1| - 1 | ", vars))

# Run mnlogit
system.time(fit.mnlogit <- mnlogit(fm, data, "choices"))  # runs on 1 proc

# Run mlogit
mdat <- mlogit.data(data[order(data$indivID), ], "response", shape="long",
    alt.var="choices")
system.time(fit.mlogit <- mlogit(fm, mdat))   # Newton-Raphson
system.time(fit.mlogit <- mlogit(fm, mdat, method='bfgs'))

# Formula for type 'Z' problems
vars <- paste("X", 1:50, sep="", collapse=" + ")
fm <- formula(paste("response ~ ", vars, "| - 1 | 1"))
# Code for running mnlogit and mlogit is the same as for type 'X' & 'Z'

# Formula for type 'YZ' problems
# 5 variables of type 'Z' and 45 variables of type 'Y'
vars <- paste("X", 1:45, sep="", collapse=" + ")
fm <- formula(paste("response ~ X46 + X47 + X48 + X49 + X50| - 1 | ", vars))
# Code for running mnlogit and mlogit is the same as for type 'X' & 'Z'

############################################################################
### Sample code for results of Table 3 of the paper (parallel execution)
############################################################################

library(mnlogit)
source('simChoiceModel.R')

# Type 'X' problems
numChoices <- 20
data <- makeModel('X', K=numChoices)  # generate data 
# Default args set: p = 50 variables, N = K * p * 20 observations

# Make formula 
vars <- paste("X", 1:50, sep="", collapse=" + ")
fm <- formula(paste("response ~ 1|", vars, " - 1 | 1"))

system.time(fit.mnlogit <- mnlogit(fm, data, "choices", ncores = 1)) 
system.time(fit.mnlogit <- mnlogit(fm, data, "choices", ncores = 2)) 
system.time(fit.mnlogit <- mnlogit(fm, data, "choices", ncores = 4)) 

############################################################################
### Code from Appendix C of the paper 
############################################################################

library(mnlogit)
source('simChoiceModel.R')

# Generate simulated data
data <- makeModel('X', K=5)
K = length(unique(data$choices))
N = nrow(data)/K
p = ncol(data) - 3
np = (K - 1) * p
cat(paste0("Number of choices in simulated data = K = ", K,
    ".\nNumber of observations in si    mulated data = N = ", N,
    ".\nNumber of variables = p = ", p,
    ".\nNumber of model parameters = (K - 1) * p = ", (K-1)*p, "."))

# Make formula for mnlogit and mlogit
vars <- paste("X", 1:50, sep="", collapse=" + ")
fm <- formula(paste("response ~ 1|", vars, " - 1 | 1"))

# Run mnlogit
system.time(fit.mnlogit <- mnlogit(fm, data, "choices"))  # runs on 1 proc

# Run mlogit
library(mlogit) 
mdat <- mlogit.data(data[order(data$indivID), ], "response", shape="long",
    alt.var="choices")
system.time(fit.mlogit <- mlogit(fm, mdat))   # Newton-Raphson
system.time(fit.mlogit <- mlogit(fm, mdat, method='bfgs'))

# Run nnet
library(nnet)
ndat <- data[which(data$response > 0), ]
ff <- paste("choices ~", vars, "- 1")   # formula for nnet
system.time(fit.nnet <- multinom(ff, ndat, reltol=1e-12))

# Run VGAM
library(VGAM)
stop.vglm <- vglm.control(epsilon = 1e-6)
system.time(fit.vglm <- vglm(ff, data=ndat, multinomial, control=stop.vglm))
#############################################################################
