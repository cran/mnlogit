### R code from vignette source 'mnlogit.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: mnlogit.Rnw:113-114
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: mnlogit.Rnw:171-172
###################################################
library(mnlogit)


###################################################
### code chunk number 3: mnlogit.Rnw:177-179
###################################################
data(Fish, package = 'mnlogit')
head(Fish, 8)


###################################################
### code chunk number 4: mnlogit.Rnw:274-275 (eval = FALSE)
###################################################
## fm <- formula(mode ~ price | income | catch)


###################################################
### code chunk number 5: mnlogit.Rnw:279-282 (eval = FALSE)
###################################################
## fm <- formula(mode ~ price | income - 1 | catch)
## fm <- formula(mode ~ price | income | catch - 1)
## fm <- formula(mode ~ 0 + price | income | catch)


###################################################
### code chunk number 6: mnlogit.Rnw:286-291 (eval = FALSE)
###################################################
## fm <- formula(mode ~ 1 | income | catch) 
## fm <- formula(mode ~ price | 1 | catch) 
## fm <- formula(mode ~ price | income | 1)
## fm <- formula(mode ~ price | 1 | 1)
## fm <- formula(mode ~ 1 | 1 | price + catch)


###################################################
### code chunk number 7: mnlogit.Rnw:295-298 (eval = FALSE)
###################################################
## fm <- formula(mode ~ price + catch | 1 | 1)
## fm <- formula(mode ~ price + catch | 1) 
## fm <- formula(mode ~ price + catch)


###################################################
### code chunk number 8: mnlogit.Rnw:308-309 (eval = FALSE)
###################################################
## ?mnlogit


###################################################
### code chunk number 9: mnlogit.Rnw:312-315 (eval = FALSE)
###################################################
## mnlogit(formula, data, choiceVar, maxiter = 50, ftol = 1e-6,
##         gtol = 1e-6, weights = NULL, ncores = 1, na.rm = TRUE,   
##         print.level = 0, linDepTol = 1e-6, ...)


###################################################
### code chunk number 10: mnlogit.Rnw:324-325
###################################################
fm <- formula(mode ~ price | income | catch)


###################################################
### code chunk number 11: mnlogit.Rnw:342-344
###################################################
fit <- mnlogit(fm, Fish, "alt", ncores=2)
class(fit)


###################################################
### code chunk number 12: mnlogit.Rnw:348-349
###################################################
print(fit$est.stats)


###################################################
### code chunk number 13: mnlogit.Rnw:357-358
###################################################
print(fit$model.size)


###################################################
### code chunk number 14: mnlogit.Rnw:856-857
###################################################
library(mlogit)


###################################################
### code chunk number 15: mnlogit.Rnw:862-864
###################################################
source("simChoiceModel.R")
data <- makeModel('X', K=5)


###################################################
### code chunk number 16: mnlogit.Rnw:867-872
###################################################
K = length(unique(data$choices))
N = nrow(data)/K
p = ncol(data) - 3   
np = (K - 1) * p
cat(paste0("Number of choices in simulated data = K = ", K, ".\nNumber of observations in simulated data = N = ", N, ".\nNumber of variables = p = ", p, ".\nNumber of model parameters = (K - 1) * p = ", (K-1)*p, "."))


###################################################
### code chunk number 17: mnlogit.Rnw:876-878
###################################################
vars <- paste("X", 1:50, sep="", collapse=" + ")
fm <- formula(paste("response ~ 1|", vars, " - 1 | 1"))


###################################################
### code chunk number 18: mnlogit.Rnw:881-882
###################################################
system.time(fit.mnlogit <- mnlogit(fm, data, "choices"))  # runs on 1 proc


###################################################
### code chunk number 19: mnlogit.Rnw:885-889
###################################################
mdat <- mlogit.data(data[order(data$indivID), ], "response", shape="long", 
alt.var="choices")
system.time(fit.mlogit <- mlogit(fm, mdat))   # Newton-Raphson
system.time(fit.mlogit <- mlogit(fm, mdat, method='bfgs')) 


###################################################
### code chunk number 20: mnlogit.Rnw:899-903
###################################################
library(nnet)
ndat <- data[which(data$response > 0), ]
ff <- paste("choices ~", vars, "- 1")   # formula for nnet
system.time(fit.nnet <- multinom(ff, ndat, reltol=1e-12)) 


###################################################
### code chunk number 21: mnlogit.Rnw:907-910
###################################################
library(VGAM)
stop.vglm <- vglm.control(epsilon = 1e-6)
system.time(fit.vglm <- vglm(ff, data=ndat, multinomial, control=stop.vglm))


