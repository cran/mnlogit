### R code from vignette source 'mnlogit.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: mnlogit.Rnw:106-107
###################################################
options(prompt = "R> ", continue = "+  ", width = 60, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: mnlogit.Rnw:158-159
###################################################
library("mnlogit")


###################################################
### code chunk number 3: mnlogit.Rnw:167-169
###################################################
data("Fish", package = 'mnlogit')
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
### code chunk number 8: mnlogit.Rnw:310-313 (eval = FALSE)
###################################################
## mnlogit(formula, data, choiceVar = NULL, maxiter = 50, ftol = 1e-6, gtol
##         = 1e-6, weights = NULL, ncores = 1, na.rm = TRUE, print.level = 0, 
##         linDepTol = 1e-6, start = NULL, alt.subset = NULL, ...)


###################################################
### code chunk number 9: mnlogit.Rnw:322-323
###################################################
fm <- formula(mode ~ price | income | catch)


###################################################
### code chunk number 10: mnlogit.Rnw:340-342 (eval = FALSE)
###################################################
## fit <- mnlogit(fm, Fish, ncores=2)
## class(fit)


###################################################
### code chunk number 11: mnlogit.Rnw:344-346
###################################################
fit <- mnlogit(fm, Fish, choiceVar="alt", ncores=2)
class(fit)


###################################################
### code chunk number 12: mnlogit.Rnw:351-352
###################################################
print(fit, what = "eststat")


###################################################
### code chunk number 13: mnlogit.Rnw:360-361
###################################################
print(fit, what = "modsize")


###################################################
### code chunk number 14: mnlogit.Rnw:374-376 (eval = FALSE)
###################################################
## library("mnlogit")
## ?lrtest


###################################################
### code chunk number 15: mnlogit.Rnw:881-882
###################################################
library("mlogit")


###################################################
### code chunk number 16: mnlogit.Rnw:887-890
###################################################
source("simChoiceModel.R")
numCovars = 50
data <- makeModel('X', K=10, numCovars=numCovars)


###################################################
### code chunk number 17: mnlogit.Rnw:893-898
###################################################
K = length(unique(data$choices))
N = nrow(data)/K
p = ncol(data) - 3   
np = (K - 1) * p
cat(paste0("Number of choices in simulated data = K = ", K, ".\nNumber of observations in simulated data = N = ", N, ".\nNumber of variables = p = ", p, ".\nNumber of model parameters = (K - 1) * p = ", (K-1)*p, "."))


###################################################
### code chunk number 18: mnlogit.Rnw:902-904
###################################################
vars <- paste("X", 1:numCovars, sep="", collapse=" + ")
fm <- formula(paste("response ~ 1|", vars, " - 1 | 1"))


###################################################
### code chunk number 19: mnlogit.Rnw:907-908
###################################################
system.time(fit.mnlogit <- mnlogit(fm, data, "choices"))  


###################################################
### code chunk number 20: mnlogit.Rnw:911-915
###################################################
mdat <- mlogit.data(data[order(data$indivID), ], "response", shape="long", 
alt.var="choices")
system.time(fit.mlogit <- mlogit(fm, mdat))   # Newton-Raphson
system.time(fit.mlogit <- mlogit(fm, mdat, method='bfgs')) 


###################################################
### code chunk number 21: mnlogit.Rnw:925-929
###################################################
library("nnet")
ndat <- data[which(data$response > 0), ]
fm.nnet <- paste("choices ~", vars, "- 1")   
system.time(fit.nnet <- multinom(fm.nnet, ndat, reltol=1e-12)) 


###################################################
### code chunk number 22: mnlogit.Rnw:933-936
###################################################
library("VGAM")
system.time(fit.vglm <- vglm(fm.nnet, data=ndat, multinomial(refLevel=1),
control=vglm.control(epsilon=1e-6)))


