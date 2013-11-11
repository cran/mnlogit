### R code from vignette source 'mnlogit.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: mnlogit.Rnw:101-102
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: mnlogit.Rnw:152-153
###################################################
library(mnlogit)


###################################################
### code chunk number 3: mnlogit.Rnw:159-161
###################################################
data(Fish, package = 'mnlogit')
head(Fish, 8)


###################################################
### code chunk number 4: mnlogit.Rnw:254-255 (eval = FALSE)
###################################################
## fm <- formula(mode ~ price | income | catch)


###################################################
### code chunk number 5: mnlogit.Rnw:259-262 (eval = FALSE)
###################################################
## fm <- formula(mode ~ price | income - 1 | catch)
## fm <- formula(mode ~ price | income | catch - 1)
## fm <- formula(mode ~ 0 + price | income | catch)


###################################################
### code chunk number 6: mnlogit.Rnw:266-271 (eval = FALSE)
###################################################
## fm <- formula(mode ~ 1 | income | catch) 
## fm <- formula(mode ~ price | 1 | catch) 
## fm <- formula(mode ~ price | income | 1)
## fm <- formula(mode ~ price | 1 | 1)
## fm <- formula(mode ~ 1 | 1 | price + catch)


###################################################
### code chunk number 7: mnlogit.Rnw:275-278 (eval = FALSE)
###################################################
## fm <- formula(mode ~ price + catch | 1 | 1)
## fm <- formula(mode ~ price + catch | 1) 
## fm <- formula(mode ~ price + catch)


###################################################
### code chunk number 8: mnlogit.Rnw:288-289 (eval = FALSE)
###################################################
## ?mnlogit


###################################################
### code chunk number 9: mnlogit.Rnw:292-295 (eval = FALSE)
###################################################
## mnlogit(formula, data, choiceVar, maxiter = 25, ftol = 1e-6,
##         gtol = 1e-6, ncores = 1, na.rm = TRUE, print.level = 0, 
##          linDepTol = 1e-6, ...)


###################################################
### code chunk number 10: mnlogit.Rnw:304-307
###################################################
fm <- formula(mode ~ price | income | catch)
fit <- mnlogit(fm, Fish, "alt", ncores=2)
class(fit)


###################################################
### code chunk number 11: mnlogit.Rnw:311-312
###################################################
print(fit$est.stats)


###################################################
### code chunk number 12: mnlogit.Rnw:320-321
###################################################
print(fit$model.size)


###################################################
### code chunk number 13: mnlogit.Rnw:346-347 (eval = FALSE)
###################################################
## fm <- formula(mode ~ 1 | income | price + catch)


###################################################
### code chunk number 14: mnlogit.Rnw:728-730
###################################################
library(nnet)
library(mlogit)


###################################################
### code chunk number 15: mnlogit.Rnw:735-738
###################################################
source("simChoiceModel.R")
data <- makeModel('X', K=5)
dim(data)


###################################################
### code chunk number 16: mnlogit.Rnw:742-744
###################################################
vars <- paste("X", 1:50, sep="", collapse=" + ")
fm <- formula(paste("response ~ 1|", vars, "| 1"))


###################################################
### code chunk number 17: mnlogit.Rnw:747-748
###################################################
system.time(fit.mnlogit <- mnlogit(fm, data, "choices"))  # runs on 1 proc


###################################################
### code chunk number 18: mnlogit.Rnw:751-755
###################################################
mdat <- mlogit.data(data[order(data$indivID), ], "response", shape="long", 
                    alt.var="choices")
system.time(fit.mlogit <- mlogit(fm, mdat))   # Newton-Raphson
system.time(fit.mlogit <- mlogit(fm, mdat, method='bfgs')) 


###################################################
### code chunk number 19: mnlogit.Rnw:762-765
###################################################
ndat <- data[which(data$response > 0), ]
ff <- paste("choices ~", vars)   # formula for nnet
system.time(fit.nnet <- multinom(ff, ndat, reltol=1e-10, abstol=1e-8)) 


