########## Normally Distributed DVs 1 ############

library(bbmle)
library(arm)

# Data from 2012 American National Election Study (~3100 non-Latino-white-identifying respondents)

anesdata <- na.omit(read.delim(
  "./2012 ANES_Economic Prefs.txt", header=TRUE))
summary(anesdata)

### Simple linear regression of social welfare support on income (coded 0-1) and union membership (dichot) ###

X <- array(NA, c(length(anesdata$econ01),3)) #Initialize the design matrix
X[ ,1] <- 1 #column of 1s for constant
X[ ,2] <- anesdata$income01 #household income, coded from 0-1
X[ ,3] <- anesdata$union #dichotomous indicator for family member union membership
y <- anesdata$econ01 #additive scale of several social welfare policy items, coded from 0-1

# Estimation via OLS #

summary(lm(y ~ X[,2] + X[,3]))

## Estimation via MLE using iterative method for maximizing LL ##

# Define the LL function in general terms (this code can be used for any # of predictors):

LL_normreg = function(params, y, X){
  B = matrix(NA, nrow = length(params) - 1, ncol = 1)
  B[,1] = params[-length(params)]
  sigma    = params[[length(params)]]
  minusll  = -sum(dnorm(y, X %*% B, sigma, log=T))
  return(minusll)
}

# Declare the names of the parameters (from B0 to B[# of predictors], and sigma):

parnames(LL_normreg) <- c("B0", "B1", "B2", "sigma")

# Fit the model using mle2 ('vecpar=TRUE' tells mle2 that the first argument passed to the
  # LL function is a vector of all parameters with names declared in 'parnames' above and in the start values):

fit <- mle2(LL_normreg, start = c(B0 = mean(y), B1 = 0, B2 = 0, sigma = sd(y)),
            data=list(y=y,X=X), vecpar = TRUE, control=list(maxit=5000))

summary(fit)
confint(fit)
vcov(fit)
sqrt(vcov(fit))


### Modeling the Variance ###

# Logistic link function example #

curve(plogis(x), from=-4, to=4)

prop <- runif(100,0,1)
hist(prop)

hist(logit(prop))

# Exponential link function example #

curve(exp(x), from=-2, to=2)


## Heteroskedastic regression
 # Model for variance: sigma^2 = exp(g0 + g1*age)

LL_hetreg <- function(b0, b1, b2, g0, g1, data){
  -sum( log(dnorm( (y - b1*x1 - b2*x2 - b0)/sqrt(exp((g0 + g1*z1))) )) - (g0 + g1*z1)/2 )
}

fit2 <- mle2(LL_hetreg, start=list(b0=mean(anesdata$econ01), b1=0, b2=0, g0=sd(anesdata$econ01), g1=0),
             data=list(y=anesdata$econ01, x1=anesdata$income01, x2=anesdata$union, z1=anesdata$age01))

summary(fit2)
confint(fit2)
vcov(fit2)

# Plot relationship of age to residual standard deviation

par(mar=c(5,4,2,2)) #Set margins of plot
age <- c(0,.25,.5,.75,1) #Create new var called 'age' that takes on 5 values from 0 to 1
plot(age, sqrt(exp(coef(fit2)[4] + coef(fit2)[5]*age)),
     type="p", ylab="Residual standard deviation", xlab="Age (min to max)", axes=F, ylim=c(.2,.25))
lines(age, sqrt(exp(coef(fit2)[4] + coef(fit2)[5]*age)), type="l")
axis(1, at=c(0,.25,.5,.75,1))
axis(2, at=c(.2,.225,.25))
title(main="The effect of age on the residual standard deviation")


## Heteroskedastic regression
 # Model for variance: sigma^2 = exp(g0 + g1*age + g2*conscientiousness)

LL_hetreg2 <- function(b0, b1, b2, g0, g1, g2, data){
  -sum( log(dnorm( (y - b1*x1 - b2*x2 - b0)/sqrt(exp((g0 + g1*z1 + g2*z2))) )) - (g0 + g1*z1 + g2*z2)/2 )
}

fit3 <- mle2(LL_hetreg2, start=list(b0=mean(anesdata$econ01), b1=0, b2=0, g0=sd(anesdata$econ01), g1=0, g2=0),
             data=list(y=anesdata$econ01, x1=anesdata$income01, x2=anesdata$union, z1=anesdata$age01, z2=anesdata$consc01))

summary(fit3)
confint(fit3)
vcov(fit3)











































