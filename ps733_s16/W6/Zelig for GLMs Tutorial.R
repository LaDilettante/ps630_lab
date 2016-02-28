################# Zelig - Basic Tutorial for Use in R w/ GLMs ################

library(Zelig)


data <- read.table(
  "C:/Users/cdj19/Documents/Teaching/Graduate/Classes/MLE/Spring 2016/Data/circuit court abortion cases.txt", header=TRUE)
summary(data)


### Example for generating probs and marginal effects using zelig and simulation from the joint asymptotic sampling distribution ###

## Fit a basic probit regression ##

z1 <- zelig(libvote ~ ideology + postcasey + caseyprov + scmedian, data=data, model="probit")
summary(z1)


## Set predictors to chosen values (by default, zelig sets non-specified variables to central tendencies) ##

x1.1 <- setx(z1, ideology=0, postcasey=1, scmedian=quantile(data$scmedian, .05)) # Everything that is not specified is set to central tendencies
x1.1


## Draw 1000 simulations from multivariate normal with mean = estimates from z1, and var-cov matrix from z1 ##

s1.1 <- sim(z1, x1.1)
summary(s1.1)


## Set predictors to a new value ##

x1.2 <- setx(z1, ideology=0, postcasey=1, scmedian=quantile(data$scmedian, .95)) # Everything that is not specified is set to central tendencies
x1.2


## Calculate CHANGE in probability Y=1 from x1.1 to x1.2 ##

s1.2 <- sim(z1, x1.1, x1.2)
summary(s1.2)


## You can also do this via bootstrapping (draw 1000 datasets w/ replacement from original data) ##

s1.2 <- sim(z1, x1.1, x1.2, bootstrap=TRUE)
summary(s1.2)


#### Working w/ quantities of interest ####

## If you want to work with the simulated quantities of interest directly,
 # you (apparently) need to use Zelig 5 framework
 # This is not very intuitive, but here is the example above using that framework
 # At the end, I have extracted the simulated probs and marginal effects
 # These can be used to make graphs, etc.

# Initialize a new zelig object for a probit regresssion #

zobject <- zprobit$new()


# Estimate the probit model using the zelig function

zobject$zelig(libvote ~ ideology + postcasey + caseyprov + scmedian, data=data)
summary(zobject)

# Set x values #

zobject$setx(ideology=0, postcasey=1, scmedian=quantile(data$scmedian, .05))
zobject$setx1(ideology=0, postcasey=1, scmedian=quantile(data$scmedian, .95))

# Run simulations #

zobject$sim()
summary(zobject)


## Extract particular quantities of interest ##

# Predicted prob for first setx command #

sims.x <- zobject$getqi(qi="ev", xvalue="x")
mean(sims.x)
quantile(sims.x, probs=c(.025,.975))

# Predicted prob for second setx command

sims.x1 <- zobject$getqi(qi="ev", xvalue="x1")
mean(sims.x1)
quantile(sims.x1, probs=c(.025,.975))

# First difference, ie, change in prob from x to x1 #

sims.fd <- zobject$getqi(qi="fd", xvalue="x1")
mean(sims.fd)
quantile(sims.fd, probs=c(.025,.975))




























