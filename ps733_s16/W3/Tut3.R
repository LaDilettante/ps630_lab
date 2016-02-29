library(arm)

# DATA

data <- read.table(
  "circuit court abortion cases.txt", header=TRUE)
summary(data)

# LOGIT

fit4 <- glm(libvote ~ ideology + postcasey + caseyprov + scmedian, data=data, family=binomial(link="logit"))

display(fit4)

# PROBIT

fit5 <- glm(libvote ~ ideology + postcasey + caseyprov + scmedian, data=data, family=binomial(link="probit"))

display(fit5)

### Predicted probabilities using the 'predict' function ###

# Hold ideology at median, and postcasey and caseyprov at 0, and vary SC median across its range #

newdata1 <- data.frame(ideology=c(rep(median(data$ideology),7)), 
                       postcasey=c(rep(0,7)), 
                       caseyprov=c(rep(0,7)), 
                       scmedian=c(-.1,-.05,0,.05,.1,.15,.2))
newdata1

p4 <- predict(fit5, newdata=newdata1, se.fit=TRUE, type="response")

p4
p4$fit
p4$se


# Hold ideology at median, and postcasey and caseyprov at 1, and vary SC median across its range #

newdata2 <- data.frame(ideology=c(rep(median(data$ideology),7)), 
                       postcasey=c(rep(1,7)), 
                       caseyprov=c(rep(1,7)), 
                       scmedian=c(-.1,-.05,0,.05,.1,.15,.2))
newdata2

p5 <- predict(fit5, newdata=newdata2, se.fit=TRUE, type="response")

p5
p5$fit
p5$se

# Differences in probability change

p4$fit[1]-p4$fit[2]
p5$fit[1]-p5$fit[2]


## Marginal effect for change in SC median from -.10 to .20 on p(Y=1) ##

# With other predictors set to 0 #

marg1 <- dnorm(coef(fit5)[1] + coef(fit5)[2]*0 + coef(fit5)[3]*0 + coef(fit5)[4]*0 + coef(fit5)[5]*-.10)*(coef(fit5)[5]*.30)
marg1

# With other predictors set to different values #

marg2 <- dnorm(coef(fit5)[1] + coef(fit5)[2]*-.5 + coef(fit5)[3]*1 + coef(fit5)[4]*1 + coef(fit5)[5]*-.10)*(coef(fit5)[5]*.30)
marg2

# With intercept set to -2 #

marg3 <- dnorm(-2 + coef(fit5)[2]*-.5 + coef(fit5)[3]*1 + coef(fit5)[4]*1 + coef(fit5)[5]*-.10)*(coef(fit5)[5]*.30)
marg3

# With intercept set to 2 #

marg4 <- dnorm(2 + coef(fit5)[2]*-.5 + coef(fit5)[3]*1 + coef(fit5)[4]*1 + coef(fit5)[5]*-.10)*(coef(fit5)[5]*.30)
marg4

## Odds for logit ##

# Example #
coef(fit4)[[1]]

odds1 <- exp(coef(fit4)[[1]] + coef(fit4)[[2]]*0 + coef(fit4)[[3]]*0 + coef(fit4)[[4]]*0 + coef(fit4)[[5]]*0)
odds1

odds2 <- exp(coef(fit4)[[1]] + coef(fit4)[[2]]*0 + coef(fit4)[[3]]*0 + coef(fit4)[[4]]*1 + coef(fit4)[[5]]*0)
odds2

odds2/odds1

exp(coef(fit4)[[4]]*1)


####### Zelig - Basic Tutorial for Use in R w/ GLMs ########

library(Zelig)

### Example for generating marginal effects using zelig and simulation from the joint asymptotic sampling distribution ###


## Fit a basic probit regression ##

z1 <- zelig(libvote ~ ideology + postcasey + caseyprov + scmedian, data=data, model="probit")
summary(z1)


## Set predictors to chosen values (by default, zelig sets non-specified variables to means) ##

x1.1 <- setx(z1, ideology=0, postcasey=1, scmedian=quantile(data$scmedian, .05)) # Everything that is not specified is set to mean
x1.1


## Draw 1000 simulations from multivariate normal with mean = estimates from z1, and var-cov matrix from z1 ##

s1.1 <- sim(z1, x1.1)
summary(s1.1)


## Set predictors to a new value ##

x1.2 <- setx(z1, ideology=0, postcasey=1, scmedian=quantile(data$scmedian, .95)) # Everything that is not specified is set to central tendencies
x1.2


## Calculate CHANGE in probability Y=1 from x1.1 to x1.2 (marginal effect of scmedian) ##

s1.2 <- sim(z1, x1.1, x1.2) #Calculate change from x1.1 to x1.2 (x1.2 - x1.1)
summary(s1.2)

summary(s1.2)$fd #fd stands for 'first difference'


## You can also do this via bootstrapping (draw 1000 datasets w/ replacement from original data) ##

s1.2 <- sim(z1, x1.1, x1.2, bootstrap=TRUE)
summary(s1.2)

summary(s1.2)$fd 