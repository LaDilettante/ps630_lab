########## Normally Distributed DVs 2 ############

library(arm)
library(geoR)

# Data from 2012 American National Election Study (~3100 white-identifying respondents) #

anesdata <- na.omit(read.delim("W2/2012 ANES_Economic Prefs.txt", header=TRUE))
summary(anesdata)


## Simulation easy example ##

# Fit bivariate OLS

fit1 <- lm(econ01 ~ male, data=anesdata)
summary(fit1)

# Extract estimates for b0 and b1 and their var-cov matrix

mu_hat <- coef(fit1)
mu_se <- vcov(fit1)

# Draw 1000, length-2 vectors from a bivariate normal with mean and var-cov as above

sims1 <- mvrnorm(1000, mu_hat, mu_se)

# Summarize simulations

hist(sims1[ ,1])
hist(sims1[ ,2])

quantile(sims1[ ,1], probs=c(.025,.50,.975))
quantile(sims1[ ,2], probs=c(.025,.50,.975))

confint(fit1)

# Estimate kernal density of that distribution

kden <- kde2d(sims1[,1], sims1[,2])

# Plot the bivariate normal

persp(kden, phi = 40, theta = 60, xlab="b0", ylab="b1", zlab="density")

# Estimate correlation of b0 and b1

vcov(fit1)[2,1]/sqrt(vcov(fit1)[1,1]*vcov(fit1)[2,2])


## OLS regression of support for social welfare on several predictors ##

fit2 <- lm(econ01 ~ age01 + male + income01 + union + unemp + openness01, data=anesdata)
display(fit2)
confint(fit2)

# Coefficient figure #

par(mar=c(5,4,2,2))
plot(coef(fit2)[2:7], 1:6, pch=19, ylim=c(.5,6.5), xlim=c(-.2,.2),
     ylab="", xlab="Marginal effect (min to max)", axes=F)
axis(1, at=c(-.2,-.1,0,.1,.2))
axis(2, at=c(1:6), labels=c("Age","Male","Income","Union","Unemp","Open"), las=2, pos=-.2)
segments(confint(fit2)[2:7,1], 1:6, confint(fit2)[2:7,2], 1:6, lty=1, lwd=1)
abline(v=0, lty=2)
title(main="Predictors of support for social welfare")

coefplot(fit2)

# Adding education via nominal operationalization #

fit3 <- lm(econ01 ~ age01 + male + income01 + union + unemp + openness01
           + educ2 + educ3 + educ4 + educ5, data=anesdata)
display(fit3)

# Take 1000 draws of beta vector #

sims1 <- sim(fit3, n.sims=1000) #returns a matrix w/ 1000 rows and K columns
coef(sims1)[1:10,]

lessBA <- (3*coef(sims1)[ ,1] + coef(sims1)[ ,8] + coef(sims1)[ ,9])/3 #average support for SW for people with <BA degree (holding other vars at 0)
BA <- coef(sims1)[ ,1] + coef(sims1)[ ,10] #average support for SW for people w/ BA degree (holding other vars at 0)

hist(lessBA, xlim=c(.25,.5))
hist(BA, add=T, col="grey")

BAtest <- lessBA - BA
hist(BAtest)

mean(BAtest)
quantile(BAtest, probs=c(.025,.50,.975))


### What is going on 'under the hood' of 'sim'? ###

fit2 <- lm(econ01 ~ age01 + male + income01 + union + unemp + openness01, data=anesdata)
summary.lm(fit2)

#Unscaled covariance matrix
cov <- summary.lm(fit2)$cov.unscaled

#Estimate for sigma
sighat <- summary.lm(fit2)$sigma

## Basic 'sim' procedure ##

sigdraws <- array(NA, c(1000))
betadraws <- array(NA, c(1000, length(coef(fit2))))

for (i in 1:1000){
  sigdraws[i] <- sighat*sqrt((summary.lm(fit2)$df[2])/rchisq(1, summary.lm(fit2)$df[2])) #Take draws of sigma
  betadraws[i, ] <- mvrnorm(1, coef(fit2), cov*(sigdraws[i])^2) #Draw beta vector
}

# The formula for sigdraws given in Gelman is equivalent to taking draws
 # from a scaled-inverse-chisquare distribution with n-k df and scale = sighat^2
 # which is the marginal posterior distribution of sigma^2 given a uniform (uninformative) prior

test1 <- sighat*sqrt(10/rchisq(100000,10))
test2 <- sqrt(rinvchisq(100000,10,sighat^2)) #requires package geoR

hist(test1)
hist(test2, add=T, col="grey")

# Summarize simulations

table_fit2 <- array(NA, c(length(coef(fit2)), 3))
for (i in 1:length(coef(fit2))){
  table_fit2[i, ] <- quantile(betadraws[ ,i], probs=c(.025,.5,.975))
}


## Model w/ interaction between household income and gender ##

fit4 <- lm(econ01 ~ age01 + male + income01 + union + unemp + openness01 + income01:age01, data=anesdata)
display(fit4)
confint(fit4)

# Conditional confidence intervals using Brambor, Clark, and Golder method

age <- c(0,.25,.5,.75,1)

marg <- coef(fit4)[4] + coef(fit4)[8]*age
se_marg <- sqrt(vcov(fit4)[4,4] + vcov(fit4)[8,8]*age^2 + 2*age*vcov(fit4)[8,4])

low <- marg - 1.96*se_marg
high <- marg + 1.96*se_marg

par(mar=c(5,4,2,2))
plot(1:5, marg,
     pch=19, ylim=c(-.3,.1), xlim=c(1,5),
     ylab="Conditional marginal effect of income", xlab="Age (min to max)", axes=F)
lines(1:5, marg, lty=1, lwd=1)
lines(1:5, low, lty=2, lwd=1)
lines(1:5, high, lty=2, lwd=1)
axis(1, at=c(1:5), labels=c(0,.25,.5,.75,1))
axis(2, at=c(-.3,-.2,-.1,0,.1))
abline(h=0, lty=3)
title(main="Interaction of Gender and Income")


# Confidence intervals for conditional marginal effects using simulation

sims2 <- sim(fit4, n.sims=1000)
coef(sims2)[1,]

inc_marg <- array(NA, c(1000,5)) #A 1000 x 5 matrix to hold the 1000 simulations for each of 5 values of age
inc_table <- array(NA, c(5,3)) #A 5 x 3 matrix to hold means of simulations (1st col) and their confidence intervals (2nd and 3rd cols)
age <- c(0,.25,.5,.75,1) #A vector that contains a range of reasonable values for the moderating variable, age

for (j in 1:5){ #Loop over values of age
  for (i in 1:1000){ #Loop over simulations
    inc_marg[i,j] <- coef(sims2)[i,4] + coef(sims2)[i,8]*age[j] #Partial derivative of regression model w.r.t to income
  }
  inc_table[j,1] <- mean(inc_marg[ ,j]) #Average the simulations for a given value of age
  inc_table[j,2:3] <- quantile(inc_marg[ ,j], probs=c(.025,.975)) #Calculate 95% CIs for a given value of age
}

inc_table

par(mar=c(5,4,2,2))
plot(1:5, inc_table[1:5,1],
     pch=19, ylim=c(-.3,.1), xlim=c(1,5),
     ylab="Conditional marginal effect of income", xlab="Age (min to max)", axes=F)
lines(1:5, inc_table[1:5,1], lty=1, lwd=1)
lines(1:5, inc_table[1:5,2], lty=2, lwd=1)
lines(1:5, inc_table[1:5,3], lty=2, lwd=1)
axis(1, at=c(1:5), labels=c(0,.25,.5,.75,1))
axis(2, at=c(-.3,-.2,-.1,0,.1))
abline(h=0, lty=3)
title(main="Interaction of Age and Income")


































