rm(list = ls())

data <- na.omit(read.delim("W6/2004 ANES_Abortion.txt", header=TRUE))
summary(data)

library(arm)

### Estimate ordered probit and logit ###

data$abort <- as.factor(data$abort) # polr requires that the DV be declared as a factor variable

fit1 <- polr(abort ~ age + male + black + hisp + income + rural + church + wrole, data=data, method="probit")
display(fit1)
summary(fit1)

# There is no p-value by default, when there is a lot of data we can assume that
# the t distribution has a lot of df and approximate the normal distribution
# then we can calculate p value by comparing the t-value against a normal dist

coefficients <- coef(summary(fit1))
p <- pnorm(abs(coefficients[ , "t value"]), lower.tail = F) * 2
coefficients <- cbind(coefficients, p_value = p)
coefficients

# Confidence interval
ci <- confint(fit1)
ci

exp(cbind(coef(fit1), ci))
