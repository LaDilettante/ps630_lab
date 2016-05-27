# Tutorial 9: In-class exercise
# Solutions

# Create an IV x with normal distribution

x=rnorm(1000,mean=10,sd=4)

# Create a DV y as linear function of x + random error

y=2+2*x+rnorm(1000,0,4)

# Plot x and y and show the regression line

plot(x,y)

reg1=lm(y ~ x)
summary(reg1)

abline(reg1)

# Censor from below by setting all values y that are 1 SD below mean of Y to this value 

y2=y
y2[y < (mean(y)-sd(y))]=mean(y)-sd(y)

summary(y)
summary(y2)

# Plot x and the new y and show the regression line

plot(x,y2, ylim=c(min(y),max(y)))

reg2=lm(y2 ~ x)
summary(reg2)

abline(reg2)

# Exclude all values of y that are 1 SD below the mean of Y

y3=y[y>(mean(y)-sd(y))]
x3=x[y>(mean(y)-sd(y))]

# Plot x and the new y and show the regression line

plot(x3,y3, ylim=c(min(y),max(y)))

reg3=lm(y3 ~ x3)
summary(reg3)

abline(reg3)

# Show all three regression lines in one plot and mark the censored points red

plot(x,y)

abline(reg1)

points(x,y2, col="red")
abline(reg2, col="red")

abline(reg3, col="blue")

# Estimate a tobit model for y2 and x, show it on a plot

library(VGAM)

reg4=vglm(y2 ~ x, tobit(Lower=(mean(y)-sd(y))))

summary(reg4)

plot(x,y2, ylim=c(min(y),max(y)))

abline(a=coef(reg4)[1], b=coef(reg4)[3], col="green")

abline(reg1)
