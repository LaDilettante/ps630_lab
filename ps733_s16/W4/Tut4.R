### TUTORIAL 4

data <- read.table(
  "circuit court abortion cases.txt", header=TRUE)
summary(data)

### First derivative as marginal effect (see lecture) ###

curve(pnorm(-.50 + .50*x), from=-5, to=5, ylab="", xlab="X") # Cumulative normal as function of x
curve(dnorm(-.50 + .50*x)*.50, from=-5, to=5, add=T) # first derivative as function of x
curve((dnorm(-.50 + .50*-2)*.50)*x + pnorm(-.50 + .50*-2) - (dnorm(-.50 + .50*-2)*.50)*-2, from=-5, to=5, add=T, lty=3) # tangent line to CDF at -2
lines(-2, pnorm(-.50 + .50*-2), type="p", pch=19) # point for tangent line at -2
curve((dnorm(-.50 + .50*0)*.50)*x + pnorm(-.50 + .50*0), from=-5, to=5, add=T, lty=3) # tangent line to CDF at 0
lines(0, pnorm(-.50), type="p", pch=19) # point for tangent line at 0


### Zelig stuff ###

library(Zelig)

### Standard method, but I can't figure out how to extract estimates from sim objects ###

z1 <- zelig(libvote ~ ideology + postcasey + caseyprov + scmedian, data=data, model="probit")
summary(z1)

x1 <- setx(z1, ideology=0, postcasey=1, scmedian=quantile(data$scmedian, .05)) # Everything that is not specified is set to central tendencies
x1

s1 <- sim(z1, x1)
summary(s1)

### New method ###

z1 <- zprobit$new() # Define a new zelig probit object

z1$zelig(libvote ~ ideology + postcasey + caseyprov + scmedian, data=data) # Run model
summary(z1)

z1$setx(ideology=0, postcasey=1, scmedian=quantile(data$scmedian, .05)) # Set variables to values
z1$setx1(ideology=0, postcasey=1, scmedian=quantile(data$scmedian, .95)) # Set variables to values
summary(z1)

z1$sim() # Generate simulations
summary(z1)

sims.x <- z1$getqi(qi="ev", xvalue="x") #Store predicted probs in new var
mean(sims.x)
quantile(sims.x, probs=c(.025,.975))

sims.x1 <- z1$getqi(qi="ev", xvalue="x1") #Store predicted probs in new var
mean(sims.x1)
quantile(sims.x1, probs=c(.025,.975))

sims.fd <- z1$getqi(qi="fd", xvalue="x1") #Store marginal effects in new var
mean(sims.fd)
quantile(sims.fd, probs=c(.025,.975))



### Estimation using 'glm' ###

library(arm)

fit4 <- glm(libvote ~ ideology + postcasey + caseyprov + scmedian, data=data, family=binomial(link="logit"))
display(fit4)


fit5 <- glm(libvote ~ ideology + postcasey + caseyprov + scmedian, data=data, family=binomial(link="probit"))
display(fit5)


### Post-estimation via averaging over data ###

# Begin with simulation of beta vector from posterior #

sim5 <- sim(fit5, n.sims=1000)
coef(sim5)[1:10,]


## Calculate predicted probabilities and confidence intervals across SC median ##

scvals <- c(-.1,-.05,0,.05,.1,.15,.2) # 7 values of SC median, from about 5th to about 95th percentile (1 probability for each value)

l1 <- array(NA, c(1000, length(data$libvote), length(scvals))) # 1 prob for each simulation for each obs for each value of SC Median
l2 <- array(NA, c(1000, length(scvals))) # average over data
l3 <- array(NA, c(7,5)) # table to hold estimates and CIs, 1st column is estimate, 2-3 is 95% CI, and 4-5 is 68% CI

l1

for (j in 1:length(scvals)){ # loop over values of scmed
  for (i in 1:1000){ # loop over simulated beta vectors
    
    l1[i, ,j] <- pnorm(coef(sim5)[i,1] + coef(sim5)[i,2]*data$ideology + coef(sim5)[i,3]*data$postcasey
                       + coef(sim5)[i,4]*data$caseyprov + coef(sim5)[i,5]*scvals[j]) # formula for calculating predicted probability
    
    l2[i,j] <- mean(l1[i, ,j]) # calculate average predicted probability across observations
    
  }
  
  l3[j,1] <- mean(l2[ ,j]) # calculate average of average probabilities across simulated beta vectors
  
  l3[j,2:5] <- quantile(l2[ ,j], probs=c(.025,.975,.16,.84)) # calculate quantiles for average probabilities
}

l3

plot(1:7, l3[ ,1], type="l", lty=1, lwd=2, ylim=c(0,1), xlab="SC median (lib -> con)", ylab="p(liberal vote)", axes=F)
lines(1:7, l3[ ,2], type="l", lty=2, lwd=1)
lines(1:7, l3[ ,3], type="l", lty=2, lwd=1)
axis(1, at=c(1,3,5,7), labels=c(-.1,0,.1,.2))
axis(2, at=c(0,.25,.5,.75,1))
title(main="SC Median and Courts of Appeals Abortion Votes - Probit")
legend("topright", c("Predicted prob", "95% CI"), lty=c(1,2), lwd=c(2,1), bty="n", inset=.05)


### Marginal Effects ###

l1 <- array(NA, c(1000,length(data$libvote))) #Generate 1 marginal effect for each observation within each of the 1000 simulations
l2 <- array(NA, c(1000)) #Average over the data to get avg. marginal effect for each simulation

for (i in 1:1000){
  l1[i, ] <- pnorm(coef(sim5)[i,1] + coef(sim5)[i,2]*data$ideology + coef(sim5)[i,3]*data$postcasey 
                   + coef(sim5)[i,4]*data$caseyprov + coef(sim5)[i,5]*quantile(data$scmedian, probs=.95)[1]) - #Predicted value at scmedian=95th percentile
    
    pnorm(coef(sim5)[i,1] + coef(sim5)[i,2]*data$ideology + coef(sim5)[i,3]*data$postcasey 
          + coef(sim5)[i,4]*data$caseyprov + coef(sim5)[i,5]*quantile(data$scmedian, probs=.05)[1]) #Predicted value at scmedian=5th percentile
  
  l2[i] <- mean(l1[i, ])
}


hist(l2)

table5 <- c(mean(l2[]), quantile(l2[], probs=c(.025,.975,.16,.84))) #Avg marginal + 95% and 68% CIs
table5


# Estimation using ppavg

library(foreign)
LDC=read.dta("../LDC_IO_replication.dta")

complete = with(LDC, complete.cases(ecris2, l1polity, gdp_pc_95d, gatt_wto_new, fdignp, lnpop, avnewtar, l1ecris2, yrsoffic, usheg))
LDC = LDC[complete,]

reg1=glm(ecris2 ~ polityiv_update2 + gdp_pc_95d + gatt_wto_new + fdignp  + lnpop + avnewtar + l1ecris2 + yrsoffic + usheg, data=LDC, family=binomial(link = "probit"))
summary(reg1)

library(arm)
sims1=sim(reg1)

summary(LDC$gdp_pc_95d)

ppavg(reg1,sims1,x1name="polityiv_update2",x1vals=c(-10,0,10))

ppavg(reg1,sims1,x1name="polityiv_update2",x1vals=c(-10,0,10),x2name="gdp_pc_95d",x2vals=c(364.10,1954.00))