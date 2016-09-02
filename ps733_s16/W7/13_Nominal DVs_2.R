########### Multinomial Logit - 2 ############

data <- na.omit(read.delim(
  "C:/Users/cdj19/Documents/Teaching/Graduate/Classes/MLE/Spring 2016/Data/2004 ANES_Abortion.txt", header=TRUE))
summary(data)

data$abort <- as.factor(data$abort)

### Using VGAM with simulation to get average probs and effects (observed value approach) ###

library(arm)

fit3 <- vglm(abort ~ male + church + income + rural, data=data, multinomial(refLevel=1))
summary(fit3)

# Simulate from sampling distribution of parameters #

sim3 <- mvrnorm(1000, coef(fit3), vcov(fit3))
sim3[1,]

# Calculate probs with church moving from 0 to 4 #

p1.1 <- array(NA, c(1000,length(data$abort),5))
p2.1 <- array(NA, c(1000,length(data$abort),5))
p3.1 <- array(NA, c(1000,length(data$abort),5))
p4.1 <- array(NA, c(1000,length(data$abort),5))

p1.2 <- array(NA, c(1000,5))
p2.2 <- array(NA, c(1000,5))
p3.2 <- array(NA, c(1000,5))
p4.2 <- array(NA, c(1000,5))

p1.3 <- array(NA, c(5,3))
p2.3 <- array(NA, c(5,3))
p3.3 <- array(NA, c(5,3))
p4.3 <- array(NA, c(5,3))

church <- c(0:4)

for (j in 1:length(church)){
  for (i in 1:1000){
    
    # Calculate exp(XB) for each of the 3 models
    eXB <- array(NA, c(3,length(data$abort)))
    for (k in 1:3){
      eXB[k, ] <- exp(sim3[i,k] + sim3[i,k+1*3]*data$male + sim3[i,k+2*3]*church[j] + sim3[i,k+3*3]*data$income + sim3[i,k+4*3]*data$rural)
    }
    
    # Calculate probs of each outcome
    p1.1[i, ,j] <- 1 / (1 + eXB[1,] + eXB[2,] + eXB[3,])
    p2.1[i, ,j] <- eXB[1,] / (1 + eXB[1,] + eXB[2,] + eXB[3,])
    p3.1[i, ,j] <- eXB[2,] / (1 + eXB[1,] + eXB[2,] + eXB[3,])
    p4.1[i, ,j] <- eXB[3,] / (1 + eXB[1,] + eXB[2,] + eXB[3,])
    
    # Average over data
    p1.2[i,j] <- mean(p1.1[i, ,j])
    p2.2[i,j] <- mean(p2.1[i, ,j])
    p3.2[i,j] <- mean(p3.1[i, ,j])
    p4.2[i,j] <- mean(p4.1[i, ,j])
  }
  # Average over sims
  p1.3[j,1] <- mean(p1.2[,j])
  p2.3[j,1] <- mean(p2.2[,j])
  p3.3[j,1] <- mean(p3.2[,j])
  p4.3[j,1] <- mean(p4.2[,j])
  
  # 95% CIs
  p1.3[j,2:3] <- quantile(p1.2[,j], probs=c(.025,.975))
  p2.3[j,2:3] <- quantile(p2.2[,j], probs=c(.025,.975))
  p3.3[j,2:3] <- quantile(p3.2[,j], probs=c(.025,.975))
  p4.3[j,2:3] <- quantile(p4.2[,j], probs=c(.025,.975))
}

p1.3
p2.3
p3.3
p4.3


# Calculate effect of church on each prob #

p.1 <- array(NA, c(1000,length(data$abort),4))
p.2 <- array(NA, c(1000,4))
p.3 <- array(NA, c(4,3))

for (i in 1:1000){
  
  # Calculate exp(XB) for each of the 3 models
  eXB1 <- array(NA, c(3,length(data$abort)))
  eXB2 <- array(NA, c(3,length(data$abort)))
  for (k in 1:3){
    eXB1[k, ] <- exp(sim3[i,k] + sim3[i,k+1*3]*data$male + sim3[i,k+2*3]*4 + sim3[i,k+3*3]*data$income + sim3[i,k+4*3]*data$rural)
    eXB2[k, ] <- exp(sim3[i,k] + sim3[i,k+1*3]*data$male + sim3[i,k+2*3]*0 + sim3[i,k+3*3]*data$income + sim3[i,k+4*3]*data$rural)
  }
  
  # Calculate change in prob of each outcome
  p.1[i, ,1] <- (1 / (1 + eXB1[1,] + eXB1[2,] + eXB1[3,])) - (1 / (1 + eXB2[1,] + eXB2[2,] + eXB2[3,]))
  p.1[i, ,2] <- (eXB1[1,] / (1 + eXB1[1,] + eXB1[2,] + eXB1[3,])) - (eXB2[1,] / (1 + eXB2[1,] + eXB2[2,] + eXB2[3,]))
  p.1[i, ,3] <- (eXB1[2,] / (1 + eXB1[1,] + eXB1[2,] + eXB1[3,])) - (eXB2[2,] / (1 + eXB2[1,] + eXB2[2,] + eXB2[3,]))
  p.1[i, ,4] <- (eXB1[3,] / (1 + eXB1[1,] + eXB1[2,] + eXB1[3,])) - (eXB2[3,] / (1 + eXB2[1,] + eXB2[2,] + eXB2[3,]))
  
  # Average over data
  p.2[i,1] <- mean(p.1[i, ,1])
  p.2[i,2] <- mean(p.1[i, ,2])
  p.2[i,3] <- mean(p.1[i, ,3])
  p.2[i,4] <- mean(p.1[i, ,4])
}
# Average over sims
p.3[1,1] <- mean(p.2[,1])
p.3[2,1] <- mean(p.2[,2])
p.3[3,1] <- mean(p.2[,3])
p.3[4,1] <- mean(p.2[,4])

# 95% CIs
p.3[1,2:3] <- quantile(p.2[,1], probs=c(.025,.975))
p.3[2,2:3] <- quantile(p.2[,2], probs=c(.025,.975))
p.3[3,2:3] <- quantile(p.2[,3], probs=c(.025,.975))
p.3[4,2:3] <- quantile(p.2[,4], probs=c(.025,.975))

p.3


### Fit stats ###

# PCP #

Y.p <- predict(fit3, type="response")
Y.hat <- apply(Y.p, 1, function(x) which.max(x))
Y <- data$abort

cp <- 0
for (i in 1:length(Y)){
  if (Y.hat[i] == Y[i]){cp <- cp + 1}
}
pcp <- cp/length(Y)
pcp


# PRE #

pre <- (pcp - max(table(Y)/length(Y)))/(1-max(table(Y)/length(Y)))
pre


# ePCP #

Y.p <- predict(fit3, type="response")
Y <- data$abort

s <- c(rep(NA,ncol(Y.p)))

for (j in 1:ncol(Y.p)){
  s[j] <- 0
  for (i in 1:length(Y)){
    if (Y[i] == j){s[j] <- s[j] + Y.p[i,j]} else{s[j] <- s[j]}
  }
}

epcp <- sum(s)/length(Y)
epcp


### IIA Assumption using package mlogit ###

library(mlogit)

# Describe data to mnlogit (wide means current data has 1 row for each unit) #

mnldata <- mlogit.data(data, choice="abort", shape="wide")
mnldata[1:10,]

# Specify formula (the 1 is indicating that there are no alternative-specific attributes) #

f1 <- mFormula(abort ~ 1 | male + church + income + rural)

# Fit model #

fit1 <- mlogit(f1, mnldata)
summary(fit1)


### Test IIA assumption w/ Hausman-McFadden test ###

# Estimate restricted model #

fit2 <- mlogit(f1, mnldata, alt.subset = c("1","2","3"))
summary(fit2)

# Test against unrestricted model #

hmftest(fit1, fit2)


### Multinomial Probit ###

## Using mlogit package [this will take awhile, and I couldn't get it to work] ##

fit3 <- mlogit(f1, mnldata, probit=T)
summary(fit3)


## Using MNP package ##

library(MNP)

fit4 <- mnp(abort ~ male + rural, data=data, base=1, n.draws=20000, verbose=T) 
summary(fit4)

newdata <- data.frame(abort=c(0,0), male=c(0,0), rural=c(0,1)) # It wanted the DV in there for some reason, doesn't change anything

p4 <- predict(fit4, newdata = newdata, type = "prob", n.draws = 100, verbose = TRUE)

apply(p4$p, c(1,2), mean) # average over draws for rows and columns of obs X categories matrix


## Full example w/ convergence checks ##

library(coda)

# Run three times with 3 different sets of starting values to check convergence [not converging unless trace is set to False!!!] #
# Getting Truncnorm error a large percentage of the time (but not always) #

fit4a <- mnp(abort ~ male + rural, data=data, n.draws=20000, verbose=T, trace=F) 
fit4b <- mnp(abort ~ male + rural, data=data, n.draws=20000, verbose=T, trace=F, 
             coef.start = c(-1,1,-1,1,-1,1,-1,1,-1), 
             cov.start = matrix(.5, ncol=3, nrow=3) + diag(.5, 3)) 
fit4c <- mnp(abort ~ male + rural, data=data, n.draws=20000, verbose=T, trace=F,
             coef.start = c(1,-1,1,-1,1,-1,1,-1,1),
             cov.start = matrix(.9, ncol=3, nrow=3) + diag(.1, 3))


fit4.coda <- mcmc.list(chain1=mcmc(fit4a$param[,-10]), #make sure to remove the first diagonal entry in the cov matrix (fixed to 1)
                       chain2=mcmc(fit4b$param[,-10]),
                       chain3=mcmc(fit4c$param[,-10]))

# Check convergence #

gelman.diag(fit4.coda, transform=T)

# summarize over chains #

summary(fit4.coda)



### Conditional Logit ###

library(mlogit)

data("TravelMode", package = "AER")
head(TravelMode)

TM <- mlogit.data(TravelMode, choice="choice", shape="long", alt.var="mode", chid.var="individual")
TM[1:20,]


# Only choice-associated variables (pure conditional logit) #

f1 <- mFormula(choice ~ vcost + travel) 

fit1 <- mlogit(f1, TM)
summary(fit1)


air <- exp(0 - .019*(59) - .004*(100)) / (  exp(0 - .019*(59) - .004*(100)) +
                                       exp(.383 - .019*(25) - .004*(417)) + 
                                       exp(.579 - .019*(10) - .004*(180)) +
                                       exp(1.33 - .019*(31) - .004*(372)) )


train <- exp(.383 - .019*(25) - .004*(417)) / (  exp(0 - .019*(59) - .004*(100)) +
                                       exp(.383 - .019*(25) - .004*(417)) + 
                                       exp(.579 - .019*(10) - .004*(180)) +
                                       exp(1.33 - .019*(31) - .004*(372)) )

bus <- exp(.579 - .019*(10) - .004*(180)) / (  exp(0 - .019*(59) - .004*(100)) +
                                          exp(.383 - .019*(25) - .004*(417)) + 
                                          exp(.579 - .019*(10) - .004*(180)) +
                                          exp(1.33 - .019*(31) - .004*(372)) )

car <- exp(1.33 - .019*(31) - .004*(372)) / (  exp(0 - .019*(59) - .004*(100)) +
                                          exp(.383 - .019*(25) - .004*(417)) + 
                                          exp(.579 - .019*(10) - .004*(180)) +
                                          exp(1.33 - .019*(31) - .004*(372)) )

cbind(air,train,bus,car)

predict(fit1, newdata=TM[1:4,])


# Add individual-associated IVs (mixed conditional and multinomial logit) #

f2 <- mFormula(choice ~ vcost + travel | income + size) 

fit2 <- mlogit(f2, TM)
summary(fit2)










































