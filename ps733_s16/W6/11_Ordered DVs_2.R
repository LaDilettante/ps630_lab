############ Ordered DVs - 2 #############

### Parallel lines assumption ###

curve(pnorm(2 + x), from=-5, to=5, ylab="p(Y>m)")
curve(pnorm(.5 + x), from=-5, to=5, add=T)
curve(pnorm(-.75 + x), from=-5, to=5, add=T)
curve(pnorm(-1.5 + x), from=-5, to=5, add=T)
abline(v=c(-2,-.5,.75,1.5), lty=3)

### Examples for probs ###

data <- na.omit(read.delim(
  "W6/2004 ANES_Abortion.txt", header=TRUE))
summary(data)

library(arm)

### Estimate ordered probit and logit ###

data$abort <- as.factor(data$abort) # polr requires that the DV be declared as a factor variable

fit1 <- polr(abort ~ age + male + black + hisp + income + rural + church + wrole, data=data, method="probit")
display(fit1)

fit2 <- polr(abort ~ age + male + black + hisp + income + rural + church + wrole, data=data, method="logistic") # note "logistic" not "logit" as with glm
display(fit2)


### Series of BRMs specification (now includes (m-1) "intercepts" with cut point = 0 ###

library(VGAM)

fit3 <- vglm(abort ~ age + male + black + hisp + income + rural + church + wrole, data=data,
            cumulative(parallel = T, reverse=T, link = "logit"))
summary(fit3)
summary(fit2)


## using predict ##

newdata <- data.frame(age=c(rep(47,7)), male=c(rep(0,7)), black=c(rep(0,7)),
                      hisp=c(rep(0,7)), income=c(rep(16,7)), rural=c(rep(0,7)),
                      church=c(rep(2,7)), wrole=c(1:7))

pp_wrole <- predict(fit4, newdata=newdata, type="probs")
pp_wrole

plot(1:7, pp_wrole[1:7,1], type="l", lty=1, ylab="p(Y=m)",
     xlab="Belief about women's role in society (lib -> con)",
     ylim=c(0,1), xlim=c(.75,7.25), axes=F)
lines(1:7, pp_wrole[1:7,2], lty=2)
lines(1:7, pp_wrole[1:7,3], lty=3)
lines(1:7, pp_wrole[1:7,4], lty=4)
axis(1, at=c(1:7))
axis(2, at=c(0,.25,.5,.75,1))
legend("topright", c("Never legal","Rape, incest,etc.","Some restrictions","Always legal"),
       lty=c(1,2,3,4), bty="n", inset=.05)
title(main="Beliefs about women's role and abortion preferences")


## using Zelig ##

library(Zelig)
library(ZeligChoice)


## using simulation and observed value (averaging over data) approach ##

detach("package:Zelig", unload = T)

sim2 <- sim(fit2, n.sims=1000)
coef(sim2)[1,]



### For all probs across women's role, an example ###

X_temp <- model.matrix(polr(abort ~ age + male + black + hisp + income + rural + church + wrole, data=data, method="logistic"))
X_temp <- X_temp[ ,2:ncol(X_temp)] # get rid of intercept

X <- array(NA, c(length(data$abort), ncol(X_temp), 7)) # 7 IV arrays, 1 each for the 7 values of wrole

for (i in 1:7){
  X[ , ,i] <- X_temp
  X[ , 8, i] <- i # set wrole column to each of 7 values
}


l1 <- array(NA, c(1000, length(X[,1,1]), 7, 4))
l2 <- array(NA, c(1000, 7, 4))
l3 <- array(NA, c(7, 3, 4))

for (j in 1:7){
  for (i in 1:1000){
    l1[i, ,j,1] <- plogis(coef(sim2)[i,1] - X[,,j] %*% coef(sim2)[i,4:ncol(coef(sim2))])
    l1[i, ,j,2] <- plogis(coef(sim2)[i,2] - X[,,j] %*% coef(sim2)[i,4:ncol(coef(sim2))]) - plogis(coef(sim2)[i,1] - X[,,j] %*% coef(sim2)[i,4:ncol(coef(sim2))])
    l1[i, ,j,3] <- plogis(coef(sim2)[i,3] - X[,,j] %*% coef(sim2)[i,4:ncol(coef(sim2))]) - plogis(coef(sim2)[i,2] - X[,,j] %*% coef(sim2)[i,4:ncol(coef(sim2))])
    l1[i, ,j,4] <- plogis(-coef(sim2)[i,3] + X[,,j] %*% coef(sim2)[i,4:ncol(coef(sim2))])

    l2[i,j,1] <- mean(l1[i,,j,1])
    l2[i,j,2] <- mean(l1[i,,j,2])
    l2[i,j,3] <- mean(l1[i,,j,3])
    l2[i,j,4] <- mean(l1[i,,j,4])
  }
  l3[j,1,1] <- mean(l2[,j,1])
  l3[j,1,2] <- mean(l2[,j,2])
  l3[j,1,3] <- mean(l2[,j,3])
  l3[j,1,4] <- mean(l2[,j,4])

  l3[j,2:3,1] <- quantile(l2[,j,1], probs=c(.025,.975))
  l3[j,2:3,2] <- quantile(l2[,j,2], probs=c(.025,.975))
  l3[j,2:3,3] <- quantile(l2[,j,3], probs=c(.025,.975))
  l3[j,2:3,4] <- quantile(l2[,j,4], probs=c(.025,.975))
}

l3[,,1] # mean probabilities and 95% CIs for Y=1
l3[,,2] # mean probabilities and 95% CIs for Y=2
l3[,,3] # mean probabilities and 95% CIs for Y=3
l3[,,4] # mean probabilities and 95% CIs for Y=4


## using post ##

library(arm)




### Fit stats ###

# PCP #

Y.hat <- predict(fit2)
Y <- data$abort

cp <- 0
for (i in 1:length(Y)){
  if (Y.hat[i] == Y[i]){cp <- cp + 1}
}
pcp <- cp/length(Y)

# PRE #

pre <- (pcp - max(table(Y)/length(Y)))/(1-max(table(Y)/length(Y)))


# ePCP #

p2 <- predict(fit2, type="probs")
Y <- data$abort

s <- c(rep(NA,length(levels(Y))))

for (j in 1:length(levels(Y))){
  s[j] <- 0
  for (i in 1:length(Y)){
    if (Y[i] == j){s[j] <- s[j] + p2[i,j]} else{s[j] <- s[j]}
  }
}

epcp <- sum(s)/length(Y)
































