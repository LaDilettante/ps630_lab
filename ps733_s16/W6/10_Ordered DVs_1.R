############ Ordered DVs - 1 #############

###### Rare events logit (need package Zelig) and Zelig command 'relogit' ######

### Bias in intercept as function of sample size and average p(Y=1) [see King and Zeng for equation] ###

n <- c(500,1000,10000,100000)

curve((x-.5)/(100*x*(1-x)), from=0, to=1, ylab="bias in B0", xlab="average p(Y=1)")
for (i in 1:length(n)){
  curve((x-.5)/(n[i]*x*(1-x)), from=0, to=1, add=T)
}

curve((x-.5)/(100*x*(1-x)), from=0, to=1, xlim=c(0,.05), ylab="bias in B0", xlab="average p(Y=1)")
for (i in 1:length(n)){
  curve((x-.5)/(n[i]*x*(1-x)), from=0, to=1, add=T)
}

### Use relogit (use standard Zelig post-estimation following the zelig command ###

library(Zelig)
data(mid)
summary(mid)

rel.conflict <- zelig(conflict ~ major + contig + power + maxdem + mindem + years, data = mid, model = "relogit", tau = .00343)
summary(rel.conflict)


###### Ordered DVs ######

data <- na.omit(read.delim(
  "C:/Users/cdj19/Documents/Teaching/Graduate/Classes/MLE/Spring 2016/Data/2004 ANES_Abortion.txt", header=TRUE))
summary(data)

detach("package:Zelig", unload=T)
library(arm)


### Latent variable approach ###

curve(dnorm(x), from=-3, to=3, ylab="density", xlab="", axes=F)
axis(1, at=c(-3,-2,-1,0,1,2,3), labels=c(-3,-2,-1,"XB",1,2,3))
axis(2)
title(main="Ordered probit for 5-category DV")
abline(v=c(-2.5,-1.5,-.5,2), lty=2)
text(-2.25,.3, "cut1")
text(-1.25,.3, "cut2")
text(-.25,.3, "cut3")
text(2.25,.3, "cut4")

curve(pnorm(x), from=-3, to=3, ylab="prob", xlab="", axes=F)
axis(1, at=c(-3,-2,-1,0,1,2,3), labels=c(-3,-2,-1,"XB",1,2,3))
axis(2)
title(main="Ordered probit for 5-category DV")
abline(v=c(-2.5,-1.5,-.5,2), lty=2)
text(-2.25,.3, "cut1")
text(-1.25,.3, "cut2")
text(-.25,.3, "cut3")
text(2.25,.3, "cut4")

### identification example ###

# p(Y=1) for 4-cat DV with cutpoints of -1, 1, and 2, and B0 = -2 and B1 = 1 #

x <- c(-1,0,1)

pnorm(-1 - (-2 + 1*x))

# Now subtract 1 from both c1 and B0 #

pnorm(-2 - (-3 + 1*x))

# Now add 2 #

pnorm(1 - (0 + 1*x))


### MLE for ordered probit using bbmle ###

library(bbmle)

y <- as.factor(data$abort)
X <- model.matrix(lm(abort ~ age + male + black + hisp + income + rural + church + wrole, data=data))
X <- X[ ,2:ncol(X)] # note the exclusion of the intercept column here

LL_oprobit <- function(params,y,X){
  y_wide <- array(NA, c(length(y), length(levels(y))))
  for (i in 1:length(y)){
    for (j in 1:length(levels(y))){
      if (y[i] == j){y_wide[i,j] <- 1} else{y_wide[i,j] <- 0}
    }
  }
  B <- params[1:ncol(X)]
  c1 <- params[ncol(X)+1]
  c2 <- params[ncol(X)+2]
  c3 <- params[ncol(X)+3]
  p1 <- pnorm(c1 - X %*% B)
  p2 <- pnorm(c2 - X %*% B) - pnorm(c1 - X %*% B)
  p3 <- pnorm(c3 - X %*% B) - pnorm(c2 - X %*% B)
  p4 <- 1 - pnorm(c3 - X %*% B)
  minusll  = -sum(y_wide[,1]*log(p1) + y_wide[,2]*log(p2) + y_wide[,3]*log(p3) + y_wide[,4]*log(p4)) #looks very similar to binomial, just add additional categories of DV
  return(minusll)
}

parnames(LL_oprobit) <- c("age","male","black","hisp","income","rural","church","wrole","c1","c2","c3")

fit1 <- mle2(LL_oprobit, start = c(age=0, male=0, black=0, hisp=0, income=0, rural=0, church=0, wrole=0, c1=0, c2=.5, c3=1),
             data=list(y=y, X=X), vecpar = TRUE)
summary(fit1)


### MLE for ordered logit using bbmle ###

y <- as.factor(data$abort)
X <- model.matrix(lm(abort ~ age + male + black + hisp + income + rural + church + wrole, data=data))
X <- X[ ,2:ncol(X)] # note the exclusion of the intercept column here

LL_ologit <- function(params,y,X){
  y_wide <- array(NA, c(length(y), length(levels(y))))
  for (i in 1:length(y)){
    for (j in 1:length(levels(y))){
      if (y[i] == j){y_wide[i,j] <- 1} else{y_wide[i,j] <- 0}
    }
  }
  B <- params[1:ncol(X)]
  c1 <- params[ncol(X)+1]
  c2 <- params[ncol(X)+2]
  c3 <- params[ncol(X)+3]
  p1 <- plogis(c1 - X %*% B)
  p2 <- plogis(c2 - X %*% B) - plogis(c1 - X %*% B)
  p3 <- plogis(c3 - X %*% B) - plogis(c2 - X %*% B)
  p4 <- 1 - plogis(c3 - X %*% B)
  minusll  = -sum(y_wide[,1]*log(p1) + y_wide[,2]*log(p2) + y_wide[,3]*log(p3) + y_wide[,4]*log(p4))
  return(minusll)
}

parnames(LL_ologit) <- c("age","male","black","hisp","income","rural","church","wrole","c1","c2","c3")

fit2 <- mle2(LL_ologit, start = c(age=0, male=0, black=0, hisp=0, income=0, rural=0, church=0, wrole=0, c1=0, c2=.5, c3=1),
             data=list(y=y, X=X), vecpar = TRUE)
summary(fit2)


### polr for oprobit and ologit (requires arm) ###

data$abort <- as.factor(data$abort) # polr requires that the DV be declared as a factor variable

fit3 <- polr(abort ~ age + male + black + hisp + income + rural + church + wrole, data=data, method="probit")
display(fit3)

fit4 <- polr(abort ~ age + male + black + hisp + income + rural + church + wrole, data=data, method="logistic") # note "logistic" not "logit" as with glm
display(fit4)


# Generate predicted probabilities of membership in each category for all observations #

pp4 <- predict(fit4, type="probs")
pp4[1:10,]


# Generate predicted probs across women's role var using average value/hold all else constant approach #

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


# What does this look like if we do it manually? #

X <- as.matrix(newdata)

p1 <- plogis(fit4$zeta[1] - X %*% coef(fit4)) #1st category
p2 <- plogis(fit4$zeta[2] - X %*% coef(fit4)) - plogis(fit4$zeta[1] - X %*% coef(fit4)) # 2nd cat
p3 <- plogis(fit4$zeta[3] - X %*% coef(fit4)) - plogis(fit4$zeta[2] - X %*% coef(fit4)) # 3rd cat
p4 <- plogis(-fit4$zeta[3] + X %*% coef(fit4)) # 4th cat

table <- cbind(p1,p2,p3,p4)
table
pp_wrole


## Using Zelig to use simulation to get CIs, etc. ##

library(Zelig)
library(ZeligChoice) # Add-on package for categorical vars

# Old Zelig approach #

z1 <- zelig(abort ~ age + male + black + hisp + income + rural + church + wrole, data=data, model="oprobit")

x1 <- setx(z1, wrole=c(1:7))

s1 <- sim(z1, x=wr1) # Get predicted probs across all 7 values of wrole

x2 <- setx(z1, wrole=1)
x3 <- setx(z1, wrole=7)

s2 <- sim(z1, x=x2, x1=x3) #Get change in probabilities moving from wrole = 1 to 7
s2

# New language in Zelig 5 #

z2 <- zoprobit$new()

z2$zelig(abort ~ age + male + black + hisp + income + rural + church + wrole, data=data)

z2$setx(wrole=1)
z2$setx1(wrole=7)

z2$sim()
summary(z2)


## Using the observed value/averaging over the data approach ##

detach("package:ZeligChoice", unload = T)

sim4 <- sim(fit4, n.sims=1000)
coef(sim4)[1,]

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
    l1[i, ,j,1] <- plogis(coef(sim4)[i,1] - X[,,j] %*% coef(sim4)[i,4:ncol(coef(sim4))]) 
    l1[i, ,j,2] <- plogis(coef(sim4)[i,2] - X[,,j] %*% coef(sim4)[i,4:ncol(coef(sim4))]) - plogis(coef(sim4)[i,1] - X[,,j] %*% coef(sim4)[i,4:ncol(coef(sim4))]) 
    l1[i, ,j,3] <- plogis(coef(sim4)[i,3] - X[,,j] %*% coef(sim4)[i,4:ncol(coef(sim4))]) - plogis(coef(sim4)[i,2] - X[,,j] %*% coef(sim4)[i,4:ncol(coef(sim4))]) 
    l1[i, ,j,4] <- plogis(-coef(sim4)[i,3] + X[,,j] %*% coef(sim4)[i,4:ncol(coef(sim4))]) 
    
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



## Marginal effects for Y > m using observed value approach ##

l1 <- array(NA, c(1000, length(X[,1,1]), 7))
l2 <- array(NA, c(1000, 7))
l3 <- array(NA, c(7,3))

for (j in 1:7){
  for (i in 1:1000){
    l1[i, ,j] <- plogis(-coef(sim4)[i,1] + X[,,j] %*% coef(sim4)[i,4:ncol(coef(sim4))]) # p(Y > 1) [replace coef(sim4)[i,1] to get different prob]
    l2[i,j] <- mean(l1[i,,j])
  }
  l3[j,1] <- mean(l2[ ,j])
  l3[j,2:3] <- quantile(l2[ ,j], probs=c(.025,.975))
}
  
l3

## Using post function ##

source("C:/Program Files/R/R-3.2.3/My Functions/post.R")

# Predicted probabilities and effect of church, 0 to 1 #

p1 <- post(fit3, "church", c(0,1))
p1$est

# p(Y > 2) #

p2 <- post(fit3, "wrole", c(1,7), "church", c(0,1), cut=2)
p2$est

























































