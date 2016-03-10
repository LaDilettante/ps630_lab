#################################
### Tutorial 8 - Count Models ###
#################################

##### Count DVs #####

data <- read.table("wilsonpiazzavars.txt", header=TRUE)
summary(data)

#DV:
#FGTDDom: count of domestic terrorist incidents in a given country in a given year
#IVs:
#chga_dem: 1 = democracy, 0 otherwise
#geddes1: 1 = personalist auth regimes, 0 otherwise
#geddes2: 1 = military auth regimes, 0 otherwise [excluded category]
#geddes3: 1 = single party-based auth regimes, 0 otherwise
#geddes4: 1 = monarchies, 0 otherwise
#geddes5: 1 = hybrid (party-military-personalist) auth regimes, 0 otherwise
#logpop = country's population
#logarea = country's area
#logGNI = gross national income per capita
#GINI = inequality measured by the Gini coefficient on a 0-100 scale
#ColdWar: 1 if 1970-1991, 0 otherwise
#Durable = durability score from Polity IV


# Introduction #

hist(data$FGTDDom, main="Count of domestic terrorist attacks in year", xlab="")

hist(rpois(1000, 1),xlab="Count", ylab="Density", main="mu = 1")
hist(rpois(1000, 5),xlab="Count", ylab="Density", main="mu = 5")
hist(rpois(1000, 50),xlab="Count", ylab="Density", main="mu = 50")


### Using GLM ###

library(arm)

fit4 <- glm(FGTDDom ~ chga_dem + geddes1 + geddes3 + geddes4 + geddes5 
            + logpop + logGNI + logarea + ColdWar, data=data, family=poisson(link="log"))
display(fit4)

exp(coef(fit4))


# Use predict with type=response to get expected count #

predict(fit4, type="response")

newdata <- data.frame(chga_dem=c(0,1), geddes1=c(0,0), geddes3=c(0,0), geddes4=c(0,0), geddes5=c(0,0),
                      logpop=c(rep(mean(data$logpop),2)), logarea=c(rep(mean(data$logarea),2)),
                      logGNI=c(rep(mean(data$logGNI),2)), ColdWar=c(1,1))

p1 <- predict(fit4, type="response", newdata=newdata)
p1

# Plot distribution of terrorist attacks per year given expected count for military v. democracy #

plot(0:30, dpois(0:30,p1[1]), type="p", pch=19, 
     ylab="p(Y=y)", xlab="# of terrorist attacks in a given year", xlim=c(0,30))
lines(0:30, dpois(0:30,p1[2]), type="p", pch=1)
lines(0:30, dpois(0:30,p1[1]), lty=2)
lines(0:30, dpois(0:30,p1[2]), lty=2)
legend("topright", c("Military Dictatorships","Democracies"), pch=c(19,1), bty="n", inset=.05)
title(main="Distribution of terrorist attacks for military and democratic regimes")


### Probs for specific counts ###

dpois(5, p1[1])
dpois(10, p1[1])
dpois(15, p1[1])

dpois(5, p1[2])
dpois(10, p1[2])
dpois(15, p1[2])

# Plot p(Y = 5,10,15) as a function of population in country with military regime *

summary(data)

lpop <- c(1,2,3,4,5)
p5 <- dpois(5, exp(1.19 + .46*lpop - .05*7 - .03*12 + .64*1))
p10 <- dpois(10, exp(1.19 + .46*lpop - .05*7 - .03*12 + .64*1))
p15 <- dpois(15, exp(1.19 + .46*lpop - .05*7 - .03*12 + .64*1))

plot(1:5, p5, type="l", ylim=c(0,.2), xlim=c(.75,5.25), xlab="Size of population", ylab="p(Y=y)")
lines(1:5, p10, lty=2)
lines(1:5, p15, lty=3)
legend("topright", c("Y=5","Y=10","Y=15"), lty=c(1,2,3), bty="n", inset=.05)
title(main="Change in p(Y=y) for military regimes as a function of population size")


### Using Zelig ###

library(Zelig)

z1 <- zelig(FGTDDom ~ chga_dem + geddes1 + geddes3 + geddes4 + geddes5 
            + logpop + logGNI + logarea + ColdWar, data=data, model="poisson")
summary(z1)

x <- setx(z1, chga_dem=0, geddes1=0, geddes3=0, geddes4=0, geddes5=0, ColdWar=1)
x1 <- setx(z1, chga_dem=1, geddes1=0, geddes3=0, geddes4=0, geddes5=0, ColdWar=1)

s1 <- sim(z1, x=x, x1=x1)
s1


### Using sim in arm ###

detach("package:Zelig", unload=T)
library(arm)

sim4 <- sim(fit4, n.sims=1000)
coef(sim4)[1,]

## expected count for democracies and military regimes averaging over data ##

X1 <- model.matrix(glm(FGTDDom ~ chga_dem + geddes1 + geddes3 + geddes4 + geddes5 
                       + logpop + logGNI + logarea + ColdWar, data=data, family=poisson(link="log")))
X2 <- model.matrix(glm(FGTDDom ~ chga_dem + geddes1 + geddes3 + geddes4 + geddes5 
                       + logpop + logGNI + logarea + ColdWar, data=data, family=poisson(link="log")))
X1[1,]

# Military regimes
X1[,"chga_dem"] <- 0
X1[,"geddes1"] <- 0
X1[,"geddes3"] <- 0
X1[,"geddes4"] <- 0
X1[,"geddes5"] <- 0

ec_mil <- apply(apply(X1, 1, function (x) exp(coef(sim4) %*% x)), 1, mean)
mean(ec_mil)
sd(ec_mil)
quantile(ec_mil, probs=c(.025,.975))

# Dem regimes
X2[,"chga_dem"] <- 1
X2[,"geddes1"] <- 0
X2[,"geddes3"] <- 0
X2[,"geddes4"] <- 0
X2[,"geddes5"] <- 0

ec_dem <- apply(apply(X2, 1, function (x) exp(coef(sim4) %*% x)), 1, mean)
mean(ec_dem)
sd(ec_dem)
quantile(ec_dem, probs=c(.025,.975))


### Using post ###

source(file.choose())

# Same thing as directly above #

out1 <- post(fit4, "chga_dem", c(0,1), holds=list(geddes1=0, geddes3=0, geddes4=0, geddes5=0))
out1$est

# All else constant at means approach #

out2 <- post(fit4, "chga_dem", c(0,1), 
             holds=list(geddes1=0, geddes3=0, geddes4=0, geddes5=0, ColdWar=1,
                        logpop=mean(data$logpop), logGNI=mean(data$logGNI), logarea=mean(data$logarea)))
out2$est

# Calculate effects of key variables and plot them #

out3a <- post(fit4, "ColdWar", c(0,1))
out3b <- post(fit4, "logpop", c(quantile(data$logpop, probs=c(.05,.95))))
out3c <- post(fit4, "logGNI", c(quantile(data$logGNI, probs=c(.05,.95))))
out3d <- post(fit4, "logarea", c(quantile(data$logarea, probs=c(.05,.95))))


plot(1:4, c(out3a$est[3,1],out3b$est[3,1],out3c$est[3,1],out3d$est[3,1]), pch=19,
     xlim=c(.75,4.25), ylim=c(-5,25), xlab="", ylab="Change in expected count", axes=F)
axis(1, at=c(1:4), labels=c("Cold War","Population","GNI","Area"))
axis(2, at=c(-5,0,5,10,15,20,25))
segments(1:4, c(out3a$est[3,2],out3b$est[3,2],out3c$est[3,2],out3d$est[3,2]),
         1:4, c(out3a$est[3,3],out3b$est[3,3],out3c$est[3,3],out3d$est[3,3]))
abline(h=0, lty=3)
title(main="Effects of predictors on expected count of terrorist attacks")

# Vary two vars at once #

out4 <- post(fit4, "chga_dem", c(0,1), "ColdWar", c(0,1), holds=list(geddes1=0, geddes3=0, geddes4=0, geddes5=0))
out4$est


### Including exposure via the offset parameter ###

data_wide <- read.table("wilsonpiazza_wide.txt", header=TRUE)
summary(data_wide)

# Without offset #

fit5 <- glm(tcount ~ as.factor(tdem) + avgpop, data=data_wide, family=poisson(link="log"))
display(fit5)

# With offset #

fit6 <- glm(tcount ~ as.factor(tdem) + avgpop, data=data_wide, family=poisson(link="log"), offset=log(years))
display(fit6)


## Exposure for police stops (baseline category is black, 2=hisp, 3=white) ##
# Using offset tells us about # of stops relative to rates of arrest (stops per arrest) #
# Essentially we are normalizing the stop rates across groups to their rates of arrest,
# and then examining whether some groups are stopped disproprtionately given arrest rates #

data2 <- read.table("C:/Users/cdj19/Documents/Teaching/Graduate/Classes/MLE/Spring 2016/Data/stopfrisk.txt", header=TRUE)
summary(data2)

fit6 <- glm(stops ~ factor(eth), data=subset(data2, past.arrests>0), family=poisson(link="log"),
            offset=log(past.arrests))
display(fit6)

# Add fixed effects (dummies) for precinct #

fit7 <- glm(stops ~ factor(eth) + factor(precinct), data=subset(data2, past.arrests>0), family=poisson(link="log"),
            offset=log(past.arrests))
display(fit7)