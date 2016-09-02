######## Homework #5 #######

data <- read.delim(
  "W7/2012 ANES PID.txt", header=T)
summary(data)


### 1 ###

library(arm)
source("post.R")

data$tpid <- as.factor(data$tpid)

fit1 <- polr(tpid ~ south + male + union + hisp + black + income01, data=data, method="logistic")
summary(fit1)

p1 <- post(fit1, "income01", c(quantile(data$income, probs=c(.05,.95))))
p1$est


library(VGAM)

fit2 <- vglm(ordered(tpid) ~ south + male + union + hisp + black + income01,
             cumulative(parallel = F, reverse=T, link = "logit"), data=data)
summary(fit2)

fit3 <- vglm(ordered(tpid) ~ south + male + union + hisp + black + income01,
             cumulative(parallel = T, reverse=T, link = "logit"), data=data)
summary(fit3)

lrtest(fit3,fit2)


fit4 <- vglm(ordered(tpid) ~ south + male + union + hisp + black + income01,
             cumulative(parallel = F~(income01), reverse=T, link = "logit"), data=data)
summary(fit4)

lrtest(fit4,fit2)


### 2 ###

detach("package:arm", unload=T)
detach("package:VGAM", unload=T)

library(Zelig)
library(ZeligChoice)


z1 <- zelig(tpid ~ south + male + union + hisp + black + income01, data=data, model="mlogit")
summary(z1)

x <- setx(z1, south=0, male=0, union=0, hisp=0, black=0, income01=.44)

s1 <- sim(z1, x=x)
s1

data$tpid <- as.numeric(data$tpid)
fit2 <- vglm(tpid ~ south + male + union + hisp + black + income01,
             data=data, multinomial)
# PCP PRE
Y.p <- predict(fit2, type="response")
Y.hat <- apply(Y.p, 1, function(x) which.max(x))
Y <- data$tpid

cp <- 0
for (i in 1:length(Y)){
  if (Y.hat[i] == Y[i]){cp <- cp + 1}
}
pcp <- cp/length(Y)
pcp

pre <- (pcp - max(table(Y)/length(Y)))/(1-max(table(Y)/length(Y)))
pre

# South #
x <- setx(z1, south=0, male=0, union=0, hisp=0, black=0, income01=.44)
x1 <- setx(z1, south=1, male=0, union=0, hisp=0, black=0, income01=.44)
s2 <- sim(z1, x=x, x1=x1)
s2

# Male #
x <- setx(z1, south=0, male=0, union=0, hisp=0, black=0, income01=.44)
x1 <- setx(z1, south=0, male=1, union=0, hisp=0, black=0, income01=.44)
s3 <- sim(z1, x=x, x1=x1)
s3

# Union #
x <- setx(z1, south=0, male=0, union=0, hisp=0, black=0, income01=.44)
x1 <- setx(z1, south=0, male=0, union=1, hisp=0, black=0, income01=.44)
s4 <- sim(z1, x=x, x1=x1)
s4

# Hispanic #
x <- setx(z1, south=0, male=0, union=0, hisp=0, black=0, income01=.44)
x1 <- setx(z1, south=0, male=0, union=0, hisp=1, black=0, income01=.44)
s5 <- sim(z1, x=x, x1=x1)
s5

# Black #
x <- setx(z1, south=0, male=0, union=0, hisp=0, black=0, income01=.44)
x1 <- setx(z1, south=0, male=0, union=0, hisp=0, black=1, income01=.44)
s6 <- sim(z1, x=x, x1=x1)
s6

# Income #
x <- setx(z1, south=0, male=0, union=0, hisp=0, black=0, income01=quantile(data$income01, probs=.05))
x1 <- setx(z1, south=0, male=0, union=0, hisp=0, black=0, income01=quantile(data$income01, probs=.95))
s7 <- sim(z1, x=x, x1=x1)
s7


detach("package:ZeligChoice", unload=T)

data$tpid <- data$tpid - 1
fit5 <- glm(tpid ~ south + male + union + hisp + black + income01,
            data=subset(data, tpid==0 | tpid==2), family=binomial(link="logit"))
summary(fit5)




































