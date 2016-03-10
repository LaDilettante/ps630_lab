################## Binary DVs - Lecture 4 #####################

library(arm)

data1 <- read.table(
  "swingvoterafrica.txt", header=TRUE)
summary(data1)


#1268 non-missing respondent surveys in Ghana in 2008, four months prior to a presidential and legislative election. 
#Paper attempts to understand why some voters are swing voters rather than "core," 
#stable party supporters. Key hypotheses contrast MP performance versus clientelism as competing explanations. 

#Variables:

#'tsv_binary': core DV; binary, '1'=swing voter, '0'=core party supporter
#'lawmaking': binary, '1'= respondent views MP's performance in Parliament (e.g., passing bills) positively 
#'cliesup': 0-5, clientelism supply, degree of exposure to small gifts from candidates ahead of elections 
#'age': 0-3, categories are 18-22, 23-35, 36-55, and 55 and up
#'male': '1' = male, '0' = female
#'edulevel': 0-4 (effectively 0-3), from no formal schooling to post-secondary education
#'party_member': binary, '1'= active in political party, '0'= not active 


###### Some Basic Diagnostics ######

#Estimate a simple model

fit1 <- glm(tsv_binary ~ lawmaking + cliesup + age + male + edulevel + party_member, data=data1, family=binomial(link = "probit"))
display(fit1)

### Residuals (see Gelman chapter on logit, p. 97) ###

#Basically, we plot average residuals within groups defined by portions of the x-axis 
#(which is the predicted probs or individual predictors). We then look for patterns

# Calculate predicted probs for model #

pp1 <- predict(fit1, type="response")

#Generate residuals, observed Y minus expected value of Y (i.e., the predicted probs)

p_res1 <- data1$tsv_binary - pp1

#You can see that if we just plot the residuals versus the predicted, this is not helpful

plot(pp1,p_res1, xlab="Predicted probabilities", ylab="Response residuals")
abline(h=0, lty=2)

#The 'arm' package contains a function for plotting binned residuals
#Grey lines are 95% bounds, within which 95% of binned resids should fall
#The second command fits a lowess (or loess) line to the residuals; lowess lines are locally weighted lines of best fit
#Essentially, they "follow the relationship" between x and y across the range of x, even when it is non-linear
#This allows us to see if any patterns emerge in the residuals, as below, one does which suggests at least some non-linearity

binnedplot(pp1,p_res1)
lines(lowess(pp1,p_res1), lwd=2)

#We can also examine plots of binned residuals against individual predictors, 
#looking for problems with each (again, see Gelman p. 98 for another e.g.)

#Education looks pretty good

binnedplot(data1$edulevel, p_res1)
lines(lowess(data1$edulevel,p_res1), lwd=2)

#Cliesup, not so much...look at the prominant pattern for the second category

binnedplot(data1$cliesup, p_res1)
lines(lowess(data1$cliesup,p_res1), lwd=2)

#We can examine this further by looking at the distribution of cliesup itself
#Clearly this variable has a weird distribution, 
#where most respondents have a value of zero and only a few have a value of 1

hist(data1$cliesup)

#The distribution of this variable suggests that it is misoperationalized. 
#That is, the influence of the variable if probably best modeled
#as some sort of categorical distinction. This could be done in multiple ways, 
#e.g., yes or no to clientalistic benefits, or something like that you could, perhaps, 
#just dummy out all categories to allow for full non-linearity. 
#Ultimately, you would have to make a judgment call on your own


## Other Residuals ##

# Pearson residuals #

pearson <- residuals(fit1, type="pearson")
binnedplot(data1$cliesup, pearson)
lines(lowess(data1$cliesup,pearson), lwd=2)

# Deviance residuals #

dev <- residuals(fit1, type="deviance")
sum(dev^2)
display(fit1)

binnedplot(data1$cliesup, dev)
lines(lowess(data1$cliesup,dev), lwd=2)


### Influence measures ###

library(car)

#point 181 looks potentially problematic

influenceIndexPlot(fit1, vars=c("Cook","hat","Studentized"), id.n=3)

#Another way to look at influence, again 181 has the largest Cook's

influencePlot(fit1, id.n=3)

#let's remove 181 and see if things change meaningfully
#This command will compare coefficients from the original model with an "updated" model that removes ('-c()') point 181

compareCoefs(fit1, update(fit1, subset=-c(181)))

#Things look the same after removing 181, what if we remove the top 3 influential points?
#Still fine, so I think we are good

compareCoefs(fit1, update(fit1, subset=-c(181,216,1091)))


### Simulating data to check model fit ###

display(fit1)

# Simulate the coefficients from the posterior #

sim1 <- sim(fit1, n.sims=1000)
coef(sim1)[1,]

# Generate 1000 new datasets #

n <- length(data1$tsv_binary)
X <- cbind(rep(1,n), data1$lawmaking, data1$cliesup, data1$age, data1$male, data1$edulevel, data1$party_member)
y_rep <- array(NA, c(1000, n))

for (i in 1:1000){
  p <- pnorm(X %*% coef(sim1)[i, ])
  y_rep[i, ] <- rbinom(n, 1, p)
}


plot(data1$cliesup, data1$tsv_binary)
lines(lowess(data1$cliesup, data1$tsv_binary), lwd=2)
for (i in 1:1000){
  lines(lowess(data1$cliesup, y_rep[i, ]), lty=2, col="grey")
}

hist(data1$cliesup)


##### Other Binomial Models #####

### normal versus t ###

curve(dnorm(x), from=-5, to=5, ylab="")
curve(dt(x, 3), from=-5, to=5, add=T, lty=2)
legend("topright", c("normal","t(3)"), lty=c(1,2), bty="n", inset=.05)


### Simulated example with data drawn from normal distribution but with a few outliers ###

set.seed(2)
simdata <- array(NA, c(1000,2))
for (i in 1:1000){
  simdata[i,2] <- rnorm(1,0,1)
  if ( (0 + 1*simdata[i,2] + rnorm(1,0,1)) > 0){simdata[i,1] <- 1} else{simdata[i,1] <- 0}
}

# Add some bad outliers to the data #

simdata[1,1] <- 1
simdata[1,2] <- -5
simdata[2,1] <- 0
simdata[2,2] <- 5


plot(simdata[ ,2], simdata[ ,1])


y1 <- simdata[,1]
X1 <- model.matrix(lm(simdata[,1] ~ simdata[,2]))


## Probit ##

fit2 <- glm(y1 ~ X1 - 1, family=binomial(link="probit"))
display(fit2)

compareCoefs(fit2, update(fit2, subset=-c(1,2)))


## Robit Regression ##

library(bbmle)

LL_robit_3 <- function(params,y,X){
  B <- params
  p <- pt(X %*% B, 3) #t link w/ 3 df
  minusll  = -sum(y*log(p) + (1-y)*log(1-p))
  return(minusll)
}

parnames(LL_robit_3) <- c("B0", "B1")

fit3a <- mle2(LL_robit_3, start = c(B0=0, B1=0),
             data=list(y=y1,X=X1), vecpar = TRUE)

fit3b <- mle2(LL_robit_3, start = c(B0=0, B1=0),
              data=list(y=y1[3:length(y1)],X=X1[3:length(y1),]), vecpar = TRUE)


summary(fit3a)
summary(fit3b)


plot(X1[ ,2], y1, xlim=c(-5,5))
curve(pnorm(0 + 1*x), from=-5, to=5, add=T)
curve(pnorm(coef(fit2)[1] + coef(fit2)[2]*x), from=-5, to=5, add=T, lty=2)
legend("topleft", c("True model","Probit fit"), lty=c(1,2), bty="n", inset=.05)


plot(X1[ ,2], y1, xlim=c(-5,5))
curve(pnorm(0 + 1*x), from=-5, to=5, add=T)
curve(pt(coef(fit3a)[1] + coef(fit3a)[2]*x, 3), from=-5, to=5, add=T, lty=2)
legend("topleft", c("True model","Robit fit"), lty=c(1,2), bty="n", inset=.05)


mean(abs(pnorm(coef(fit2)[1] + coef(fit2)[2]*X1[,2]) - pnorm(0 + 1*X1[,2])))
mean(abs(pt(coef(fit3a)[1] + coef(fit3a)[2]*X1[,2], 3) - pnorm(0 + 1*X1[,2])))


### Robit w/ Inf df = Probit ###

LL_robit_Inf <- function(params,y,X){
  B <- params
  p <- pt(X %*% B, Inf) #t link w/ Inf df
  minusll  = -sum(y*log(p) + (1-y)*log(1-p))
  return(minusll)
}

parnames(LL_robit_Inf) <- c("B0", "B1")

fit4 <- mle2(LL_robit_Inf, start = c(B0=0, B1=0),
              data=list(y=y1,X=X1), vecpar = TRUE)

summary(fit4)
summary(fit2)



###Asymmetric links
#package 'VGAM' can be useful here, see documentation; but 'glm' can estimate cloglog

library(VGAM)

#cloglog

curve(cloglog(x, inverse=TRUE), from=-10, to=5, ylab="cumulative probability", xlab="linear predictor", main="CDF")
curve(1-exp(-exp(x)), from=-10, to=5, ylab="cumulative probability", xlab="linear predictor", main="CDF")

curve(log(-log(1-x)), from=0, to=1, ylab="linear predictor", xlab="p(Y=1)")


cll <- glm(tsv_binary ~ lawmaking + cliesup + age + male + edulevel + party_member, data=data1, family=binomial(link="cloglog"))
display(cll)


#Example calculation of predicted probs, hold all at zero and vary party_member

curve(cloglog(coef(cll)[1] + coef(cll)[7]*x, inverse=TRUE), 
      from=0, to=1, ylim=c(0,1), xlab="party member", ylab="pr(swing)")

p1 <- cloglog(coef(cll)[1] + coef(cll)[7]*0, inverse=T)
p2 <- cloglog(coef(cll)[1] + coef(cll)[7]*1, inverse=T)

p1
p2


####### Dealing w/ Separation #######

#There are a few steps you need to do here. You will use the function 'bayesglm' in the 'arm' package

#First, recode all binary inputs (only IVs) to have a mean = 0 and differ by 1 unit in lower and upper categories
#Second, recode all non-binary inputs to have a mean = 0 and an SD = .5

#This coding of the vars is set for the default prior distribution utilized by Bayes' GLM (i.e., a t with 1 df, i.e., the Cauchy)
#You can see that this prior distribution basically says: 'most coefs will be within 5 units of 0, with occasional large effects'

curve(dcauchy(x, location = 0, scale = 2.5), from=-10, to=10)

#You can compare to other versions with different scales; the larger the scale, the more diffuse the prior

curve(dcauchy(x, location = 0, scale = 2.5), from=-10, to=10)
curve(dcauchy(x, location = 0, scale = 5), from=-10, to=10, add=TRUE, lty=2)

#Note also, that the defaults are set for a logit link not a probit, so use logit

law.sd <- data1$lawmaking - mean(data1$lawmaking)
clie.sd <- (data1$cliesup - mean(data1$cliesup))/2*sd(data1$cliesup)
age.sd <- (data1$age - mean(data1$age))/2*sd(data1$age)
edu.sd <- (data1$edulevel- mean(data1$edulevel))/2*sd(data1$edulevel)
male.sd <- data1$male- mean(data1$male)
party.sd <- data1$party_member- mean(data1$party_member)


#Run the model; the default prior is the one you want most of the time, but you can adjust if needed

#Compare normal logit to Bayes logit

l1 <- glm(tsv_binary ~ law.sd + clie.sd + age.sd + male.sd + edu.sd + party.sd, data=data1, family=binomial(link = "logit"))

bl1 <- bayesglm(data1$tsv_binary ~ law.sd + clie.sd + age.sd + male.sd + edu.sd + party.sd, family=binomial(link = "logit"))

compareCoefs(l1,bl1)

##As you can see, things look very similar, as they should: there are no separation problems, and the priors are weak

#If you set a non-informative prior then the logit and the Bayes logit will be truly identical
#They are the exact same model, because the posterior in Bayes is proportional to the prior times the likelihood

bl2 <- bayesglm(data1$tsv_binary ~ law.sd + clie.sd + age.sd + male.sd + edu.sd + party.sd, 
                prior.scale=Inf, prior.df=Inf, family=binomial(link = "logit"))

compareCoefs(l1,bl2)

##What if we use informative priors? Things should change more. 

bl3 <- bayesglm(data1$tsv_binary ~ law.sd + clie.sd + age.sd + male.sd + edu.sd + party.sd, 
                prior.scale=.01, prior.df=1, family=binomial(link = "logit"))

compareCoefs(l1,bl3)

##As you can see, all of the estimated parameters in the Bayes model are now close to 0
 #This is because I placed a very strong prior with a mean of 0 on each coefficient


###Now let's create a problem with separation and then show how this method solves it

#We will set swing voter DV to '0' whenever party_member = 1, thus creating perfect separation

swing.new <- data1$tsv_binary
for (i in 1:length(data1$tsv_binary)){
  if (data1$party_member[i] == 1){swing.new[i] <- 0}
  if (data1$party_member[i] == 0){swing.new[i] <- 1}
}

#Demonstrate perfect separation

cor(swing.new, data1$party_member)

#Try to estimate logit with and without perfect separation, you will get an error message for lack of convergence in the latter

l1 <- glm(data1$tsv_binary ~ law.sd + clie.sd + age.sd + male.sd + edu.sd + party.sd, family=binomial(link = "logit"))

l2 <- glm(swing.new ~ law.sd + clie.sd + age.sd + male.sd + edu.sd + party.sd, family=binomial(link = "logit"))

#Fit the same model with weakly informative priors

bl2 <- bayesglm(swing.new ~ law.sd + clie.sd + age.sd + male.sd + edu.sd + party.sd, family=binomial(link = "logit"))
display(bl2)

## There are still major problems here obviously, but the model estimates. Ultimately, this is just a contrived example


###### Rare events logit (need package Zelig) and Zelig command 'relogit' ######

### Bias in intercept as function of sample size and average p(Y=1) ###

n <- c(500,1000,10000,100000)

curve((x-.5)/(100*x*(1-x)), from=0, to=1, ylab="bias in B0", xlab="average p(Y=1)")
for (i in 1:length(n)){
  curve((x-.5)/(n[i]*x*(1-x)), from=0, to=1, add=T)
}

curve((x-.5)/(100*x*(1-x)), from=0, to=1, xlim=c(0,.05), ylab="bias in B0", xlab="average p(Y=1)")
for (i in 1:length(n)){
  curve((x-.5)/(n[i]*x*(1-x)), from=0, to=1, add=T)
}

### Use relogit ###

library(Zelig)

data(mid)
summary(mid)

l.conflict <- glm(conflict ~ major + contig + power + maxdem + mindem + years, data=mid, family=binomial(link = "logit"))
rel.conflict <- zelig(conflict ~ major + contig + power + maxdem + mindem + years, data = mid, model="relogit")

summary(l.conflict)
summary(rel.conflict)



















































