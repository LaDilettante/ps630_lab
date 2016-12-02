###########################
#
# RD analysis: Replication of Fouirnaies and Hall (2012), ``The Financial Incumbency Advantage"
# Last update: Nov 10 2015
# RD Software Website:  https://sites.google.com/site/rdpackages/
# Illustration for "A Practical Guide to Regression Discontinuity Designs in Political Science", by Christopher Skovron and Rocio Titiunik
#########################

library(rdrobust)
library(foreign)
library(ggplot2)
library(plyr)
library(dplyr)
source("rddensity_fun.R")
source("rdbwdensity.R")
source("rddensity.R")

set.seed(48104)

# Set your directory here, or just use your current directory
## setwd('~/Dropbox/RD-review/data-examples/Fouirnaies and Hall')

############################
#
# Open data
#
###########################
dat = read.dta('fouirnaies_hall_financial_incumbency_advantage.dta')
names(dat)
table(dat$statelevel)

# Filter to only state legislative level - just for illustrative purposes in our example
dat = filter(dat, statelevel==1)

# summarize running variable: vote share
summary(dat$rv)

# summarize outcome variable
summary(dat$dv_money)

# visualize outcome variable 
p = qplot(rv,dv_money,data=dat)+xlab("Democratic margin of victory at t")+ylab("Democratic share of contributions at t+1")
p

############################
#
#  Continuity-based analysis
#
###########################

#############
# Validation 
##############
#  Density test 
rddensity(X = dat$rv, vce="jackknife")

# histogram of density test
p = ggplot(dat,aes(x=rv, fill = factor(dat$rv>0)))+geom_histogram(binwidth=0.5)+xlim(-25,25)+geom_vline(xintercept = 0)+xlab("Democratic vote share at t")+scale_colour_manual(values = c("red","blue"))+theme_bw()+theme(legend.position='none')
p
ggsave("mccrary.pdf")

# Covariate as outcome: total money in race
rdrobust(dat$total_race_money,dat$rv,all=TRUE)

pdf('total-money.pdf')
rdplot(dat$total_race_money,dat$rv,x.lim = c(-10,10),y.lim = c(0,400000),x.lab="Democratic margin of victory at t",y.lab="Total money in race at t+1", title = "")
dev.off()

# Covariate as outcome: total votes in race
rdrobust(dat$total_votes,dat$rv,all=TRUE)

pdf('total-votes.pdf')
rdplot(dat$total_votes,dat$rv,x.lim = c(-10,10),y.lim = c(0,50000),x.lab="Democratic margin of victory at t",y.lab="Total votes in race at t+1", title = "")
dev.off()

# Covariate as outcome: Democratic incumbent
rdrobust(dat$dem_inc,dat$rv,all=TRUE)

pdf('dem-inc.pdf')
rdplot(dat$dem_inc,dat$rv,x.lim = c(-10,10),y.lim = c(0,1),x.lab="Democratic margin of victory at t",y.lab="Democratic incumbent at t", title = "")
dev.off()

# Covariate as outcome: Republican incumbent
rdrobust(dat$rep_inc,dat$rv,all=TRUE)

pdf('rep-inc.pdf')
rdplot(dat$rep_inc,dat$rv,x.lim = c(-10,10),y.lim = c(0,1),x.lab="Democratic margin of victory at t",y.lab="Republican incumbent at t", title='')
dev.off()

# Covariate as outcome: Total group money 
rdrobust(dat$total_group_money,dat$rv,all=TRUE)

pdf('total-group-money.pdf')
rdplot(dat$total_group_money,dat$rv,x.lim = c(-10,10),y.lim = c(0,200000),x.lab="Democratic margin of victory at t",y.lab="Total group money in race at t+1", title='')
dev.off()

#############
#  Estimation and inference for RD effect
##############
# RD Effect for main outcome variable ==> Local linear polynomial estimation with optimal bandwidth
rdrobust(dat$dv_money,dat$rv,all=TRUE)
# optimal h is 8.9743

# Fouirnaies and Hall use local linear and bandwidths of 1, 2, and 3 percentage points. This is what that looks like.
rdrobust(dat$dv_money,dat$rv,h=1,all=TRUE)
rdrobust(dat$dv_money,dat$rv,h=2,all=TRUE)
rdrobust(dat$dv_money,dat$rv,h=3,all=TRUE)

# RD Effect for main outcome variable ==> Local linear polynomial estimation with optimal bandwidth
rdrobust(dat$dv_money,dat$rv, p = 2, q =3, all=TRUE)

# MSE-optimal bandwidth for local linear estimation 
h = rdbwselect(dat$dv_money,dat$rv, p =1)$bws[1]

# Run local quadratic estimation using the optimal h from the local linear estimated above==> 
# note that conventional local quadratic estimation 
rdrobust(dat$dv_money,dat$rv, p = 2, q =3, h = h, rho = 1, all=TRUE)
# is equivalent to robust local linear estimation 
rdrobust(dat$dv_money,dat$rv, p = 1, rho = 1, all=TRUE)

#  RD plot
# CCT bandwidth in main estimation was about 9, restricting x axis  to (-10,10)
pdf('estimate-10-points.pdf')
rdplot(dat$dv_money,dat$rv,x.lim=c(-10,10),x.lab="Democratic margin of victory at t",y.lab="Democratic share of contributions at t+1", title = "RDD Estimate")
dev.off()

# RD plot using whole data range
pdf('estimate-entire-range.pdf')
rdplot(dat$dv_money,dat$rv,x.lab="Democratic margin of victory at t",y.lab="Democratic share of contributions at t+1", title = "RDD Estimate")
dev.off()

#########################
# Placebo cutoffs
##########################

# c = 1
rdrobust(dat$dv_money,dat$rv, c = 1 ,all=TRUE)

# c = -3
rdrobust(dat$dv_money,dat$rv, c = -3 ,all=TRUE)

############################
#
#  Randomization-based analysis: Analysis in the paper is implemented in Stata. R implementation coming soon.#
#  See Stata code fouirnaies-hall-replication-final.do
###########################

