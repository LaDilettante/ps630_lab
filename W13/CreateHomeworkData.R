# Script to create the data for the RDD
# SocialSecurity

set.seed(2015)
x=rnorm(1000,mean=50,sd=10)
set.seed(2)
y=rep(100,1000)+2.5*x+rnorm(1000,mean=0,sd=10)

summary(x)
summary(y)

for (i in 1:1000){
  if (x[i] >= 50.00){
    y[i] = y[i] + 20
  }
}

SocialSecurity=cbind(x,y)
colnames(SocialSecurity)=c("LibPartyVoteShare","SocialSecurityExp")
SocialSecurity=as.data.frame(SocialSecurity)

save(SocialSecurity, file="SocialSecurity.Rdata")