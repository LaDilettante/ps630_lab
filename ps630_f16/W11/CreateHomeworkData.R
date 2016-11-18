##########################################
### Script to create the data for the RDD
##########################################

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



##########################################
### Script to create the data for the DiD
##########################################

set.seed(2015)
year1=rep(2015,1000)

set.seed(2016)
year2=rep(2016,1000)

set.seed(2015)
restaurant_owner_2015=rbinom(1000, size=1, prob=0.5)

set.seed(2016)
restaurant_owner_2016=rbinom(1000, size=1, prob=0.5)

set.seed(2015)
private_expenditures_2015=rnorm(1000,mean=25000,sd=5000)

set.seed(2016)
private_expenditures_2016=rnorm(1000,mean=27500,sd=5000) + restaurant_owner_2016 * 2000

Data2015=cbind(year1,restaurant_owner_2015,private_expenditures_2015)
colnames(Data2015)=c("Year","RestOwn","PrivExp")
Data2015=as.data.frame(Data2015)
save(Data2015, file="Data2015.Rdata")

Data2016=cbind(year2,restaurant_owner_2016,private_expenditures_2016)
colnames(Data2016)=c("Year","RestOwn","PrivExp")
Data2016=as.data.frame(Data2016)
save(Data2016, file="Data2016.Rdata")


