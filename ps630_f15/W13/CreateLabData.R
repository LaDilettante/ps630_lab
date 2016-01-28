# Script to create the data for the RDD for the lab
# Medicare

set.seed(2013)
x=rnorm(1000,mean=65,sd=10)
y=rep(800,1000)+4*x
plot(x,y1)

set.seed(2)
y=y+rnorm(1000,mean=0,sd=200)
plot(x,y)

summary(x)
summary(y)

for (i in 1:1000){
  if (x[i] >= 65.00){
    y[i] = y[i] + 350 + rnorm(1,mean=0,sd=100)
  }
}

x=round(x)

set.seed(3)
a=rnorm(1000,mean=100,sd=10)
set.seed(1000)
b=rnorm(1000,mean=100,sd=10)

Medicare=cbind(x,y,a,b)
colnames(Medicare)=c("Age","MedicalExp","covariate1","covariate2")
Medicare=as.data.frame(Medicare)

save(Medicare, file="Medicare.Rdata")