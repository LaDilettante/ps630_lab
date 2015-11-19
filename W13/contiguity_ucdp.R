library(cshapes) # using cshapes 0.4-1
library(car)

setwd("~/Documents/Teaching/630-F15/data")



#load ucd intrastate conflict data
ucdp<-read.csv("ucdp_1onset2012.csv")
ucdp<-ucdp[ucdp$incidencev412==1,]

########form neighborhood connectivity matrices and associated data for 2005
#use the cshapes distmatrix function to form a matrix of distances in 2005
dmat <- distmatrix(as.Date("2005-1-1"), type="capdist")
#define neihbors as having distance less than 900km
adjmat <- ifelse(dmat > 900, 0, 1)
#get into an NXN matrix by using the diag function
diag(adjmat) <- 0
#adjmat is now the connectivity matrix
#Row standardize by taking each element and dividing by the row sums
adjmat_rs<-adjmat
for (i in 1:dim(adjmat)[1]){
  adjmat_rs[i,]<-adjmat_rs[i,]/colSums(adjmat, 1)[i]
}
#adjmat_rs is now the row-standardized connectivity matrix (W) for use in the spatial lag model
#Note that values are missing for rows in which there are no neighbors

########form country-level data with neighbor information in 2005
#get the country codes
codes<-row.names(adjmat)
nbr<-data.frame(codes)
#create a variable for whether the country experienced an intrastate armed conflict in 2005
nbr$intra<-ifelse(codes %in% ucdp$gwno[ucdp$year==2005]==TRUE,1,0)
#create a temporal lag variable, which will be used to form the spatial lags and avoid simulteneity bias
nbr$intra_lag<-ifelse(codes %in% ucdp$gwno[ucdp$year==2004]==TRUE,1,0)
#compute the number of conflicts that all the neighbors had in the previous year, by multiplying the non-standarized connectivity matrix times the lagged conflict variable 
nbr$nbr_intra<-adjmat %*% nbr$intra_lag
#Calculate the number of neighbors each state has
ones<-matrix(1,ncol=1,nrow=dim(nbr)[1])
nbr$borders<-adjmat %*% ones
#Put in the year
nbr$year<-matrix(2005,ncol=1,nrow=dim(nbr)[1])

#compute Wy (actually Wy_[t-1])
#Note that this is the same as the number of neighboring conflicts divided by the number neighbors (nbr_intra/borders)
nbr$nbr_intra_rs<-adjmat_rs %*% nbr$intra_lag

#Run the spatial lag model, controlling for number of borders
#Note that because the DV is binary, this is a linear probability model (not ideal, but useful for the example)
spat_lag<-lm(intra~nbr_intra_rs + borders, data=nbr)
summary(spat_lag)
