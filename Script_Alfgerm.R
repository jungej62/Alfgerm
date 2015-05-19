#Alfalfa germination paper. Analysis initiated 5/18/2015
library(reshape2); library(plyr); library(lattice)
dat<-read.csv("ConstantTempsAll.csv", header=T)
str(dat)
dat$fRep<-factor(dat$Rep)
dat$fPEG<-factor(dat$PEG)
dat$fTemp<-factor(dat$Temp)
#make a new data frame that includes that original counts of germinated
#seeds observed each day, as well as a new column for each day 
#that includes the proportion germinated. This will be cumulative,
#so for day 5, P5 will be the total number of germinated seeds up to that
#day divided by the number of seeds in the initial test
dat2<-cbind(dat, matrix(data=NA, nrow=1386, ncol=31, 
                        dimnames=list(dimnames(dat)[[1]],c(paste(c(rep("P",31)),c(1:31), sep="")))))
#calculate the proportion for the first day by hand
dat2[40]<-dat2[6]/dat2[5]
#Loop to calculate proportions from day 7-end
for(i in 7:36){
  tsamp<-apply(dat2[,6:i],1, sum, na.rm=T)
  dat2[i+34]<-tsamp/dat2[5]
  
}
#Add NA to proportion columns where counts were NA
for(i in 1:length(dat2[,1])){
  for(j in 6:36){
    if(is.na(dat2[i,j])==TRUE){
      dat2[i,j+34]<-NA
    }
  }
}

dat2<-cbind(dat2, apply(dat2[,6:36],1, sum, na.rm=T)/dat2[5])
colnames(dat2)[71]<-"FinalGerm"
longdat2<-melt(dat2[,c(1:5,40:70)], id=c("Rep","SeedID", "PEG", "Temp","StartSeeds"))
longdat2<-arrange(longdat2, SeedID, PEG, Temp, Rep)
longdat2$Day<-substr(as.character(longdat2$variable), 2,3)



