#Alfalfa germination paper. Analysis initiated 5/18/2015
library(reshape2); library(plyr); library(lattice); library(devtools)
library(ggplot2);library(nlme)
source_gist('https://gist.github.com/jungej62/5a76b72bcd0b12a7ce8a')

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
longdat2$Day<-as.numeric(substr(as.character(longdat2$variable), 2,3))
colnames(longdat2)[7]<-"germ"
longdat3<-droplevels(subset(longdat2, PEG!=0.2))

#Looking at Germ rate by temp within seed ID, note that PEG is held constant
m1dat<-summarySE(subset(longdat3, longdat2$PEG=="0"),
                 measurevar="germ", groupvars=c("SeedID", "Temp", "Day"))
m1dat<-na.omit(m1dat)
ggplot(m1dat, aes(x=Day, y=germ, color=factor(Temp)))+
  facet_wrap(~SeedID)+
  geom_point()+
  geom_line()

#This is what we want to do for each seed id, PEG, and Temp combination
#Here the data are subsetted, but need to find an automated way to do 
#this through all levels.
a2dat<-subset(longdat3, SeedID=="A-2")


#ggplot can do what we need, but it doesn't provide the coefficients.
#We need to find a way to do that outside ggplot
ggplot(a2dat, aes(y=germ, x=Day, color=factor(Rep)))+
       facet_grid(PEG~Temp)+
       geom_point()+
       geom_smooth(method="nls", formula= y ~ (max(y, na.rm=T))*(1-exp((beta-x)/alpha)),
              se=FALSE,
              start = list(alpha = 2, beta = 0.002),
              control = list(maxiter=100000))

#Note that we need to make the function spit out NA if the subsetted data set does
#not have enough values to analyze. That's where the if statement comes in.
#It has to be added as a data frame for the dlply and ldply functions to work
fungerm<-function(data2){
  if(sum(data2$germ, na.rm=T)<0.01){return(as.data.frame(NA))}
  else{
  nls(germ ~ (max(germ, na.rm=T))*(1-exp((beta-Day)/alpha)),
      data=data2, 
      start = list(alpha = 2, beta = 0.002),
      control = list(maxiter=100000))
}
}

#Troubleshooting

fungerm2<-function(data2){
  outdat<-try(nls(germ ~ (max(germ, na.rm=T))*(1-exp((beta-Day)/alpha)),
          data=data2, 
          start = list(alpha = 2, beta = 0.002),
          control = list(maxiter=100000)))
  if(typeof(outdat)=="character"){return(as.data.frame(NA))}
  else{
    nls(germ ~ (max(germ, na.rm=T))*(1-exp((beta-Day)/alpha)),
        data=data2, 
        start = list(alpha = 2, beta = 0.002),
        control = list(maxiter=100000))
  }
}

fungerm2(subset(longdat3, Rep==4&SeedID=="A-2"&PEG==0.1&Temp==1.1))

#This runs the nls model on all SeedID/PEG/Temp/Rep subsettind combinations.
fitList<-dlply(longdat3, .(SeedID, PEG, Temp, Rep), fungerm2)

#This extracts the coefficients from the modesls from the list and puts them in a table.
fitOut<-ldply(fitList, function(x){if(is.data.frame(x)==T){return(cbind(NA,NA))} else{coef(x)}})
colnames(fitOut)[5]<-"alpha";colnames(fitOut)[6]<-"beta"
#Make a new veriable, D, based on equation in paper. This is the day where
#50% of seeds should have germinated
fitOut$D<-fitOut$beta+(fitOut$alpha*log(2))

#Now regress Temp on the inverse of D for each SeedID and PEG combo
germList<-dlply(fitOut, .(SeedID,PEG), function(x) lm((1/D)~Temp, data=x))
germOut<-ldply(germList, function(x){coef(x)})
colnames(germOut)[3]<-"Y_intercept"
#Calculating x intercept, which is base temp for germination
#germOut$X_intercept<-(0-germOut$Y_intercept)/germOut$Temp
#should predict it from the model so we can also get a variance
germOut$X_intercept<-ldply(germList, function(x){predict(x, data.frame("Temp"=0), se.fit = TRUE)$fit})[,3]
germOut$X_int_SE<-ldply(germList, function(x){predict(x, data.frame("Temp"=0), se.fit = TRUE)$se.fit})[,3]


ggplot(germOut,
       aes(y=(X_intercept), x=PEG))+
  facet_wrap(~SeedID)+
  geom_point()+
  geom_pointrange(aes(ymax = X_intercept+X_int_SE, ymin=X_intercept-X_int_SE))

#conducting the anova on x-intercept and slope
#This doesn't work because there are no reps. Need to calculate the regression for
#each rep before trying to do the anova. 
byrep_germList<-dlply(fitOut, .(SeedID,PEG,Rep), function(x) lm((1/D)~Temp, data=x))
byrep_germOut<-ldply(byrep_germList, function(x){coef(x)})
colnames(byrep_germOut)[4]<-"Y_intercept"
byrep_germOut$X_intercept<-(0-byrep_germOut$Y_intercept)/byrep_germOut$Temp

baseanova<-lm(X_intercept~SeedID*factor(PEG), data=byrep_germOut)
slopeanova<-lm(Temp~SeedID*factor(PEG), data=byrep_germOut)

anova(baseanova)
anova(slopeanova)



