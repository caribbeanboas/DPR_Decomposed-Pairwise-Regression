#Package: DPR 1.0
#Type: Package
#Title: Implement the decomposed pairwise regression approach from Koizumi et al. (2006) using an Fst matrix
#Version: 1.0
#Date: 2011-03-03
#Author: Benjamin M. Fitzpatrick and R. Graham Reynolds
#Maintainer: R. Graham Reynolds <rgraham@utk.edu>
#Description: This package contains two funtions. The first DPR1, is for creating the decomposed pairwise regressions for the dataset and for ranking the slope and #intercept of each deme relative to other demes in the dataset to identify non-equilibrium demes ("putative outliers," Koizumi et al. 2006).
#The secong function, DPR2, is to implement the ad hoc AIC analysis to identify non-equilibrium demes in a dataset ("true outliers" Koizumi et al. 2006)
#License: What license is it under?
#LazyLoad: yes

#This program works for 12 demes. If you have more or less, then you have to modify the function “make.dm” in DPR1. For instance, if you had 10 demes, then delete #the last two lines of “make.dm:”
#c(x[56:66],rep(NA,1)),
#c(x[67:78],rep(NA,0))
#)-  be sure to leave this last parenthesis

#The first thing you must do is create your square distance and Fst matrices, named GM.txt and Fst.txt in the directory. Then load the function DRP1, and execute by #calling DPR1(12,1,12,1)

#Then load DPR2 and execute by calling DPR2(1,12) to run the Koizumi et al. 2006 analysis

#For more help see the annotated pdf files or contact me


______________
#\name{DPR 1.0-package}
#\alias{DPR 1.0-package}
#\alias{DPR 1.0}
#\docType{package}
#\title{
#What the package does (short line)
#Implement the decomposed pairwise regression approach from Koizumi et al. (2006) using an Fst matrix}
#\description{
#More about what it does (maybe more than one line)
#This package contains two funtions. The first DPR1, is for creating the decomposed pairwise regressions for the dataset and for ranking the slope and intercept of #each deme relative to other demes in the dataset to identify non-equilibrium demes ("putative outliers," Koizumi et al. 2006).
#The second function, DPR2, is to implement the ad hoc AIC analysis to identify non-equilibrium demes in a dataset ("true outliers" Koizumi et al. 2006)
#}
#\details{
#\tabular{ll}{
#Package: \tab DPR 1.0\cr
#Type: \tab Package\cr
#Version: \tab 1.0\cr
#Date: \tab 2011-03-03\cr
#License: \tab What license is it under?\cr
#LazyLoad: \tab yes\cr
#}
#Step1) establish a dataset as an Fst (or distance) symmetrical matrix and label as "Fdist.txt", save as a tab delimited text file.
#It is important to note that the default works for 12-deme datasets.
#If you have fewer than 12 demes in your dataset, simply add "NA" to the other rows and columns to make a 12x12 square matrix. Then use n=actual pops, y=12.
#If you have more than 12 demes in your dataset, then you will need to modify the function "make.dm", see this function file in the "man" folder.
#Step2) Create a geographic distance matrix as label as "GM", save as a tab delimited text file.
#Make sure that this is a 12x12 matrix, see above
#Step3) Function #1 is "DPR1", and is called by: DPR1(n,x,y,z); where n is the actual number of populations in your dataset, x is the number of independant datasets #(i.e. simulations), y is the number of demes per dataset (should be 12, see above), and z is the number of the deme suspected to be an outlier, which will show up #in the IBD plot
#Step4) Function #2 is "DPR2", and is called by: DPR2(GM,Fdist,x,y); where the variables are the same as in DPR1
#}
#\author{
#Benjamin M. Fitzpatrick and R. Graham Reynolds
#
#Maintainer: Who to complain to <yourfault@somewhere.net>
#R. Graham Reynolds (rgraham@utk.edu)
#}
#\references{
#Reynolds, R.G. 2011. Islands, metapopulations, and archipelagos: genetic equilibrium and non-equilibrium in the context of conservation. Ph.D. dissertation. #University of Tennessee, Knoxville.
#Koizumi, I., Yamamoto, S., and Maekawa, K. 2006. Decomposed pairwise regression analysis of genetic and geographic distances reveals a metapopulation structure of #stream-dwelling Dolly Varden charr. Molecular Ecology 15: 3175-3189.
#}
#Decomposed Pairwise Regression
#DPR 1.0/man
#\keyword{ package }
#\seealso{
#}
#\examples{
#GM=matrix(c(0:6,5:1,1),ncol=12,nrow=12)
#write.table(GM,sep="\t","GM")
#
#fdist=matrix(rbeta(144,.2,.5),ncol=12,nrow=12)
#write.table(fdist,sep="\t","fdist")
#
#DPR1(12,1,12,1)
#DPR2(1,12)
#
#}


___________________________
DPR1 <-
function(n,x,y,z){
#Begin by inputting the geographic distance matrix
#Example Geographic distance matrix
#Warning message OK, refers to brute force wrap-around for this example
#Example is a stepping-stone matrix for 12 demes arranged in a circle with ordination by 1 at each step (i.e. #0,1,2,3,4,5,6,5,4,3,2,1)
GM<-read.table("GM.txt")
GM=as.matrix(as.dist(GM))

#View Geographic Distance Matrix
cat("View Geographic Distance Matrix","\n")
print(GM)
cat("","\n")
#Input is a symmetrical FST Distance Matrix between demes across all simulations
#Example data is for 200 simulations of a 12 deme metapopulation
Fdist<-read.table("fdist.txt")
#Create a data file as a list of FST values for each simulation
fst.mat=as.matrix(Fdist)
#Input the number of rows to correspond to the number of simulations
fst=matrix(t(fst.mat),nrow=x,byrow=TRUE)
#Check that matrix is properly formatted
write.table(fst, "fst.mat.txt")

#Input FST Matrix and convert to a full square matrix
cat("This works for exactly twelve populations","\n")
cat("If your dataset contains fewer, fill in the extra rows and columns of both the Fst and the Geographic distance matrices with NA to create a full twelve by twelve square matrix","\n")
cat("Alternatively, or if your dataset contains  more than twelve populations, check the documentation to refit the function make.dm","\n")
Data<-read.table("fst.mat.txt")
make.dm=function(x,ncol=y){
dm=rbind(
c(x[1],rep(NA,11)),
c(x[2:3],rep(NA,10)),
c(x[4:6],rep(NA,9)),
c(x[7:10],rep(NA,8)),
c(x[11:15],rep(NA,7)),
c(x[16:21],rep(NA,6)),
c(x[22:28],rep(NA,5)),
c(x[29:36],rep(NA,4)),
c(x[37:45],rep(NA,3)),
c(x[46:55],rep(NA,2)),
c(x[56:66],rep(NA,1)),
c(x[67:78],rep(NA,0)))
}
make.dm(Data[1,])

#Extract Simulation #1 for the decomposed IBD Regression
cat("Extract Simulation #1 for the decomposed IBD Regression","\n")
cat("","\n") 
sim1=as.matrix(as.dist(make.dm(Data[1,])))

#Transform FST Matrix to RoussetÂ’s (1997) Distance Measure FST /(1-FST)
sim1.fst<-sim1/(1-sim1)

#Isolation by Distance Regressions
#as.dist takes matrix and turns it into a distance object
#res1 is the linear model (lm) for the regression
res1=lm(as.dist(sim1.fst) ~ as.dist(GM))
cat("","\n")
cat("res1 is the linear model (lm) for the regression","\n")
print(summary(res1))

# Plot the IBD Regression of simulation #1
#Red regression line is the mean IBD regression for all 12 demes in the simulation
plot(GM,sim1.fst, ylab="Genetic Distance",xlab="Geographic Distance", main="Decomposed 
Pairwise Regression for Simulation 1")

#Decomposed IBD regressions for each deme in simulation #1
# Dotted line is the regression of the focal population (Deme 0)
abline(lm(sim1.fst[1,]~GM[1,]), lwd=2, lty=3)
for(i in 2:n){
abline(lm(sim1.fst[i,]~GM[i,]))
}
title(sub="Red regression line is the mean IBD regression for all 12 demes in the simulation
Dotted line is the regression of the focal population (Deme 0)",cex.sub=0.5)

# Now compute the decomposed IBD slopes and intercepts
#Recompose into means
#Draw a table of the following for each simulation: 1)mean intercept 2)mean slope 3)variance of the #intercept 4)variance of the slope.
#Input # of simulations and the # of populations in each simulation
sims=x
pops=y

#Create an empty matrix with nrow = the number of demes and ncol = the number of simulations
Slope=Intercept=matrix(NA,nrow=y,ncol=x)
# Draw an empty table of the following for each simulation: 1)mean intercept 2)mean slope 3)variance of the #intercept 4)variance of the slope.
Means=cbind(mean.I=rep(NA,sims),var.I=rep(NA,sims),mean.S=rep(NA,sims),var.S=rep(NA,sims))

#Loop to fill Slope and Intercept matrices as well as Means table
for(i in 1:sims){
sim=as.matrix(as.dist(make.dm(Data[i,])))
sim=sim/(1-sim)
coefs=lm(as.dist(sim)~as.dist(GM))$coef
for(j in 1:pops){
coefs=rbind(coefs,lm(sim[j,]~GM[j,])$coef)
}
Means[i,1]=mean(coefs[2:y+1,1])
Means[i,2]=var(coefs[2:y+1,1])
Means[i,3]=mean(coefs[2:y+1,2])
Means[i,4]=var(coefs[2:y+1,2])
Slope[,i]=coefs[2:(y+1),2]
Intercept[,i]=coefs[2:(y+1),1]
}

cat("Barplots are only useful for multiple datasets (simulations)","\n")

#Rank the Intercept for each simulation and place in a matrix with nrow = number of demes per simulation 
#and ncol = number of simulations
rank.I=matrix(NA,nrow=y,ncol=x)
vectorranks.I = c()
for (i in 1:sims){
rank.I[,i]=order(Intercept[,i])
}
#Rank the Intercept of the focal deme (deme 1)
vectorranks.I.2=rank.I[z,]

#View graph of Intercept Ranks
windows()
barplot(table(vectorranks.I.2), ylim=c(0,x), col="white", xlab="Rank of Focal Deme Relative to Other Demes in Each Simulation", ylab="Frequency", 
main="Intercept Rank of Focal Deme", sub="1 is highest, 12 is lowest",axis.lty=1)


#Rank the Slope for each simulation and place in a matrix with nrow = number of demes per simulation 
#and ncol = number of simulations
rank.S=matrix(NA,nrow=y,ncol=x)

vectorranks.S = c()
for (i in 1:sims){
rank.S[,i]=order(Slope[,i])
}
#Rank the Slope of the focal deme (deme 1)
vectorranks.S.2=rank.S[z,]

#View graph of Slope Ranks
windows()
barplot(table(vectorranks.S.2), ylim=c(0,x), col="white", xlab="Rank of Focal Deme Relative to Other Demes in Each Simulation", ylab="Frequency", 
main="Slope Rank of Focal Deme",sub="1 is highest, 12 is lowest", axis.lty=1)

#View tables of Slope and Intercept Ranks
#Top row is the Rank (from 1 to the maximum # of demes in simulations; example is 1-12) where 1 is the highest #and 12 is the lowest
#Bottom row is the number of focal demes that fall under each rank class across all simulations
#Recall that each focal deme is ranked relative to others in its simulation (metapopulation)
cat("Table of Interept Ranks","\n")
cat("Top row is the Rank (from 1 to the maximum # of demes in simulations; example is 1-12) where 1 is the highest and 12 is the lowest","\n")
cat("Bottom row is the number of focal demesthat fall under each rank class across all simulations","\n")
cat("Recall that each focal deme is ranked relative to others in its simulation (metapopulation)","\n")
print(table(vectorranks.I.2))
cat("","\n")
cat("Table of Slope Ranks")
cat("Top row is the Rank (from 1 to the maximum # of demes in simulations; example is 1-12) where 1 is the highest and 12 is the lowest","\n")
cat("Bottom row is the number of focal demes that fall under each rank class across all simulations","\n")
cat("Recall that each focal deme is ranked relative to others in its simulation (metapopulation)","\n")
print(table(vectorranks.S.2))
cat("","\n")

#Write tables to text files
cat("","\n")
cat("Two text files of each table are written to the directory","\n")
write.table(vectorranks.S.2,sep="\t","Slope.txt",col.names=NA)
write.table(vectorranks.I.2,sep="\t","Intercept.txt",col.names=NA)
cat("","\n")
cat("Three graphs should be generated, one regression plot and two barplots","\n")
cat("","\n")
}



______________________



DPR2 <-
function(x,y){
#Begin by inputting the geographic distance matrix
#Example Geographic distance matrix
#Warning message OK, refers to brute force wrap-around for this example
#Example is a stepping-stone matrix for 12 demes arranged in a circle with ordination by 1 at each step (i.e. #0,1,2,3,4,5,6,5,4,3,2,1)
GM<-read.table("GM.txt")
GM=as.matrix(as.dist(GM))

#View Geographic Distance Matrix
cat("View Geographic Distance Matrix","\n")
print(GM)
cat("","\n")

#Input is a symmetrical FST Distance Matrix between demes across all simulations
#Data example is for 200 simulations of a 12 deme metapopulation
Data<-read.table("fdist.txt")

#Create a data file as a list of FST values for each simulation
fst.mat=as.matrix(Data)

#Input the number of rows to correspond to the number of simulations
fst=matrix(t(fst.mat),nrow=x,byrow=TRUE)

#Check that matrix is properly formatted
write.table(fst, "fst.mat.txt")

#Input the number of simulations in the FST matrix and the number of demes in each simulation
sims=x
pops=y

#Input the number of rows and columns in the FST matrix file ("fst.mat")
#A file of 12 demes will have nrow=13
sim.k=sim.l=matrix(NA,nrow=y+1,ncol=x)

#best.k is to identify the demes from each simulation that are outliers in the Koizumi et al. (2006) analysis
#best.c is to identify the demes from each simulation that are outliers in the AICC model selection from #Burnham #and Anderson (2004)
best.k=best.c=NA

#Label demes
deme=0:y

#Vr is an empty matrix with ncol= # of demes plus 1, and nrow= # of simulations
#Vr will be filled by the following loop
#LL.K is the Log Likelihood from the linear models
#AIC.K is the AIC for the Koizumi analysis, defined as AIC = 2K + n ln(RSS/n)
#where K is the # of parameters, n is the # of populations, and RSS is the residual sum of squares
#AIC.C is the AIC analysis with a correction for small sample size from Burnham and Anderson (2004). 
Vr=LL.K=AIC.K=AIC.C=matrix(NA,ncol=pops+1,nrow=sims)

#Loop analysis for all simulations
for (i in 1:sims){
sim=matrix(fst[i,],ncol=y)
#Transform FST Matrix to RoussetÂ’s (1997) Distance Measure FST /(1-FST)
sim.fst<-sim/(1-sim)
Fit.0=lm(as.dist(sim.fst)~as.dist(GM))
vr=var(Fit.0$res)
ll.koiz=logLik(Fit.0)
aic.koiz=2*2+pops*log(sum(Fit.0$res^2)/pops)+2*2*(2+1)/(pops-2-1)
aic.c=-2*logLik(Fit.0)*2/(pops-1)+4+12/(pops-3)
for (j in 1:pops){
fit.d=lm(as.dist(sim.fst[-j,-j])~as.dist(GM[-j,-j]))
vr=c(vr,var(fit.d$res))
ll.koiz=c(ll.koiz,logLik(fit.d))
aic.koiz=c(aic.koiz,4+(pops-1)*log(sum(fit.d$res^2)/(pops-1))+12/((pops-1)-3))
aic.c=c(aic.c,-2*logLik(fit.d)*2/(pops-2)+12/(pops-4))
}
sim.k[,i]=order(aic.koiz)
best.k[i]=deme[order(aic.koiz)][1]
best.c[i]=deme[order(aic.c)][1]
LL.K[i,]=ll.koiz
AIC.K[i,]=aic.koiz
AIC.C[i,]=aic.c
Vr[i,]=vr
}

cat("Print AIC.K scores","\n")
print(AIC.K)
cat("","\n")
#Best.k is a table showing the deme number on the top row and the number of simulations where that deme was #identified as an outlier in the bottom row
#In the example, there are 12 demes, identified as 1-12, and 200 simulations
#Demes which have undergone an historical event should be identified as outliers
#For example, if deme #1 underwent a dramatic bottleneck in 200 independent simulations, this table will show how #many times deme 1 was identified as an outlier relative to other demes in the simulation
cat("Best.k is a table showing the deme number on the top row and the number of simulations where that deme was identified as an outlier in the bottom row","\n")
cat("Demes which have undergone an historical event should be identified as outliers","\n")
cat("For example, if deme #1 underwent a dramatic bottleneck in 200 independent simulations, this table will show how many times deme 1 was identified as an outlier relative to other demes in the simulation","\n")
print(table(best.k))
cat("","\n")

#Write the table of outliers as a text file
cat("Table of outliers written as a text file (best.k)","\n")
write.table (best.k, "best.k.txt")
cat("","\n")
}

