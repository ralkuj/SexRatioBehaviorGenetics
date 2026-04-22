################################################################################
###
### File: 04.SexRatio_simulations_Normal_liability_20250225.R
### Author: Ralf Kuja-Halkola
### Purpose: Simulate different levels of power for detecting non-null heritability
###          for sex ratio.
### Date: 2024-11-07
### Updated: 2025-02-25
###
################################################################################


################################################################################
### The references:
# Song&Zhang:
#   Song S, Zhang J. In search of the genetic variants of human sex ratio at birth: was Fisher wrong about sex ratio evolution? Proc Biol Sci. Oct 2024;291(2033):20241876. doi:10.1098/rspb.2024.1876
# Zietsch et al:
#   Zietsch BP, Walum H, Lichtenstein P, Verweij KJH, Kuja-Halkola R. No genetic contribution to variation in human offspring sex ratio: a total population study of 4.7 million births. Proc Biol Sci. 2020;287(1921):20192849. DOI: 10.1098/rspb.2019.2849
################################################################################


################################################################################
### Libraries used
require(MASS)
library(drgee)
library(parallel)

### Nr cores used
nrcores <- 50

################################################################################


################################################################################
### Data pre-steps
### Use the sibling structure in Zietsch et al, use only men
### In this simulation the spouses of the men are ignored. 

### Load
load('~/Sex ratio/Familiality_of_offspring_sex_cousinpairs_20190116.Rdata')

# Remove twins in offspring generation
dat <- dat[ dat$nrbornsamedate==1 & dat$nrbornsamedate2==1 , ]
# Restrict data to full sibs only
dat <- dat[ dat$sibtype==0 , c('LOPNR','LOPNR2','KON','KON2','LOPNRSPOUSE','LOPNRSPOUSE2','LOPNRMOR','LOPNRFAR','LOPNRBARN','LOPNRBARN2','konbarn','konbarn2' ) ]
# Only men
dat <- dat[ dat$KON==0 & dat$KON2==0 , ]
# Get all index-generation individuals (including spouses) in the data
datInd <- dat[ !duplicated(dat$LOPNR) , c('LOPNR','LOPNRMOR','LOPNRFAR') ]
# Create child dataset
datChild <- dat[ , c('LOPNRBARN','LOPNR') ]
datChild <- datChild[ !duplicated(dat$LOPNRBARN) , ]
# Create grandparent dataset
datMOR <- data.frame( LOPNRMOR=unique(datInd$LOPNRMOR) )
datFAR <- data.frame( LOPNRFAR=unique(datInd$LOPNRFAR) )
# A separate copy to be used in simulation
datNew <- dat
################################################################################


################################################################################
### Define function for simulating offspring sex based on heritable sex allocation
### The analytic approach is the same as in Zietsch et al
SexRatioSimulation <- function( X=1 ){
### 1. Give the unique grand parents a random normal variable
  datMOR$Vmor <- rnorm(length(unique(datInd$LOPNRMOR)))
  datFAR$Vfar <- rnorm(length(unique(datInd$LOPNRFAR)))
### 2. Merge to index generation's data
  datInd <- merge( x=datInd , y=datMOR , by='LOPNRMOR' , all=T )
  datInd <- merge( x=datInd , y=datFAR , by='LOPNRFAR' , all=T )
### 3. Calculate individuals ratio as random with expected mean of parents
  datInd$Vind <- (datInd$Vmor+datInd$Vfar)/2 + sqrt(.5)*rnorm( n=dim(datInd)[1] )
  # Note that this setup will ensure .5 correlation between individuals and their parents and siblings
### 4. Merge the individual to the original data
  datChild <- merge( x=datChild , y=datInd[,c('LOPNR','Vind')] , by='LOPNR' , all.x=T )
### 5. Create simulation probability according to Song&Zhang and simulate sex of offspring
  datChild$sexsim <- (datChild$Vind-mean(datChild$Vind))/sd(datChild$Vind)
  datChild$sexsim <- datChild$sexsim*.025+0.5
  datChild$sexSZ <- rbinom(n=dim(datChild)[1],1,datChild$sexsim)
### 6. Use the same simulation, but change the standard deviation compared to assumed by Song&Zhang
  datChild$sexSZ_080sd <- rbinom(n=dim(datChild)[1],1,(datChild$sexsim-.5)*0.8+.5)
  datChild$sexSZ_090sd <- rbinom(n=dim(datChild)[1],1,(datChild$sexsim-.5)*0.9+.5)
  datChild$sexSZ_110sd <- rbinom(n=dim(datChild)[1],1,(datChild$sexsim-.5)*1.1+.5)
  datChild$sexSZ_120sd <- rbinom(n=dim(datChild)[1],1,(datChild$sexsim-.5)*1.2+.5)
  datChild$sexSZ_130sd <- rbinom(n=dim(datChild)[1],1,(datChild$sexsim-.5)*1.3+.5)
  datChild$sexSZ_150sd <- rbinom(n=dim(datChild)[1],1,(datChild$sexsim-.5)*1.5+.5)
  datChild$sexSZ_175sd <- rbinom(n=dim(datChild)[1],1,(datChild$sexsim-.5)*1.75+.5)
  
### 7. Re-merge with original data
  datNew <- merge( x=datNew , y=datChild[,c('LOPNRBARN','sexsim','sexSZ','sexSZ_080sd','sexSZ_090sd','sexSZ_110sd','sexSZ_120sd','sexSZ_130sd','sexSZ_150sd','sexSZ_175sd') ] , by='LOPNRBARN' , all.x=T )
  colnames(datChild) <- paste0(colnames(datChild),'2')
  datNew <- merge( x=datNew , y=datChild[,c('LOPNRBARN2','sexsim2','sexSZ2','sexSZ_080sd2','sexSZ_090sd2','sexSZ_110sd2','sexSZ_120sd2','sexSZ_130sd2','sexSZ_150sd2','sexSZ_175sd2') ] , by='LOPNRBARN2' , all.x=T )
### 8. Analyses
  # Simulated: Song&Zhang
  fit_SZ <- summary(gee( sexSZ ~ sexSZ2 , data=datNew , link = 'logit' , clusterid = 'LOPNRMOR' ) )$coef[2,]
  # Simulated: Song&Zhang other standard errors
  fit_SZ080sd <- summary(gee( sexSZ_080sd ~ sexSZ_080sd2 , data=datNew , link = 'logit' , clusterid = 'LOPNRMOR' ) )$coef[2,]
  fit_SZ090sd <- summary(gee( sexSZ_090sd ~ sexSZ_090sd2 , data=datNew , link = 'logit' , clusterid = 'LOPNRMOR' ) )$coef[2,]
  fit_SZ110sd <- summary(gee( sexSZ_110sd ~ sexSZ_110sd2 , data=datNew , link = 'logit' , clusterid = 'LOPNRMOR' ) )$coef[2,]
  fit_SZ120sd <- summary(gee( sexSZ_120sd ~ sexSZ_120sd2 , data=datNew , link = 'logit' , clusterid = 'LOPNRMOR' ) )$coef[2,]
  fit_SZ130sd <- summary(gee( sexSZ_130sd ~ sexSZ_130sd2 , data=datNew , link = 'logit' , clusterid = 'LOPNRMOR' ) )$coef[2,]
  fit_SZ150sd <- summary(gee( sexSZ_150sd ~ sexSZ_150sd2 , data=datNew , link = 'logit' , clusterid = 'LOPNRMOR' ) )$coef[2,]
  fit_SZ175sd <- summary(gee( sexSZ_175sd ~ sexSZ_175sd2 , data=datNew , link = 'logit' , clusterid = 'LOPNRMOR' ) )$coef[2,]
  return(data.frame( rbind( fit_SZ , fit_SZ080sd , fit_SZ090sd , fit_SZ110sd , fit_SZ120sd , fit_SZ130sd , fit_SZ150sd , fit_SZ175sd ) ) )
}

# Test
# tt0 <- Sys.time()
# testsim <- mclapply(X=rep(1,nrcores) , SexRatioSimulation , mc.cores=nrcores )
# tt1 <- Sys.time()
# tt1-tt0

### Run parallel process for simulation
nrep <- 1000
tt0 <- Sys.time()
set.seed(12341234)
maleSimulation <- mclapply(X=rep(1,nrep) , SexRatioSimulation , mc.cores=nrcores )
tt1 <- Sys.time()
tt1-tt0
length(maleSimulation)
maleSimulation[[10]]
# Power
pvalsMale <- data.frame( sexSZ=1:nrep,sexSZ_080sd=NA,sexSZ_090sd=NA,sexSZ_110sd=NA,sexSZ_120sd=NA,sexSZ_130sd=NA,sexSZ_150sd=NA,sexSZ_175sd=NA )
for(i in 1:dim(pvalsMale)[1]){ pvalsMale[i,] <- maleSimulation[[i]][,4] }
str(pvalsMale)
summary(pvalsMale)
apply( pvalsMale , 2 , function(x){ mean(x<0.05,na.rm=T) } )
# Save
saveRDS( pvalsMale ,  '~/Sex ratio/Simulation_Realdata_Male_1000repeats_20241107.Rds')

### Final run 20241107:
# sexSZ sexSZ_080sd sexSZ_090sd sexSZ_110sd sexSZ_120sd sexSZ_130sd sexSZ_150sd sexSZ_175sd
# 0.538       0.276       0.422       0.713       0.860       0.938       0.997       1.000

################################################################################


################################################################################
### Create a plot of densities for the selected standard deviations
### Using final run values from above, hard coded so no need to run above 
### simulation to plot

# Fix densities to plot
x <- seq(.3,.7,len=400)
vals <- c(.8,.9,1,1.1,1.2,1.3,1.5,1.75)
y <- matrix( NA , length(x) , length(vals) )
for( i in 1: length(vals) ){ 
  y[,i] <- dnorm( x=x , mean=.5 , sd=.025*vals[i] ) 
  y[,i] <- y[,i]/max(y[,i])
}
colnames(y) <- c(vals*.025)

# Colors and line types
cols <- c( rgb(0.25, 0.25, 0.25),rgb(0.5, 0.5, 0.5) , rgb(0, 0, 0) , rgb(0.5, 0.3, 0.3),rgb(0.65, 0.4, 0.4),rgb(0.8, 0.55, 0.55),rgb(0.9, 0.7, 0.7),rgb(0.95, 0.8, 0.8) )
ltys <- c(3,2,1,2,3,4,5,6)


# Plot
pdf(file='~/Sex ratio/Densities_and_power_20241107.pdf',width=29/2.54,height=14/2.54)

par( oma=c(1,1,1,1) )
par( mar=c(2,3,1,1) )

plot(NULL , xlim=c(.3,.7) , ylim=c(0,1) , axes=F , xlab='Sex ratio' , ylab=NA , xaxs='i' , yaxs='i' )
for(i in 1:dim(y)[2]){
  lines( x , y[,i] , col=cols[i] , lty=ltys[i] , lwd=2 )
}
#
par( xpd=NA )
text( x=.28,y=1.015 , 'Standard deviation in relation to Song&Zhang' , adj=0 )
legend( x=.28 , y=1 , legend=c(paste0('Song&Zhang - ',c('20','10'),'%'),'Song&Zhang',paste0('Song&Zhang + ',c('10','20','30','50','75'),'%')) , col=cols , lty=ltys , lwd=2 , bty='n' )
text( x=.275,y=.51 , 'Coefficient of variation' , adj=0 , cex=.8 )
legend( x=.28 , y=.5 , legend=vals*.025/.5 , col=cols , lty=ltys , lwd=2 , bty='n' , cex=.8 )
text( x=.34,y=.51 , 'Standard deviation' , adj=0 , cex=.8 )
legend( x=.34 , y=.5 , legend=vals*.025 , col=cols , lty=ltys , lwd=2 , bty='n' , cex=.8 )
text( x=.62,y=.842 , 'Power' , adj=0 , cex=1.6 )
legend( x=.6 , y=.83 , legend=paste0(c('28','42','54','71','86','94','>99','>99'),'%') , col=cols , lty=ltys , lwd=2 , bty='n' , cex=1.4 )
par( xpd=F )

axis( side=1 , at=seq(.3,.7,by=.1) )

dev.off()

################################################################################


################################################################################
################################################################################
############################## END OF FILE #####################################
################################################################################
################################################################################
