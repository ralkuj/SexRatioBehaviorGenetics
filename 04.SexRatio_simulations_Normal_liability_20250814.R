################################################################################
###
### File: 04.SexRatio_simulations_Normal_liability_20250814.R
### Author: Ralf Kuja-Halkola
### Purpose: Simulate different levels of power for detecting non-null heritability
###          for sex ratio.
### Date: 2024-11-07
### Updated: 2025-02-25, 2025-08-14 (Changing simulation approach), 2026-05-07 (commenting)
###
################################################################################


################################################################################
### The references:
# Song&Zhang:
#   Song S, Zhang J. In search of the genetic variants of human sex ratio at birth: was Fisher wrong about sex ratio evolution? Proc Biol Sci. Oct 2024;291(2033):20241876. doi:10.1098/rspb.2024.1876
# Zietsch et al:
#   Zietsch BP, Walum H, Lichtenstein P, Verweij KJH, Kuja-Halkola R. When theory cannot explain data, the theory needs rethinking. Invited replies to: Orzack SH, Hardy ICW. 2021, and Lehtonen J. 2021. Proc Biol Sci. Mar 31 2021;288(1947):20210304. doi:10.1098/rspb.2021.0304
################################################################################


################################################################################
### Libraries
require(MASS)
library(drgee)
library(parallel)
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
### Define function

SexRatioSimulation <- function( X=1 , h2 = 0.43 ){
### 1. Give the unique grand parents a random normal variable, separated into additive genetics and residual
  datMOR$VmorA <- rnorm(length(unique(datInd$LOPNRMOR))) # Mother
  datFAR$VfarA <- rnorm(length(unique(datInd$LOPNRFAR))) # Father
### 2. Merge to index generation's data
  datInd <- merge( x=datInd , y=datMOR , by='LOPNRMOR' , all=T )
  datInd <- merge( x=datInd , y=datFAR , by='LOPNRFAR' , all=T )
### 3. Calculate individuals sex ratio as random with expected mean of parents for additive genetics
  datInd$VindA <- (datInd$VmorA+datInd$VfarA)/2 + rnorm( n=dim(datInd)[1] )*sqrt(.5)
  # But overall correlation 0.5*h2
  datInd$Vind <- datInd$VindA*sqrt(h2) + rnorm(dim(datInd)[1])*sqrt(1-h2)
  # Note that this setup will ensure .5 correlation between individuals and their parents and siblings additive genetics, while the total correlation will be 0.5*h2

### 4. Merge the individual to the original data
  datChild <- merge( x=datChild , y=datInd[,c('LOPNR','Vind')] , by='LOPNR' , all.x=T )

### 5. Create simulation probability according to Song&Zhang and simulate sex of offspring
  datChild$sexsim <- (datChild$Vind-mean(datChild$Vind))/sd(datChild$Vind)
  datChild$sexsim <- datChild$sexsim*.025+0.5
  datChild$sexSZ <- rbinom(n=dim(datChild)[1],1,datChild$sexsim)
### 6. Use the same simulation, but change the standard deviation compared to assumed by Song&Zhang
  datChild$sexSZ_080sd <- rbinom(n=dim(datChild)[1],1,(datChild$sexsim-.5)*0.8+.5)
  datChild$sexSZ_120sd <- rbinom(n=dim(datChild)[1],1,(datChild$sexsim-.5)*1.2+.5)
  datChild$sexSZ_150sd <- rbinom(n=dim(datChild)[1],1,(datChild$sexsim-.5)*1.5+.5)
  datChild$sexSZ_200sd <- rbinom(n=dim(datChild)[1],1,(datChild$sexsim-.5)*2  +.5)
  
### 7. Re-merge with original data
  datNew <- merge( x=datNew , y=datChild[,c('LOPNRBARN','sexsim','sexSZ','sexSZ_080sd','sexSZ_120sd','sexSZ_150sd','sexSZ_200sd') ] , by='LOPNRBARN' , all.x=T )
  colnames(datChild) <- paste0(colnames(datChild),'2')
  datNew <- merge( x=datNew , y=datChild[,c('LOPNRBARN2','sexsim2','sexSZ2','sexSZ_080sd2','sexSZ_120sd2','sexSZ_150sd2','sexSZ_200sd2') ] , by='LOPNRBARN2' , all.x=T )
### 8. Analyses
  # Simulated: Song&Zhang
  fit_SZ <- summary(gee( sexSZ ~ sexSZ2 , data=datNew , link = 'logit' , clusterid = 'LOPNRMOR' ) )$coef[2,]
  # Simulated: Song&Zhang other standard errors
  fit_SZ080sd <- summary(gee( sexSZ_080sd ~ sexSZ_080sd2 , data=datNew , link = 'logit' , clusterid = 'LOPNRMOR' ) )$coef[2,]
  fit_SZ120sd <- summary(gee( sexSZ_120sd ~ sexSZ_120sd2 , data=datNew , link = 'logit' , clusterid = 'LOPNRMOR' ) )$coef[2,]
  fit_SZ150sd <- summary(gee( sexSZ_150sd ~ sexSZ_150sd2 , data=datNew , link = 'logit' , clusterid = 'LOPNRMOR' ) )$coef[2,]
  fit_SZ200sd <- summary(gee( sexSZ_200sd ~ sexSZ_200sd2 , data=datNew , link = 'logit' , clusterid = 'LOPNRMOR' ) )$coef[2,]
  return(data.frame( rbind( fit_SZ , fit_SZ080sd , fit_SZ120sd , fit_SZ150sd , fit_SZ200sd ) ) )
}

# Test
# tt0 <- Sys.time()
# testsim <- mclapply(X=rep(1,4) , SexRatioSimulation , mc.cores=2 , h2=0.43 )
# tt1 <- Sys.time()
# tt1-tt0

### Run parallel process for simulation
h2run  <- c(.2,.4,.6,.8,1)
h2name <- paste0('_',round(h2run*100))
tt0 <- Sys.time()
set.seed(12341234)
for(i in 1:length(h2run)){
  mSimulation <- mclapply(X=rep(1,1000) , SexRatioSimulation , mc.cores=50 , h2=h2run[i] )
#  length(mSimulation)
  # Power
  pvalsM <- data.frame( sexSZ=1:1000,sexSZ_080sd=NA,sexSZ_120sd=NA,sexSZ_150sd=NA,sexSZ_200sd=NA )
#  str(pvalsM)
  for(j in 1:dim(pvalsM)[1]){ pvalsM[j,] <- mSimulation[[j]][,4] }
#  str(pvalsM)
#  summary(pvalsM)
  print(paste0('Heritability',h2name[i],'%'))
  print(apply( pvalsM , 2 , function(x){ mean(x<0.05,na.rm=T) } ))
  print(Sys.time()-tt0)
  # Save
  saveRDS( pvalsM ,  paste0('~/Sex ratio/Simulation_Realdata_Male_1000repeats_h2',h2name[i],'percent_20250814.Rds') )
}
tt1 <- Sys.time()
tt1-tt0

#####################


################################################################################


################################################################################
### Create a table for power
library(tinytable)
simh2_20 <- readRDS('~/Sex ratio/Simulation_Realdata_Male_1000repeats_h2_20percent_20250814.Rds')
simh2_40 <- readRDS('~/Sex ratio/Simulation_Realdata_Male_1000repeats_h2_40percent_20250814.Rds')
simh2_60 <- readRDS('~/Sex ratio/Simulation_Realdata_Male_1000repeats_h2_60percent_20250814.Rds')
simh2_80 <- readRDS('~/Sex ratio/Simulation_Realdata_Male_1000repeats_h2_80percent_20250814.Rds')
simh2_100 <- readRDS('~/Sex ratio/Simulation_Realdata_Male_1000repeats_h2_100percent_20250814.Rds')
pow_20  <- apply( simh2_20 , 2 , function(x){ mean( x<0.05 ) } )[ c(2,1,3,4,5) ]
pow_40  <- apply( simh2_40 , 2 , function(x){ mean( x<0.05 ) } )[ c(2,1,3,4,5) ]
pow_60  <- apply( simh2_60 , 2 , function(x){ mean( x<0.05 ) } )[ c(2,1,3,4,5) ]
pow_80  <- apply( simh2_80 , 2 , function(x){ mean( x<0.05 ) } )[ c(2,1,3,4,5) ]
pow_100 <- apply( simh2_100 , 2 , function(x){ mean( x<0.05 ) } )[ c(2,1,3,4,5) ]

outTab <- data.frame( cbind( pow_20 , pow_40 , pow_60 , pow_80 , pow_100 ) )
rownames(outTab) <- c(paste0('Song&Zhang - ',c('20'),'%'),'Song&Zhang',paste0('Song&Zhang + ',c('20','50','100'),'%'))
colnames(outTab) <- paste0( 'Heritability ',c(20,40,60,80,100),'%')
outTab
save_tt( x=tt(x=outTab , rownames=T) , 'Z:/Project/Finished projects/SexRatio karin Verwei/2024-10-25 Possible letter/Simulation normal variance/SimulationResults_NormalPower_202050814.docx' , overwrite=T )
################################################################################



################################################################################
################################################################################
############################## END OF FILE #####################################
################################################################################
################################################################################