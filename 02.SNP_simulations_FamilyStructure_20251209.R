################################################################################
###
### File: 02.SNP_simulations_FamilyStructure_20251209.R
### Author: Ralf Kuja-Halkola
### Purpose: Implement the power simulation for Zietsch et all to detect familiaö
###          aggregation (heritability) based on one SNP with differing MAF and 
###          effect size.
### Date: 20251209
### Updated: 20260422, commenting and cleaning up code.
###
################################################################################


############################################
### The references:
# Song&Zhang:
#   Song S, Zhang J. In search of the genetic variants of human sex ratio at birth: was Fisher wrong about sex ratio evolution? Proc Biol Sci. Oct 2024;291(2033):20241876. doi:10.1098/rspb.2024.1876
# Zietsch et al:
#   Zietsch BP, Walum H, Lichtenstein P, Verweij KJH, Kuja-Halkola R. No genetic contribution to variation in human offspring sex ratio: a total population study of 4.7 million births. Proc Biol Sci. 2020;287(1921):20192849. DOI: 10.1098/rspb.2019.2849
############################################


################################################################################
### Libraries used
library(drgee)
library(parallel)
################################################################################


################################################################################
### Data pre-steps
### Use the sibling structure in Zietsch et al, use only men
### In this simulation the spouses of the men are ignored. 

### Load data
load('~/Sex ratio/Familiality_of_offspring_sex_cousinpairs_20190116.Rdata')

# Remove twins in offspring generation
dat <- dat[ dat$nrbornsamedate==1 & dat$nrbornsamedate2==1 , ]
# Restrict data to full sibs only
dat <- dat[ dat$sibtype==0 , c('LOPNR','LOPNR2','KON','KON2','LOPNRSPOUSE','LOPNRSPOUSE2','LOPNRMOR','LOPNRFAR','LOPNRBARN','LOPNRBARN2','konbarn','konbarn2','bnr','bnr2' ) ]
# Only men
dat <- dat[ dat$KON==0 & dat$KON2==0 , ]
# Get all index-generation individuals (including spouses) in the data
datInd <- dat[ !duplicated(dat$LOPNR) , c('LOPNR','LOPNRMOR','LOPNRFAR') ]
# Create child dataset
datChild <- dat[ , c('LOPNRBARN','LOPNR') ]
datChild <- datChild[ !duplicated(dat$LOPNRBARN) , ]
datChildOne <- datChild[ !duplicated(datChild$LOPNR) , ]
# Create grandparent dataset
datMOR <- data.frame( LOPNRMOR=unique(datInd$LOPNRMOR) )
datFAR <- data.frame( LOPNRFAR=unique(datInd$LOPNRFAR) )
# A separate copy to be used in simulation
datNew <- dat

################################################################################


################################################################################
### Simulation Function

famfun <- function( X=1 , MAF , bet , SR=0.514 , first=FALSE ){
### 1. Simulate alleles at the different chromosomes in the parents
  Mother <- data.frame( LOPNRMOR=datMOR$LOPNRMOR , m1=rbinom(length(datMOR$LOPNRMOR),1,MAF) ,  m2=rbinom(length(datMOR$LOPNRMOR),1,MAF) )
  Father <- data.frame( LOPNRFAR=datFAR$LOPNRFAR , p1=rbinom(length(datFAR$LOPNRFAR),1,MAF) ,  p2=rbinom(length(datFAR$LOPNRFAR),1,MAF) )
### 2. Merge to index generation's data
  Individual <- merge( x=datInd , y=Mother , by='LOPNRMOR' , all=T )
  Individual <- merge( x=Individual , y=Father , by='LOPNRFAR' , all=T )
### 3. Draw alleles for individuals
  Individual$a1 <- apply( X=cbind(Individual$m1,Individual$m2) , MARGIN=1 , FUN=function(x){sample(x,1)} )
  Individual$a2 <- apply( X=cbind(Individual$p1,Individual$p2) , MARGIN=1 , FUN=function(x){sample(x,1)} )
  Individual$genotype <- Individual$a1+Individual$a2
  if( !first ){
### 4v1. Merge child
    Children <- merge( x=datChild , y=Individual[,c('LOPNR','genotype')] , by='LOPNR' , all.x=T )
  }
  if( first ){
### 4v2. Merge child, only one child per parent
    Children <- merge( x=datChildOne , y=Individual[,c('LOPNR','genotype')] , by='LOPNR' , all.x=T )
  }
### 5. Simulation of offspring sex
  Children$sexsim <- rbinom(dim(Children)[1],1,SR+bet*Children$genotype-2*bet*MAF)
  NewData <- merge( x=datNew , y=Children[,c('LOPNRBARN','sexsim') ] , by='LOPNRBARN' , all.x=T )
  colnames(Children) <- paste0(colnames(Children),'2')
  NewData <- merge( x=NewData , y=Children[,c('LOPNRBARN2','sexsim2') ] , by='LOPNRBARN2' , all.x=T )
  
### 6. Return p-value of analysis
  return(summary(gee(sexsim~sexsim2,data=NewData,link='logit',clusterid='LOPNRMOR'))$coef[2,4])
}

################################################################################


################################################################################
### Simulate for fitst-born offspring only
# Number of repeats for each simulation
n <- 1000
# Set of parameters to investigate
MAFseq <- c(.0025,.005,.01,.025,.05,.1)
betseq <- seq(.01,.21,by=.02)

# Total number of repeats
n*length(MAFseq)*length(betseq)

### Run - One offspring
simul <- matrix(NA,length(MAFseq),length(betseq))
tt0 <- Sys.time()
for(i in 1:length(MAFseq) ){
  for(j in 1:length(betseq) ){
    simul[i,j] <- mean( unlist( mclapply( X=rep(1,n) , FUN=famfun , MAF=MAFseq[i] , bet=betseq[j] , SR=0.514 , first=T , mc.cores=50 ) ) < 0.05 )
  }
  print( paste('MAF',MAFseq[i] ,'done') )
}
tt1 <- Sys.time()
tt1-tt0
rownames(simul) <- paste( 'MAF =', MAFseq )
colnames(simul) <- paste( 'bet =', betseq )
simul

saveRDS(object=simul , file='~/Sex ratio/SimulationSNP_FamilyStructure_male_First_1000repeats_20251209.Rds' )

### Adding more points since too crude steps and plot looks bad
# Number of repeats for each simulation
n <- 1000
# Set of parameters to investigate
MAFseq <- c(.075)
betseq <- seq(.01,.21,by=.02)

# Total number of repeats
n*length(MAFseq)*length(betseq)
#n*length(MAFseq)*length(betseq)*8.35/10 /50 /60

### Run - One offspring
simul <- matrix(NA,length(MAFseq),length(betseq))
tt0 <- Sys.time()
for(i in 1:length(MAFseq) ){
  for(j in 1:length(betseq) ){
    simul[i,j] <- mean( unlist( mclapply( X=rep(1,n) , FUN=famfun , MAF=MAFseq[i] , bet=betseq[j] , SR=0.514 , first=T , mc.cores=50 ) ) < 0.05 )
    print( paste('MAF',MAFseq[i] ,', beta',betseq[j],'done.') )
  }
  print( paste('MAF',MAFseq[i] ,'done.') )
  print( paste('Time' , tt1) )
}
tt1 <- Sys.time()
tt1-tt0
rownames(simul) <- paste( 'MAF =', MAFseq )
colnames(simul) <- paste( 'bet =', betseq )
simul

saveRDS(object=simul , file='~/Sex ratio/SimulationSNP_FamilyStructure_male_First_Xtra0075_1000repeats_20251209.Rds' )

################################################################################


################################################################################
### Simulate for all offspring
# Number of repeats for each simulation
n <- 1000
# Set of parameters to investigate
MAFseq <- c(.0025,.005,.01,.025,.05,.1)
betseq <- seq(.01,.21,by=.02)

# Total number of repeats
n*length(MAFseq)*length(betseq)

### Run - One offspring
simul <- matrix(NA,length(MAFseq),length(betseq))
tt0 <- Sys.time()
for(i in 1:length(MAFseq) ){
  for(j in 1:length(betseq) ){
    simul[i,j] <- mean( unlist( mclapply( X=rep(1,n) , FUN=famfun , MAF=MAFseq[i] , bet=betseq[j] , SR=0.514 , first=F , mc.cores=50 ) ) < 0.05 )
    print( paste('MAF',MAFseq[i] ,', beta',betseq[j],'done.') )
  }
  print( paste('MAF',MAFseq[i] ,'done.') )
  print( paste('Time' , tt1) )
}
tt1 <- Sys.time()
tt1-tt0
rownames(simul) <- paste( 'MAF =', MAFseq )
colnames(simul) <- paste( 'bet =', betseq )
simul

saveRDS(object=simul , file='~/Sex ratio/SimulationSNP_FamilyStructure_male_All_1000repeats_20251209.Rds' )


### Adding more points since too crude steps and plot looks bad
# Number of repeats for each simulation
n <- 1000
# Set of parameters to investigate
MAFseq <- c(.075)
betseq <- seq(.01,.21,by=.02)

# Total number of repeats
n*length(MAFseq)*length(betseq)
#n*length(MAFseq)*length(betseq)*11.77/10 /20 /60

### Run - One offspring
simul <- matrix(NA,length(MAFseq),length(betseq))
tt0 <- Sys.time()
for(i in 1:length(MAFseq) ){
  for(j in 1:length(betseq) ){
    simul[i,j] <- mean( unlist( mclapply( X=rep(1,n) , FUN=famfun , MAF=MAFseq[i] , bet=betseq[j] , SR=0.514 , first=F , mc.cores=20 ) ) < 0.05 )
    print( paste('MAF',MAFseq[i] ,', beta',betseq[j],'done.') )
  }
  print( paste('MAF',MAFseq[i] ,'done.') )
  print( paste('Time' , tt1) )
}
tt1 <- Sys.time()
tt1-tt0
rownames(simul) <- paste( 'MAF =', MAFseq )
colnames(simul) <- paste( 'bet =', betseq )
simul

saveRDS(object=simul , file='~/Sex ratio/SimulationSNP_FamilyStructure_male_All_Xtra0075_1000repeats_20251209.Rds' )


################################################################################


################################################################################
################################################################################
############################## END OF FILE #####################################
################################################################################
################################################################################