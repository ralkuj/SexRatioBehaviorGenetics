################################################################################
###
### File: 01.SNP_simulations_Song&Zhang_20250225.R
### Author: Ralf Kuja-Halkola
### Purpose: Implement the power simulation by Song & Zhang with incorrect and 
###          correct p-value thresholds.
### Date: 2024-11-05
### Updated: 2025-02-25
###
################################################################################

############################################
### The references:
# Song&Zhang:
#   Song S, Zhang J. In search of the genetic variants of human sex ratio at birth: was Fisher wrong about sex ratio evolution? Proc Biol Sci. Oct 2024;291(2033):20241876. doi:10.1098/rspb.2024.1876
# Zietsch et al:
#   Zietsch BP, Walum H, Lichtenstein P, Verweij KJH, Kuja-Halkola R. No genetic contribution to variation in human offspring sex ratio: a total population study of 4.7 million births. Proc Biol Sci. 2020;287(1921):20192849. DOI: 10.1098/rspb.2019.2849
############################################

############################################
### Libraries used
library(parallel)

### Nr cores used
nrcores <- 50

############################################


################################################################################
### Simulations

############################################
### The simulation from Song&Zhang translated into R functions
### The code is adapted from accompanying code to Song & Zhang, accessed 2025-02-25:
### https://github.com/song88180/Human_sex_ratio/blob/main/02.2_Zietsch_power_analysis.ipynb
### Some comments below kept from Song&Zhang

### Function for calculating cousin pair probabilities for sex combinations
getprobfun <- function(beta=.1,MAF=.1,SR_WT=.52){
  # Calculate sex ratios based on the effect size
  SR_1M <- SR_WT + beta
  SR_2M <- SR_WT + 2*beta
  a <- MAF
  A <- 1-MAF
  # Genotype probabilities
  P_AA <- A^2
  P_Aa <- 2*A*a
  P_aa <- a^2
  # Proportion of the grandparents having genotype XX and XX. This is talking about the shared grandparents of the cousin pair.
  P_GP_AA_AA <- A^4
  P_GP_AA_Aa <- 4*(A^3)*a
  P_GP_AA_aa <- 2*(A^2)*(a^2)
  P_GP_Aa_Aa <- 4*(A^2)*(a^2)
  P_GP_Aa_aa <- 4*A*(a^3)
  P_GP_aa_aa <- a^4
  
  # Proportion of the parents having genotype XX and XX. This is talking about the parents of the cousin pair that are children of the grandparents.
  P_Sib_AA_AA <- P_GP_AA_AA + (0.5^2)*P_GP_AA_Aa + (0.5^4)*P_GP_Aa_Aa
  P_Sib_AA_Aa <- 0.5*P_GP_AA_Aa + (0.5^2)*P_GP_Aa_Aa
  P_Sib_AA_aa <- (0.5^3)*P_GP_Aa_Aa
  P_Sib_Aa_Aa <- (0.5^2)*P_GP_AA_Aa + P_GP_AA_aa + (0.5^2)*P_GP_Aa_Aa + (0.5^2)*P_GP_Aa_aa
  P_Sib_Aa_aa <- (0.5^2)*P_GP_Aa_Aa + 0.5*P_GP_Aa_aa
  P_Sib_aa_aa <- (0.5^4)*P_GP_Aa_Aa + (0.5^2)*P_GP_Aa_aa + P_GP_aa_aa
  
  # Proportion of the mothers having genotype XX and XX. This is talking about the mother of the cousin pair that are not necessarily the children of the grandparents.
  P_mothers_AA_AA = P_Sib_AA_AA * (1/4)*(1 + P_AA + P_AA + P_AA^2) + 
                    P_Sib_AA_Aa * (1/4)*(0 + 0 + P_AA + P_AA^2) + 
                    P_Sib_AA_aa * (1/4)*(0 + 0 + P_AA + P_AA^2) + 
                    P_Sib_Aa_Aa * (1/4)*(0 + 0 + 0 + P_AA^2) + 
                    P_Sib_Aa_aa * (1/4)*(0 + 0 + 0 + P_AA^2) + 
                    P_Sib_aa_aa * (1/4)*(0 + 0 + 0 + P_AA^2)
  
  P_mothers_AA_Aa = P_Sib_AA_AA * (1/4)*(0 + P_Aa + P_Aa + 2*P_Aa*P_AA) + 
                    P_Sib_AA_Aa * (1/4)*(1 + P_AA + P_Aa + 2*P_Aa*P_AA) + 
                    P_Sib_AA_aa * (1/4)*(0 + 0 + P_Aa + 2*P_Aa*P_AA) + 
                    P_Sib_Aa_Aa * (1/4)*(0 + P_AA + P_AA + 2*P_Aa*P_AA) + 
                    P_Sib_Aa_aa * (1/4)*(0 + 0 + P_AA + 2*P_Aa*P_AA) + 
                    P_Sib_aa_aa * (1/4)*(0 + 0 + 0 + 2*P_Aa*P_AA)
  
  P_mothers_AA_aa = P_Sib_AA_AA * (1/4)*(0 + P_aa + P_aa + 2*P_AA*P_aa) + 
                    P_Sib_AA_Aa * (1/4)*(0 + 0 + P_aa + 2*P_AA*P_aa) + 
                    P_Sib_AA_aa * (1/4)*(1 + P_AA + P_aa + 2*P_AA*P_aa) + 
                    P_Sib_Aa_Aa * (1/4)*(0 + 0 + 0 + 2*P_AA*P_aa) + 
                    P_Sib_Aa_aa * (1/4)*(0 + P_AA + 0 + 2*P_AA*P_aa) + 
                    P_Sib_aa_aa * (1/4)*(0 + P_AA + P_AA + 2*P_AA*P_aa)
  
  P_mothers_Aa_Aa = P_Sib_AA_AA * (1/4)*(0 + 0 + 0 + P_Aa^2) + 
                    P_Sib_AA_Aa * (1/4)*(0 + P_Aa + 0 + P_Aa^2) + 
                    P_Sib_AA_aa * (1/4)*(0 + 0 + 0 + P_Aa^2) + 
                    P_Sib_Aa_Aa * (1/4)*(1 + P_Aa + P_Aa + P_Aa^2) + 
                    P_Sib_Aa_aa * (1/4)*(0 + 0 + P_Aa + P_Aa^2) + 
                    P_Sib_aa_aa * (1/4)*(0 + 0 + 0 + P_Aa^2)
  
  P_mothers_Aa_aa = P_Sib_AA_AA * (1/4)*(0 + 0 + 0 + 2*P_Aa*P_aa) + 
                    P_Sib_AA_Aa * (1/4)*(0 + P_aa + 0 + 2*P_Aa*P_aa) + 
                    P_Sib_AA_aa * (1/4)*(0 + P_Aa + 0 + 2*P_Aa*P_aa) + 
                    P_Sib_Aa_Aa * (1/4)*(0 + P_aa + P_aa + 2*P_Aa*P_aa) + 
                    P_Sib_Aa_aa * (1/4)*(1 + P_Aa + P_aa + 2*P_Aa*P_aa) + 
                    P_Sib_aa_aa * (1/4)*(0 + P_Aa + P_Aa + 2*P_Aa*P_aa)
  
  P_mothers_aa_aa = P_Sib_AA_AA * (1/4)*(0 + 0 + 0 + P_aa^2) + 
                    P_Sib_AA_Aa * (1/4)*(0 + 0 + 0 + P_aa^2) + 
                    P_Sib_AA_aa * (1/4)*(0 + P_aa + 0 + P_aa^2) + 
                    P_Sib_Aa_Aa * (1/4)*(0 + 0 + 0 + P_aa^2) + 
                    P_Sib_Aa_aa * (1/4)*(0 + P_aa + 0 + P_aa^2) + 
                    P_Sib_aa_aa * (1/4)*(1 + P_aa + P_aa + P_aa^2)
  
  # Proportion of cousion pairs with different combinations of sexes:
  # Male-Male, Male-Female, Female-Male, Female-Female
  P_Cos_M_M = P_mothers_AA_AA * (SR_WT^2) + P_mothers_AA_Aa * (SR_WT*SR_1M) + 
              P_mothers_AA_aa * (SR_WT*SR_2M) + P_mothers_Aa_Aa * (SR_1M*SR_1M) + 
              P_mothers_Aa_aa * (SR_1M*SR_2M) + P_mothers_aa_aa * (SR_2M*SR_2M)
  
  P_Cos_M_F = P_mothers_AA_AA * (SR_WT*(1-SR_WT)) + 
              0.5 * P_mothers_AA_Aa * ( SR_WT*(1-SR_1M) + (1-SR_WT)*SR_1M ) + 
              0.5 * P_mothers_AA_aa * ( SR_WT*(1-SR_2M) + (1-SR_WT)*SR_2M ) + 
              P_mothers_Aa_Aa * (SR_1M*(1-SR_1M)) + 
              0.5 * P_mothers_Aa_aa * ( SR_1M*(1-SR_2M) + (1-SR_1M)*SR_2M ) + 
              P_mothers_aa_aa * (SR_2M*(1-SR_2M))
  
  P_Cos_F_M = P_Cos_M_F
  
  P_Cos_F_F = P_mothers_AA_AA * (1-SR_WT)^2 + P_mothers_AA_Aa * (1-SR_WT)*(1-SR_1M) + 
              P_mothers_AA_aa * (1-SR_WT)*(1-SR_2M) + P_mothers_Aa_Aa * (1-SR_1M)*(1-SR_1M) + 
              P_mothers_Aa_aa * (1-SR_1M)*(1-SR_2M) + P_mothers_aa_aa * (1-SR_2M)*(1-SR_2M)
  return(c(P_Cos_M_M,P_Cos_M_F,P_Cos_F_M,P_Cos_F_F))
}


# # Checking produced sexratio (off from original)
# testfunSR <- function(N=100,beta=.1,MAF=.1,SR_WT=.52){
#   fortest <- apply( rmultinom( n=N , size=rep(1,4) , prob=getprobfun(beta=beta,MAF=MAF,SR_WT=SR_WT) ) , 1 , sum )
#   return((fortest[1]*2+fortest[2]+fortest[3])/2/N )
# }
# testfunSR(1000000)
# # To get a decided sex ratio in offspring generation, subtract 2*beta*MAF
# testfunSR(N=1000000,beta=.1,MAF=.1,SR_WT=.52-2*.1*.1)

### The simulation function
sexratioFun <- function(N=100,beta=.1,MAF=.1,SR_WT=.52-2*.1*.1){
  fortest <- apply( rmultinom( n=N , size=rep(1,4) , prob=getprobfun(beta=beta,MAF=MAF,SR_WT=SR_WT) ) , 1 , sum )
  return(chisq.test(matrix(fortest,2,2))$p.value)
}

# Test using the found effect size and MAF in GWAS from Song&Zhang for sample size of Zietsch et al
sexratioFun(14000000,beta=0.097,MAF=0.0025,SR_WT=.514-2*.097*.0025)

############################################


############################################
### Testing the example from Song&Zhang, MAF=0.1 and beta=0.1
N <- 14015421
n <- 1000
tt0 <- Sys.time()
pvalsSong <- unlist(mclapply( X=rep(N,n) , FUN=sexratioFun , mc.cores=nrcores ,beta=0.1,MAF=0.1,SR_WT=.52 ))
tt1 <- Sys.time()
tt1-tt0
# The correct p-value threshold
mean( pvalsSong < 0.05 ) # 0.938
# The incorrect p-value threshold used by Song&Zhang
mean( pvalsSong < 5e-8 ) # 0.017

### Testing the found SNP from Song&Zhang, MAF=0.0025 and beta=0.097
N <- 14015421
n <- 1000
tt0 <- Sys.time()
pvalsSong <- unlist(mclapply( X=rep(N,n) , FUN=sexratioFun , mc.cores=nrcores ,beta=0.097,MAF=0.0025,SR_WT=.52 ))
tt1 <- Sys.time()
tt1-tt0
# The correct p-value threshold
mean( pvalsSong < 0.05 ) # 0.063
# The incorrect p-value threshold used by Song&Zhang
mean( pvalsSong < 5e-8 ) # 0.000

### The incorrect p-value threshold code in Song & Zhang code for calculating 
### the true positive rate:
# TPR = (np.array(P_list) < 5e-8).sum() / len(P_list)
# return TPR

############################################


############################################
### Parallel many MAF and beta for plotting
# Number of cousin pairs
N <- 14015421
# Number of repeats for each simulation
n <- 1000
# Sequence to simulate over
MAFseq <- seq(.0025,.1,len=21)
betseq <- seq(.005,0.2,len=21)

# Store correct p-value threshold and incorrect p-value threshold for each repeat
simulSong <- matrix(NA,length(MAFseq),length(betseq))
simulSong_e10m8 <- matrix(NA,length(MAFseq),length(betseq))
### Run the simulation (make sure you have allocated sufficient cores and time)
tt0 <- Sys.time()
set.seed(12341234)
for(i in 1:length(MAFseq) ){
  for(j in 1:length(betseq) ){
    pvals <- unlist( mclapply( X=rep(N,n) , FUN=sexratioFun , MAF=MAFseq[i] , bet=betseq[j] , SR_WT=.52 , mc.cores=nrcores ) ) 
    simulSong[i,j] <- mean( pvals < 0.05 )
    simulSong_e10m8[i,j] <- mean( pvals < 5e-8 )
}}
tt1 <- Sys.time()
tt1-tt0

rownames(simulSong) <- rownames(simulSong_e10m8) <- paste( 'MAF =', MAFseq )
colnames(simulSong) <- colnames(simulSong_e10m8) <- paste( 'bet =', betseq )
simulSong
simulSong_e10m8


saveRDS(object=simulSong , file='~/Sex ratio/SimulationSNP_1000repeats_SongOriginal_correct_pval_20241105.Rds' )
saveRDS(object=simulSong_e10m8 , file='~/Sex ratio/SimulationSNP_1000repeats_SongOriginal_wrong_pval_20241105.Rds' )
############################################


################################################################################




################################################################################
################################################################################
############################## END OF FILE #####################################
################################################################################
################################################################################
