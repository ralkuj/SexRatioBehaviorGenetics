################################################################################
###
### File: 03.SNP_simulations_Plotting_20251210.R
### Author: Ralf Kuja-Halkola
### Purpose: Plotting of simulations for sex ratio of one SNP, using Song&Zhang
###          and Zietsch et al settings.
### Date: 20241105
### Updated: 20251210, change Zietsch et al analyses to using real family structure.
###          20260422, commenting and cleaning up code.
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
### Repeat of Song&Zhang plot
simulSong_e10m8 <- readRDS(file='~/Sex ratio/SimulationSNP_1000repeats_SongOriginal_wrong_pval_20241105.Rds' )

### Create a heatmap plot
plotdata <- simulSong_e10m8
MAFplot <- c(0,as.numeric( substr(rownames(plotdata),7,17) ))
BETAplot <- c(0,as.numeric( substr(colnames(plotdata),7,17) ))
# Add 0 for beta=0, MAF=0
plotdata <- cbind( 0, rbind( 0 , plotdata ))

# The plot
pdf(file='~/Sex ratio/Original_incorrect_pval_20250214.pdf',width=18/2.54,height=16/2.54)
par(oma=c(1,1,1,1))
filled.contour(z=plotdata,x=MAFplot,y=BETAplot,
  levels=c(0,seq(0.1,1,by=.05)),xlim=c(0,.1),ylim=c(0,.2),axes=T,frame.plot=F,
  plot.axes={ axis(1,seq(0,.1,by=.01),labels=0:10) ; axis(2,seq(0,.2,by=.02)) } ,
  key.title=title(main = "Power",cex.main=.8) , plot.title = title(main='Song&Zhang, incorrect p-value threshold',xlab='MAF (%)' , cex.lab=1.4 ,
# Add points
    plot.axes={
      points(0.1,0.1 , pch=8 , col=1 , xpd=NA )
      text(0.1,0.092 , bquote('MAF=10%, '* beta * '=0.10') , cex=0.8 , xpd=NA , adj=1  )
      points(0.0025,0.097 , pch=8 , col=2 )
      text(0.0025,0.089 , bquote('MAF=0.25%, '* beta * '=0.097') , cex=0.8 , xpd=NA , col=2, adj=0  )
      }
    )
)
mtext(bquote(beta), side = 2, line = 3, cex = 1.4, las = 1 , font=2 )

dev.off()

############################################


############################################
### Repeat of Song&Zhang plot, but correct p-value threshold
simulSong <- readRDS(file='~/Sex ratio/SimulationSNP_1000repeats_SongOriginal_correct_pval_20241105.Rds' )

### Create a heatmap plot
plotdata <- simulSong
MAFplot <- c(0,as.numeric( substr(rownames(plotdata),7,17) ))
BETAplot <- c(0,as.numeric( substr(colnames(plotdata),7,17) ))
# Add 0 for beta=0, MAF=0
plotdata <- cbind( 0, rbind( 0 , plotdata ))

# The plot
pdf(file='~/Sex ratio/Original_correct_pval_20250214.pdf',width=18/2.54,height=16/2.54)
par(oma=c(1,1,1,1))
filled.contour(z=plotdata,x=MAFplot,y=BETAplot,
  levels=c(0,seq(.1,1,by=.05)),xlim=c(0,.1),ylim=c(0,.2),axes=T,frame.plot=F,
  plot.axes={ axis(1,seq(0,.1,by=.01),labels=0:10) ; axis(2,seq(0,.2,by=.02)) } ,
  key.title=title(main = "Power",cex.main=.8) , plot.title = title(main='Song&Zhang, corrected p-value threshold',xlab='MAF (%)' , cex.lab=1.4 ,
# Add points
    plot.axes={
      points(0.1,0.1 , pch=8 , col=1 , xpd=NA )
      text(0.1,0.092 , bquote('MAF=10%, '* beta * '=0.10') , cex=0.8 , xpd=NA , adj=1 )
      points(0.0025,0.097 , pch=8 , col=2 )
      text(0.0025,0.089 , bquote('MAF=0.25%, '* beta * '=0.097') , cex=0.8 , xpd=NA , col=2, adj=0  )
      }
    )
)
mtext(bquote(beta), side = 2, line = 3, cex = 1.4, las = 1 , font=2 )

dev.off()

############################################


############################################
### Simulation males using family structure - One offspring
simulSongMalesFirstFamily <- readRDS(file='~/Sex ratio/SimulationSNP_FamilyStructure_male_First_1000repeats_20251209.Rds' )
simulSongMalesFirstFamilyXtra <- readRDS(file='~/Sex ratio/SimulationSNP_FamilyStructure_male_First_Xtra0075_1000repeats_20251209.Rds' )

### Create a heatmap plot
plotdata <- rbind(simulSongMalesFirstFamily[1:5,],simulSongMalesFirstFamilyXtra,simulSongMalesFirstFamily[6,]) # Link them together, mind the ordering
rownames(plotdata) <- c(rownames(simulSongMalesFirstFamily)[1:5],rownames(simulSongMalesFirstFamilyXtra),rownames(simulSongMalesFirstFamily)[6])
MAFplot <- c(0,as.numeric( substr(rownames(plotdata),7,17) ))
BETAplot <- c(0,as.numeric( substr(colnames(plotdata),7,17) ))
# Add 0 for beta=0, MAF=0
plotdata <- cbind( 0, rbind( 0 , plotdata ))

# The plot
pdf(file='~/Sex ratio/MalesOnlyFamily_oneoffspring_20251010.pdf',width=18/2.54,height=16/2.54)
par(oma=c(1,1,1,1))
filled.contour(z=plotdata,x=MAFplot,y=BETAplot,
  levels=c(0,seq(.1,1,by=.05)),xlim=c(0,.1),ylim=c(0,.2),axes=T,frame.plot=F,
  plot.axes={ axis(1,seq(0,.1,by=.01),labels=0:10) ; axis(2,seq(0,.2,by=.02)) } ,
  key.title=title(main = "Power",cex.main=.8) , plot.title = title(main='Males only, one offspring \naccounting for family structure',xlab='MAF (%)' , cex.lab=1.4 ,
# Add points
    plot.axes={
      points(0.1,0.1 , pch=8 , col=1 , xpd=NA )
      text(0.1,0.092 , bquote('MAF=10%, '* beta * '=0.10') , cex=0.8 , xpd=NA , adj=1  )
      points(0.0025,0.097 , pch=8 , col=2 )
      text(0.0025,0.089 , bquote('MAF=0.25%, '* beta * '=0.097') , cex=0.8 , xpd=NA , col=2, adj=0  )
      }
    )
)
mtext(bquote(beta), side = 2, line = 3, cex = 1.4, las = 1 , font=2 )

dev.off()

############################################


############################################
### Simulation males using family structure - All offspring
simulSongMalesAllFamily <- readRDS(file='~/Sex ratio/SimulationSNP_FamilyStructure_male_All_1000repeats_20251209.Rds' )
simulSongMalesAllFamilyXtra <- readRDS(file='~/Sex ratio/SimulationSNP_FamilyStructure_male_All_Xtra0075_1000repeats_20251209.Rds' )

### Create a heatmap plot
plotdata <- rbind(simulSongMalesAllFamily[1:5,],simulSongMalesAllFamilyXtra,simulSongMalesAllFamily[6,]) # Link them together, mind the ordering
rownames(plotdata) <- c(rownames(simulSongMalesAllFamily)[1:5],rownames(simulSongMalesAllFamilyXtra),rownames(simulSongMalesAllFamily)[6])
MAFplot <- c(0,as.numeric( substr(rownames(plotdata),7,17) ))
BETAplot <- c(0,as.numeric( substr(colnames(plotdata),7,17) ))
# Add 0 for beta=0, MAF=0
plotdata <- cbind( 0, rbind( 0 , plotdata ))

# The plot
pdf(file='~/Sex ratio/MalesOnlyFamily_all_20251010.pdf',width=18/2.54,height=16/2.54)
par(oma=c(1,1,1,1))
filled.contour(z=plotdata,x=MAFplot,y=BETAplot,
  levels=c(0,seq(.1,1,by=.05)),xlim=c(0,.1),ylim=c(0,.2),axes=T,frame.plot=F,
  plot.axes={ axis(1,seq(0,.1,by=.01),labels=0:10) ; axis(2,seq(0,.2,by=.02)) } ,
  key.title=title(main = "Power",cex.main=.8) , plot.title = title(main='Males only, all offspring \naccounting for family structure',xlab='MAF (%)' , cex.lab=1.4 ,
# Add points
    plot.axes={
      points(0.1,0.1 , pch=8 , col=1 , xpd=NA )
      text(0.1,0.092 , bquote('MAF=10%, '* beta * '=0.10') , cex=0.8 , xpd=NA , adj=1  )
      points(0.0025,0.097 , pch=8 , col=2 )
      text(0.0025,0.089 , bquote('MAF=0.25%, '* beta * '=0.097') , cex=0.8 , xpd=NA , col=2, adj=0  )
      }
    )
)
mtext(bquote(beta), side = 2, line = 3, cex = 1.4, las = 1 , font=2 )

dev.off()

############################################


################################################################################
################################################################################
############################## END OF FILE #####################################
################################################################################
################################################################################