#devtools::install_github("stan-dev/cmdstanr")
#library()
library(rethinking)
library(jtools)
library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(ggplot2)
library(ppcor) ## To compute partial correlation
library(MTS)
set.seed(20211118)
bet.a = -.25
bet.b = .75  
lb.a =  .5
lb.binv = -bet.a/(bet.b*lb.a) ## Inverse of lb.b
lb.b <- lb.binv
rho.yb = bet.b + lb.binv*lb.a*bet.a
rho.ya = bet.a + lb.binv*lb.a*bet.b
N = 500
## Condition for lb.a 
# lb.aval =seq(0.5,1, length = 1000)
# lb.bval <- -bet.b/(bet.a*lb.aval)
# plot(lb.aval, 1/lb.bval)

SimMisconcept11 <- function(N, bet.a, bet.b, lb.a, lb.b){
  
  X  = rnorm(n=N) #belief Dieatary fat is ok 
  X <- X/sd(X)
  FatA <- rnorm(n=N, mean=lb.a*X, sd = sqrt(1 - lb.a^2)); #FatA = FatA/sd(FatA) #
  FatB <- rnorm(n=N, mean = lb.b*X, sd = sqrt(1 - lb.b^2)); #FatB = FatB/sd(FatB) #sqrt(1 - lb.b^2)
  
  var.y = bet.a^2 + bet.b^2 + 2*lb.a*lb.b*bet.a*(bet.b) 
  mu = bet.a*FatA + bet.b*FatB
  Y <- rnorm(N, mean = mu, sd = sqrt(1 - var.y)) # 1 );  #Y <- Y/sd(Y)
  
  #cat("\n---  Correlation values ---\n")
  
  #cat("Cor(Y, FATA)[empirical] - Cor(Y, FATA)[Theoritical] \n")
  c(cor(Y, FatA), rho.ya) 
  
  #cat("Cor(Y, FATB)[empirical] - Cor(Y, FATB)[Theoritical]\n")  
  c(cor(Y, FatB), rho.yb)
  
  #c(cor(X, FatA), cor(X, FatB))
  #c(cor(Y, FatA), rho.ya, cor(Y, FatB), rho.yb)
  #All_Models= list()
  mAB <- (lm(Y ~ FatA + FatB))
  mB <- (lm(Y ~ FatB))
  mA <- (lm(Y ~ FatA))
  
  # Res <- c( c((cbind(summary(mAB)$coefficients[,1], confint.lm(mAB)))), #vector of length 9 (mean - 2.5% - 97.5%)
  #           c((cbind(summary(mA)$coefficients[,1], confint.lm(mA)))), #vector of length 6 (mean - 2.5% - 97.5%)
  #           c((cbind(summary(mB)$coefficients[,1], confint.lm(mB)))), #vector of length 6 (mean - 2.5% - 97.5%)
  #           cor(FatA, Y), # crude correlations
  #           cor(FatB, Y), # crude correlations
  #           pcor(data.frame(Y = Y, FatB = FatB, FatA = FatA))$estimate[1,-1] # partial correlations vector of length 2
  #           )
  Parm = c(0.0, bet.a, bet.b)
  Res <- c( c( c(summary(mAB)$coefficients[,1], (confint.lm(mAB)[,2] > Parm)*(confint.lm(mAB)[,1] < Parm)) ), #vector of length 6 (mean - cover in the 95%(1) or mot (0) )
            c(summary(mA)$coefficients[,1], (confint.lm(mA)[,2] > Parm[-3])*(confint.lm(mA)[,1] < Parm[-3])), #(mean - cover in the 95%(1) or mot (0) )
            c(summary(mB)$coefficients[,1], (confint.lm(mB)[,2] > Parm[-2])*(confint.lm(mB)[,1] < Parm[-2])), #(mean - cover in the 95%(1) or mot (0) )
            cor(FatA, Y), # crude correlations
            cor(FatB, Y), # crude correlations
            pcor(data.frame(Y = Y, FatA = FatA, FatB = FatB))$estimate[1,-1] # partial correlations vector of length 2
            #pcor(data.frame(Y = Y, FatA = FatA, FatB = FatB, X = X))$estimate[1,-1]
  )
  
  (Res)
}


SimMisconcept11(N, bet.a, bet.b, lb.a, lb.b)

ResOUt <- t(sapply(rep(N, 5000), SimMisconcept11,  bet.a = bet.a, bet.b =bet.b, lb.a=lb.a, lb.b=lb.b))
colNam <- colnames(ResOUt)
colNam[15:18] <- c("rho.YFatA", "rho.YFatB", "prho.YFatA", "prho.YFatB")
colnames(ResOUt) <- colNam


cbind(colMeans(ResOUt), t(apply(ResOUt,2, quantile, probs =c(0.025, .975))))



#pcor(data.frame(Y = Y, FatB = FatB, FatA = FatA, X =X))$estimate
#pcor(data.frame(Y = Y, FatB = FatB, FatA = FatA))$estimate

#cor(data.frame(Y = Y, FatB = FatB, FatA = FatA, X =X))
#--- Generate summaries


#SimMisconceptPlot <- function(N, bet.a, bet.b, lb.a, lb.b){
  
  X  = rnorm(n=N) #belief Dieatary fat is ok 
  X <- X/sd(X)
  FatA <- rnorm(n=N, mean=lb.a*X, sd = sqrt(1 - lb.a^2)); #FatA = FatA/sd(FatA) #
  FatB <- rnorm(n=N, mean = lb.b*X, sd = sqrt(1 - lb.b^2)); #FatB = FatB/sd(FatB) #sqrt(1 - lb.b^2)
  
  var.y = bet.a^2 + bet.b^2 + 2*lb.a*lb.b*bet.a*(bet.b) 
  mu = bet.a*FatA + bet.b*FatB
  Y <- rnorm(N, mean = mu, sd = sqrt(1 - var.y)) # 1 );  #Y <- Y/sd(Y)
  
  cat("\n---  Correlation values ---\n")
  
  cat("Cor(Y, FATA)[empirical] - Cor(Y, FATA)[Theoritical] \n")
  c(cor(Y, FatA), rho.ya) 
  
  cat("Cor(Y, FATB)[empirical] - Cor(Y, FATB)[Theoritical]\n")  
  c(cor(Y, FatB), rho.yb)
  
  #c(cor(X, FatA), cor(X, FatB))
  #c(cor(Y, FatA), rho.ya, cor(Y, FatB), rho.yb)
  #All_Models= list()
  mAB <- (lm(Y ~ FatA + FatB))
  mB <- (lm(Y ~ FatB))
  mA <- (lm(Y ~ FatA))
  plot_models(mAB,mA,mB, show.values = T)
  
  Res <- c( c((cbind(summary(mAB)$coefficients[,1], confint.lm(mAB)))), #vector of length 9 (mean - 2.5% - 97.5%)
            c((cbind(summary(mA)$coefficients[,1], confint.lm(mA)))), #vector of length 6 (mean - 2.5% - 97.5%)
            c((cbind(summary(mB)$coefficients[,1], confint.lm(mB)))), #vector of length 6 (mean - 2.5% - 97.5%)
            cor(FatA, Y), # crude correlations
            cor(FatB, Y), # crude correlations
            pcor(data.frame(Y = Y, FatB = FatB, FatA = FatA))$estimate[1,-1] # partial correlations vector of length 2
  )
  
  Res
}


#--- No ill formeed problem
#library(rethinking)
#cbind(lm(Y ~ FatA + FatB)$coefficients, confint.lm(lm(Y ~ FatA + FatB)))
#cbind(lm(Y ~ FatB)$coefficients, confint.lm(lm(Y ~  FatB))  )
#cbind(lm(Y ~ FatA)$coefficients, confint.lm(lm(Y ~  FatA))  )
#plot_models(mAB,mA,mB, show.values = T) #, value.offset=.3)

# summary(lm(Y ~ FatA + FatB))
# summary(lm(Y ~ FatB))
# summary(lm(Y~ FatA))
# coeftab(lm(Y~ FatA), lm(Y~ FatB), lm(Y~ FatA + FatB))
#plot( coeftab( lm(Y~ FatA), lm(Y~ FatA + FatB)), pars =c("FatA"))
#plot( coeftab( lm(Y~ FatB), lm(Y~ FatA + FatB)), pars =c("FatB"))

#plot_summs(lm(Y~ FatA), lm(Y~ FatB), lm(Y~ FatA + FatB), scale=T,model.names = c("Y~FATA","Y~FATB", "Y~FATA+FATB"))
#coeftab_plot( coeftab( lm(Y~ FatA), lm(Y~ FatB), lm(Y~ FatA + FatB) ), pars ="M")



#--- case of Misconception 10

library(rethinking)
library(jtools)
SimRo <- function(taucut){
  #runif(1,0,1)
  set.seed(2008)
  N = 1000000
  bet.z <- .5     #runif(1,0,1)
  bet.x <- .45      #runif(1,0,1)
  lb.z <- .4
  al.x <- 0.0     #runif(1,0,1)
  al.z <- .6     #runif(1,0,1)
  
  #--- simulate
  Z <- rnorm(n=N)
  #--- Simulate X 
  varx.z <- 1 - lb.z^2
  X = lb.z*Z + rnorm(n=N, mean = 0, sd = sqrt(varx.z))
  #X <- X/sd(X) 
  #--- Simulate Y2  
  vary2.z = 1 - (al.z^2 + al.x^2)
  Y2 <- al.z*Z + al.x*X + rnorm(n=N, mean = 0, sd = sqrt(vary2.z))
  #Y2 <- Y2/sd(Y2)
  #-- Simulate Y1
  vary1.xz <- 1 - (bet.z^2 +  (bet.x^2)*(lb.z^2) + 2*bet.z*bet.x*lb.z + (bet.x^2))
  Y1 <-  bet.z*Z + bet.x*X + rnorm(n=N, mean = 0, sd= sqrt(vary1.xz))
  #Y1 <- Y1/sd(Y1)
  
  id <- which(Y1 > 2.366733) #taucut*sd(Y1) + mean(Y1))
  R0 <- lm(Y1[id] ~  X[id])
  
  R00 <- lm(Y1[-id] ~  X[-id])
  
  R1 <- lm(Y1[id] ~  X[id] + Z[id])
  R10 <- lm(Y1[-id] ~  X[-id] + Z[-id])
  
  R20 <- lm(Y2 ~ X)
  R21 <- lm(Y2 ~ Z)
  R22 <- lm(Y2 ~ X + Z)
  
  R2id0 <- lm(Y2[id] ~ X[id] )
  R2id1 <- lm(Y2[id] ~ Z[id])
  R2id2 <- lm(Y2[id] ~ X[id] + Z[id])
  
  
  R2idn0 <- lm(Y2[-id] ~ X[-id])
  R2idn1 <- lm(Y2[-id] ~ Z[-id])
  R2idn2 <- lm(Y2[-id] ~ X[-id] + Z[-id])
  
  Res <- c(mean(Y1), sd(Y1), cor(X, Z) ,cor(X[id], Z[id]), cor(X[-id], Z[-id]), 
           as.numeric(R0$coefficients)[2], as.numeric(R00$coefficients)[2],
           as.numeric(R1$coefficients)[2:3],as.numeric(R10$coefficients)[2:3],
           as.numeric(R20$coefficients)[2], as.numeric(R21$coefficients)[2],as.numeric(R22$coefficients)[2:3],
           as.numeric(R2id0$coefficients)[2], as.numeric(R2id1$coefficients)[2], as.numeric(R2id2$coefficients)[2:3],
           as.numeric(R2idn0$coefficients)[2], as.numeric(R2idn1$coefficients)[2], as.numeric(R2idn2$coefficients)[2:3])
  
}   
SimRo <- Vectorize(SimRo,"taucut")#

#--- Run code
SimRo(1.30)

taucut <- seq(.1, 2, length.out = 500)
Res <- t(SimRo(taucut))
colnames(Res)  <- c("meanY1","sdY1", "rhoXZ(all.)" ,"rhoXZ(abov.)","rhoXZ(blw.)",
                    "betaxhat_Y1X(abov.)", "betaxhat_Y1X(blw.)",
                    "betaxhat_Y1XZ(abov)","betazhat_Y1Xz(blw)", "betaxhat_Y1XZ(blw)","betazhat_Y1Xz(blw)",
                    "alphax_Y2X(all)", "alphaz_Y2Z(all)", "alphax_Y2XZ(all)", "alphaz_Y2XZ(all)",
                    "alphax_Y2X(abov)", "alphaz_Y2Z(abov)", "alphax_Y2XZ(abov)", "alphaz_Y2XZ(abov)",
                    "alphax_Y2X(blw)","alphaz_Y2Z(blw)", "alphax_Y2XZ(blw)" , "alphaz_Y2XZ(blw)")


#--- Find the zero values
Idp <- which(Res[,4] > 0)
IdN <- which(Res[,4] < 0)

cat("\n ----- Chosen case ------ \n")
Res[c(Idp[length(Idp)],IdN[1]),]
taucut[c(Idp[length(Idp)],IdN[1])]
taucut0 = taucut[IdN[1]]
taucut0

Res[c(IdN[1]),]
#Res0Tep <- SimRoPlot(taucut[IdN[1]])
#dim(Res0Tep$Mod)

#coeftab(lm(Res0Tep$Mod[,1] ~ Res0Tep$Mod[,2]) )
plot(taucut,Res[,4], type="p", ylab = "cor(X,Z)")
abline(h=0, lwd=4, col="orange")
abline(v=taucut0, lwd=4, col="blue")


#--- Type error estimates
CorXZ <- function(a, Y1, X, Z){
  
  id <- which(Y1 > mean(Y1) + sd(Y1)*a)
  cor(X[id], Z[id]) 
}

CorXZ <- Vectorize(CorXZ,"a")



#--- Looking at Y1
SimMisconcept12 <- function(N){
  
  bet.z <- .5     #runif(1,0,1)
  bet.x <- .45      #runif(1,0,1)
  lb.z <- .4
  al.x <- 0.0     #runif(1,0,1)
  al.z <- .6     #runif(1,0,1)
  
  #--- simulate
  Z <- rnorm(n=N)
  #--- Simulate X 
  varx.z <- 1 - lb.z^2
  X = lb.z*Z + rnorm(n=N, mean = 0, sd = sqrt(varx.z))
  #X <- X/sd(X) 
  #--- Simulate Y2  
  vary2.z = 1 - al.z^2
  Y2 <- al.z*Z + rnorm(n=N, mean = 0, sd = sqrt(vary2.z))
  #Y2 <- Y2/sd(Y2)
  #-- Simulate Y1
  vary1.xz <- 1 - (bet.z^2 +  (bet.x^2)*(lb.z^2) + 2*bet.z*bet.x*lb.z + (bet.x^2))
  Y1 <-  bet.z*Z + bet.x*X + rnorm(n=N, mean = 0, sd= sqrt(vary1.xz))
  #Y1 <- Y1/sd(Y1)
  # CutVec <- seq(.1,2,length.out =300)
  # CoVec <- CorXZ(CutVec, Y1, X, Z)
  # Idp <- which(CoVec > 0)
  # IdN <- which(CoVec < 0)
  # 
  # id00 <- c(Idp[length(Idp)],IdN[1])[which.min(abs(CoVec[c(Idp[length(Idp)],IdN[1])]))]
  # taucut0 = CutVec[id00]
  taucut0 = 2.366733
  id <- which(Y1 >  2.366733)#taucut0*sd(Y1) + mean(Y1))
  
  #cat("\n---  Correlation values ---\n")
  
  #cat("Cor(Y, FATA)[empirical] - Cor(Y, FATA)[Theoritical] \n")
  #c(cor(Y, FatA), rho.ya) 
  
  #cat("Cor(Y, FATB)[empirical] - Cor(Y, FATB)[Theoritical]\n")  
  #c(cor(Y, FatB), rho.yb)
  
  #c(cor(X, FatA), cor(X, FatB))
  #c(cor(Y, FatA), rho.ya, cor(Y, FatB), rho.yb)
  #All_Models= list()
  mAB <- (lm(Y1[id] ~ X[id] + Z[id]))
  mB <- (lm(Y1[id] ~ X[id]))
  mA <- (lm(Y1[id] ~ Z[id]))
  
  # Res <- c( c((cbind(summary(mAB)$coefficients[,1], confint.lm(mAB)))), #vector of length 9 (mean - 2.5% - 97.5%)
  #           c((cbind(summary(mA)$coefficients[,1], confint.lm(mA)))), #vector of length 6 (mean - 2.5% - 97.5%)
  #           c((cbind(summary(mB)$coefficients[,1], confint.lm(mB)))), #vector of length 6 (mean - 2.5% - 97.5%)
  #           cor(FatA, Y), # crude correlations
  #           cor(FatB, Y), # crude correlations
  #           pcor(data.frame(Y = Y, FatB = FatB, FatA = FatA))$estimate[1,-1] # partial correlations vector of length 2
  #           )
  Parm = c(0.0, bet.x, bet.z)
  Res <- c( c(cor(X[id], Z[id]),taucut0, c(summary(mAB)$coefficients[,1], (confint.lm(mAB)[,2] > Parm)*(confint.lm(mAB)[,1] < Parm)) ), #vector of length 6 (mean - cover in the 95%(1) or mot (0) )
            c(summary(mA)$coefficients[,1], (confint.lm(mA)[,2] > Parm[-3])*(confint.lm(mA)[,1] < Parm[-3])), #(mean - cover in the 95%(1) or mot (0) )
            c(summary(mB)$coefficients[,1], (confint.lm(mB)[,2] > Parm[-2])*(confint.lm(mB)[,1] < Parm[-2])), #(mean - cover in the 95%(1) or mot (0) )
            cor(X[id], Y1[id]), # crude correlations
            cor(Z[id], Y1[id]), # crude correlations
            pcor(data.frame(Y1 = Y1[id], X = X[id], Z = Z[id]))$estimate[1,-1] # partial correlations vector of length 2
            #pcor(data.frame(Y = Y, FatA = FatA, FatB = FatB, X = X))$estimate[1,-1]
  )
  
  (Res)
}

SimMisconcept12(500000)

#ResOUt_500 <- t(sapply(rep(500, 10000), SimMisconcept12))
#ResOUt_1000 <- t(sapply(rep(1000, 10000), SimMisconcept12))
ResOUt_2000 <- t(sapply(rep(500000, 10000), SimMisconcept12))


Temp <- ResOUt_2000 
colNam <- colnames(Temp)
colNam[c(1:2,17:20)] <- c("rho.ZX[Y1>tau]", "Tau", "rho.Y1X", "rho.Y1Z", "prho.Y1X", "prho.Y1Z")
colnames(Temp) <- colNam

#head(ResOUt)
Temp <- (Temp[abs(Temp[,1]) < 0.0006, ]);
dim(Temp)


xtable::xtable(round(cbind(colMeans(Temp), t(apply(Temp,2, quantile, probs =c(0.025, .975)))), 4))

#---------------------------------------------------
#--- Looking that Y2
#---------------------------------------------------

SimMisconcept12 <- function(N){
  
  bet.z <- .5     #runif(1,0,1)
  bet.x <- .45      #runif(1,0,1)
  lb.z <- .4
  al.x <- 0.     #runif(1,0,1)
  al.z <- .6     #runif(1,0,1)
  
  #--- simulate
  Z <- rnorm(n=N)
  #--- Simulate X 
  varx.z <- 1 - lb.z^2
  X = lb.z*Z + rnorm(n=N, mean = 0, sd = sqrt(varx.z))
  #X <- X/sd(X) 
  #--- Simulate Y2  
  vary2.z = 1 - al.z^2
  Y2 <- al.z*Z + rnorm(n=N, mean = 0, sd = sqrt(vary2.z))
  #Y2 <- Y2/sd(Y2)
  #-- Simulate Y1
  vary1.xz <- 1 - (bet.z^2 +  (bet.x^2)*(lb.z^2) + 2*bet.z*bet.x*lb.z + (bet.x^2))
  Y1 <-  bet.z*Z + bet.x*X + rnorm(n=N, mean = 0, sd= sqrt(vary1.xz))
  #Y1 <- Y1/sd(Y1)
  # CutVec <- seq(.1,2,length.out =300)
  # CoVec <- CorXZ(CutVec, Y1, X, Z)
  # Idp <- which(CoVec > 0)
  # IdN <- which(CoVec < 0)
  # 
  # id00 <- c(Idp[length(Idp)],IdN[1])[which.min(abs(CoVec[c(Idp[length(Idp)],IdN[1])]))]
  # taucut0 = CutVec[id00]
  OpRes <- optim(1.6, CorrelValZX, al.z=al.z, bet.z = bet.z, bet.x = bet.x, lb.z = lb.z)
  #taucut0 = 2.366733
  taucut0 = OpRes$par
  id <- which(Y1 > taucut0 )#taucut0*sd(Y1) + mean(Y1))
  
  #cat("\n---  Correlation values ---\n")
  
  #cat("Cor(Y, FATA)[empirical] - Cor(Y, FATA)[Theoritical] \n")
  #c(cor(Y, FatA), rho.ya) 
  
  #cat("Cor(Y, FATB)[empirical] - Cor(Y, FATB)[Theoritical]\n")  
  #c(cor(Y, FatB), rho.yb)
  
  #c(cor(X, FatA), cor(X, FatB))
  #c(cor(Y, FatA), rho.ya, cor(Y, FatB), rho.yb)
  #All_Models= list()
  mAB <- (lm(Y2[id] ~ X[id] + Z[id]))
  mB <- (lm(Y2[id] ~ X[id]))
  mA <- (lm(Y2[id] ~ Z[id]))
  
  # Res <- c( c((cbind(summary(mAB)$coefficients[,1], confint.lm(mAB)))), #vector of length 9 (mean - 2.5% - 97.5%)
  #           c((cbind(summary(mA)$coefficients[,1], confint.lm(mA)))), #vector of length 6 (mean - 2.5% - 97.5%)
  #           c((cbind(summary(mB)$coefficients[,1], confint.lm(mB)))), #vector of length 6 (mean - 2.5% - 97.5%)
  #           cor(FatA, Y), # crude correlations
  #           cor(FatB, Y), # crude correlations
  #           pcor(data.frame(Y = Y, FatB = FatB, FatA = FatA))$estimate[1,-1] # partial correlations vector of length 2
  #           )
  Parm = c(0.0, al.x, al.z)
  Res <- c( c(cor(X[id], Z[id]),taucut0, c(summary(mAB)$coefficients[,1], (confint.lm(mAB)[,2] > Parm)*(confint.lm(mAB)[,1] < Parm)) ), #vector of length 6 (mean - cover in the 95%(1) or mot (0) )
            c(summary(mA)$coefficients[,1], (confint.lm(mA)[,2] > Parm[-3])*(confint.lm(mA)[,1] < Parm[-3])), #(mean - cover in the 95%(1) or mot (0) )
            c(summary(mB)$coefficients[,1], (confint.lm(mB)[,2] > Parm[-2])*(confint.lm(mB)[,1] < Parm[-2])), #(mean - cover in the 95%(1) or mot (0) )
            cor(X[id], Y2[id]), # crude correlations
            cor(Z[id], Y2[id]), # crude correlations
            pcor(data.frame(Y2 = Y2[id], X = X[id], Z = Z[id]))$estimate[1,-1] # partial correlations vector of length 2
            #pcor(data.frame(Y = Y, FatA = FatA, FatB = FatB, X = X))$estimate[1,-1]
  )
  
  (Res)
}


CorrelValZX <- function(al.z, bet.z, bet.x, lb.z, a){
  
  Sig <- CovMat(al.z, bet.z, bet.x, lb.z)
  m1 = MTS::msqrt(Sig[1:2,1:2])
  
  
  #---------------
  Al = c(1/sqrt( 1 - (t(Sig[1:2,3])%*%solve(Sig[1:2,1:2])%*%(Sig[1:2,3]) ) ))*solve(Sig[1:2,1:2])%*%Sig[1:2,3]
  Del = (1/sqrt(c(1 + t(Al)%*%Sig[1:2,1:2]%*%Al)))*Sig[1:2,1:2]%*%Al
  
  (1/sqrt( 1 - (t(Sig[1:2,3])%*%solve(Sig[1:2,1:2])%*%(Sig[1:2,3]) ) ))*(1.6)
  1.6*sqrt(c(1 + t(Al)%*%Sig[1:2,1:2]%*%Al))
  
  #delta = c(msqrt(1 + t(Lb1)%*%Sig[1:2,1:2]%*%(Lb1))$invsqrt)*Sig[1:2,1:2]%*%Lb1
  mu = (Del*(dnorm(-a)/pnorm(-a)))
  
  # Sig[1:2,1:2] - diag(1/c(mu))%*%(tcrossprod(mu))%*%diag(1/c(mu)) #m1$invsqrt%*%(tcrossprod(mu))%*%m1$invsqrt
  SigF =  Sig[1:2,1:2] + (-(dnorm(-a)/pnorm(-a))^2 - (-a*dnorm(-a)/pnorm(-a)))*(tcrossprod(Del)) #m1$invsqrt%*%(tcrossprod(mu))%*%m1$invsqrt
  
  Res = list(Mu = mu, Sig = SigF, CorSig = cov2cor(SigF)) 
  (SigF[1,2])^2
  
}


CovMat <- function(al.z, bet.z, bet.x, lb.z){
  
  Sig <- matrix(0, nrow=3,ncol=3)
  Sig[1,2] = Sig[2,1] =  lb.z  # cor(X,Z)
  Sig[1,3] = Sig[3,1] = bet.z + lb.z*bet.x # cor(Y1, Z)
  Sig[2,3] = Sig[3,2] = bet.x + bet.z*lb.z  #cor(Y1, X)
  diag(Sig) <- 1
  
  Sig
  
}

CovMat(al.z, bet.z, bet.x, lb.z)


optim(1.6, CorrelValZX, al.z=al.z, bet.z = bet.z, bet.x = bet.x, lb.z = lb.z)

SimMisconcept12(50000)

#ResOUt2_500 <- t(sapply(rep(500, 10000), SimMisconcept12))
#ResOUt2_1000 <- t(sapply(rep(1000, 10000), SimMisconcept12))
ResOUt2_2000 <- t(sapply(rep(500000, 10000), SimMisconcept12))


Temp <- ResOUt2_2000 
colNam <- colnames(Temp)
colNam[c(1:2,17:20)] <- c("rho.ZX[Y1>tau]", "Tau", "rho.Y2X", "rho.Y2Z", "prho.Y2X", "prho.Y2Z")
colnames(Temp) <- colNam

#head(ResOUt)
sum(abs(Temp[,1]) < 0.0006)
Temp <- (Temp[abs(Temp[,1]) < 0.0006, ]);
dim(Temp)

round(cbind(colMeans(Temp), t(apply(Temp,2, quantile, probs =c(0.025, .975)))), 4)



#----------------------------------------------------------------------
#------ Including Y2 
#----------------------------------------------------------------------
bet.z <- .5     #runif(1,0,1)
bet.x <- .45      #runif(1,0,1)
lb.z <- .4  # old
al.x <- 0.0     #runif(1,0,1)
al.z <- .6     #runif(1,0,1)

#--- New
bet.z <- .2   #runif(1,0,1)
bet.x <- .6      #runif(1,0,1)
lb.z <- .4  # old


#--- CHoosing paramaters

ParmSelc <- function(lb.z){
  
  betax <- c(0, sqrt(1/(1+lb.z^2)), runif(n=1, 0, sqrt(1/(1+lb.z^2))) )
  betaz <- c(0, sqrt(1- betax[3]^2) - lb.z*betax[3], runif(n=1, 0, sqrt(1- betax[3]^2) - lb.z*betax[3]))
  cbind(betax, betaz, (betaz[3] + lb.z*betax[3])*(betax[3] + lb.z*betaz[3]))
}

ParmSelc(.4)

ParmSelc <- function(lb.z, a,b){
  
  betax <-  a*sqrt(1/(1+lb.z^2))
  betaz <- (b/a)* betax*(sqrt(1+lb.z*lb.z - a*a) - lb.z*a)
  c(betax, betaz, lb.z, (betaz + lb.z*betax)*(betax + lb.z*betaz) - lb.z, lb.z/((betaz + lb.z*betax)*(betax + lb.z*betaz)) )
}
ParmSelc <- Vectorize(ParmSelc,"lb.z")


ParmSelcMin <- function(lb.z, a,b){
  
  betax <-  a*sqrt(1/(1+lb.z^2))
  betaz <- (b/a)* betax*(sqrt(1+lb.z*lb.z - a*a) - lb.z*a)
  res = c(betax, betaz, lb.z, (betaz + lb.z*betax)*(betax + lb.z*betaz) - lb.z, lb.z/((betaz + lb.z*betax)*(betax + lb.z*betaz)) )
  (res[4])^2
}
ParmSelcMin <- Vectorize(ParmSelcMin,"lb.z")

a <- c(.2,.9)
b <- c(.2,.9)


lb.z0 <- seq(0, 1, length.out = 1000)
par(mfrow=c(2,2), font=12, cex.axis=1.6, font.axis=12, cex.lab=1.3, font.lab= 12)

LambdaZFin <- c(0.02, .035, 0.27,.4, 0.12, 0.19, 0.26, 0.45)  ## These are the value of lb.z we will use

cp =0
names <-c("N","betx","betz","lbz","tau",
          "corSel",
          "corFull",
          "PcorSel",
          "PcorFull",
          
          "mb_1","mb_2",
          "Fmb_1","Fmb_2",
          "mab_1","mab_2", "mab_3",
          "Fmab_1","Fmab_2","Fmab_3",
          #texting
          "mab_intp","mab_alx","mab_alz","mab_pb1","mab_pb2","mab_pb3", 
          "mb_intp","mb_alx","mb_pb1","mb_pb2",
          "Fmab_intp","Fmab_alx","Fmab_alz","Fmab_pb1","Fmab_pb2","Fmab_pb3", 
          "Fmb_intp","Fmb_alx","Fmb_pb1","Fmb_pb2")
Final_Res <- list() 
OutFinal <- NULL
for(i in 1:length(a)){
  Final_Res[[i]] <- list()
  for(j in 1:length(b)){
    Final_Res[[i]][[j]] <- list()
    #-Cat generate plots
    Re0 <- optimize(ParmSelcMin, c(0, .99) , a= a[i], b= b[j])
    lb.z00 <- Re0$minimum
    V00 <- ParmSelc(lb.z0, a[i], b[j])
    plot(lb.z0, V00[4,], xlab = expression(lambda[z]), ylab= "Objective", main = paste("a=",a[i],";","b=",b[j], sep=""), type="l", lwd=3)
    abline(h = 0,lwd=2)
    head(V00[,1:10])
    x0 <- lb.z0[lb.z0 < lb.z00]
    y0<- V00[4, lb.z0 < lb.z00] 
    polygon(c(x0, rev(x0)), c(y0 ,rep(0, length(y0))),
            col = "#6BD7AF")
    abline(v =lb.z00, col="grey",lwd=3,lty=2)
    abline(v =0, col="grey",lwd=3,lty=2)
    #- Perform Simulation Choosing two cases
    #-- case
    cp <- cp+1
    Parm0 <- ParmSelc(LambdaZFin[c((2*(cp-1)+1):(2*cp))], a[i],b[j])
    for(m in 1:2){
      cat("case i, j m===>", i,j, m,"\n")
      Final_Res[[i]][[j]][[m]] <- list()
      bet.x = Parm0[1,m]
      bet.z = Parm0[2,m]
      lb.z <- Parm0[3,m] 
      OT <- t(replicate(1000, SimMisconcept12V2(N, bet.x, bet.z, lb.z)))
      colnames(OT) <- names
      
      
      Final_Res[[i]][[j]][[m]] <- c(a=a[i], b= b[j], lb.z00,colMeans(OT))
      OutFinal <- rbind(OutFinal, c(a=a[i], b= b[j], lb.z00,colMeans(OT)))
    }
  }
}


#---- Summarise the simulation 
rbind(Final_Res[[1]][[1]][[1]],Final_Res[[1]][[1]][[2]])
OutFinal



V00 <- ParmSelc(lb.z0, .2,.2)
par(mfrow=c(2,2), font=12, cex=1.5)
plot(lb.z0, V00[4,],xlab = expression(lambda[z]))
abline(h = 0,lwd=3)
head(V00[,1:10])


plot(lb.z0, V00[5,])
abline(h = 1,lwd=3)
head(V00[,1:10])
a = seq(-6,6,length.out=1000)
Y00 <- (-(dnorm(-a)/pnorm(-a))^2 - (-a*dnorm(-a)/pnorm(-a)))
plot(a,-Y00)

which.min(abs(V00[4,]))
id = 400
V00[,id]
bet.x <- V00[1,id]
bet.z <- V00[2,id]
lb.z <- V00[3,id]

#@SImulate data

SimMisconcept12V2 <- function(N, bet.x, bet.z, lb.z){
  
  #bet.z <- .5     #runif(1,0,1)
  #bet.x <- .45      #runif(1,0,1)
  #lb.z <- .4
  al.x <- 0.     #runif(1,0,1)
  al.z <- .6     #runif(1,0,1)
  
  #--- simulate
  Z <- rnorm(n=N)
  #--- Simulate X 
  #varx.z <- 1 - lb.z^2
  varx.z <- 1 - var(Z)*lb.z^2
  X = lb.z*Z + rnorm(n=N, mean = 0, sd = sqrt(varx.z))
  X <- X/sd(X) 
  #--- Simulate Y2  
  #vary2.z = 1 - al.z^2
  vary2.z = 1 - (var(Z)*al.z^2 + (al.x^2)*var(X) + 2*(al.x)*al.z*cov(X, Z))
  Y2 <- al.z*Z + rnorm(n=N, mean = 0, sd = sqrt(vary2.z))
  Y2 <- Y2/sd(Y2)
  #-- Simulate Y1
  #vary1.xz <- 1 - (bet.z^2 +  (bet.x^2)*(lb.z^2) + 2*bet.z*bet.x*lb.z + (bet.x^2))
  vary1.xz <- 1 - (var(Z)*bet.z^2 +  var(X)*(bet.x^2) + 2*bet.z*bet.x*cov(X, Z))
  Y1 <-  bet.z*Z + bet.x*X + rnorm(n=N, mean = 0, sd= sqrt(vary1.xz))
  Y1 <- Y1/sd(Y1)
  # CutVec <- seq(.1,2,length.out =300)
  # CoVec <- CorXZ(CutVec, Y1, X, Z)
  # Idp <- which(CoVec > 0)
  # IdN <- which(CoVec < 0)
  # 
  # FInd cut-off pount a for selecting on Y1
  #OpRes <- optim(1.6, CorrelValZX, al.z=al.z, bet.z = bet.z, bet.x = bet.x, lb.z = lb.z)
  OpRes = optimize(CorrelValObjV2, c(-7, 7), al.z=al.z, bet.z = bet.z, bet.x = bet.x, lb.z = lb.z)
  #taucut0 = 2.366733
  taucut0 = OpRes$minimum
  id <- which(Y1 > taucut0) #taucut0*sd(Y1) + mean(Y1))
  #All_Models= list()
  mAB <- (lm(Y2[id] ~ X[id] + Z[id])) #(4)
  mB <- (lm(Y2[id] ~ X[id])) # (3)
  ## mA <- (lm(Y2[id] ~ Z[id]))
  FmAB <-  (lm(Y2 ~ X + Z)) #(4)
  FmB <-  (lm(Y2 ~ X )) # (3)
  #FmAB <-  (lm(Y2 ~ X + Z))
  # Res <- c( c((cbind(summary(mAB)$coefficients[,1], confint.lm(mAB)))), #vector of length 9 (mean - 2.5% - 97.5%)
  #           c((cbind(summary(mA)$coefficients[,1], confint.lm(mA)))), #vector of length 6 (mean - 2.5% - 97.5%)
  #           c((cbind(summary(mB)$coefficients[,1], confint.lm(mB)))), #vector of length 6 (mean - 2.5% - 97.5%)
  #           cor(FatA, Y), # crude correlations
  #           cor(FatB, Y), # crude correlations
  #           pcor(data.frame(Y = Y, FatB = FatB, FatA = FatA))$estimate[1,-1] # partial correlations vector of length 2
  #           )
  Parm = c(0.0, al.x, al.z)
  #Out <- NULL
  
  Res <- c( N, bet.x, bet.z, lb.z, taucut0,
            cor(X[id], Y2[id]), # crude correlations (1)
            cor(X, Y2), # crude correlations   (1)
            pcor(data.frame(Y2 = Y2[id], X = X[id], Z = Z[id]))$estimate[1,2], #-1], # partial correlations vector of length 2 (2)
            pcor(data.frame(Y2 = Y2, X = X, Z = Z))$estimate[1, 2], #-1] # partial correlations vector of length 2  (2)
            
            c(summary(mB)$coefficients[,1] - Parm[-3]), #(3)
            c(summary(FmB)$coefficients[,1] - Parm[-3]), #(3)
            c(summary(mAB)$coefficients[,1] - Parm), # (4)
            c(summary(FmAB)$coefficients[,1] - Parm),# (4)
            #-- Testing
            c(c(summary(mAB)$coefficients[,1], (confint.lm(mAB)[,2] > Parm)*(confint.lm(mAB)[,1] < Parm)) ), #vector of length 6 (mean - cover in the 95%(1) or mot (0) )
            #c(summary(mA)$coefficients[,1], (confint.lm(mA)[,2] > Parm[-3])*(confint.lm(mA)[,1] < Parm[-3])), #(mean - cover in the 95%(1) or mot (0) )
            c(summary(mB)$coefficients[,1], (confint.lm(mB)[,2] > Parm[-3])*(confint.lm(mB)[,1] < Parm[-3])), #(mean - cover in the 95%(1) or mot (0) )
            
            c(summary(FmAB)$coefficients[,1], (confint.lm(FmAB)[,2] > Parm)*(confint.lm(FmAB)[,1] < Parm)),
            c(summary(FmB)$coefficients[,1], (confint.lm(FmB)[,2] > Parm[-3])*(confint.lm(FmB)[,1] < Parm[-3])) #(mean - cover in the 95%(1) or mot (0) )
            
            
            #pcor(data.frame(Y = Y, FatA = FatA, FatB = FatB, X = X))$estimate[1,-1]
  )
  #Out <- rbind(Out, Res)
  return(Res)
}

#@ estimate partial Slope - will take a covariance matrix with the last row and column for Y 
PartialSlope <- function(Cov){
  J = ncol(Cov)
  sdt = (sqrt(diag(Cov)))
  CorMat <- cov2cor(Cov)
  Bet <- numeric(J-1)
  Bet[1] = (sdt[J]/sdt[1])*(CorMat[1,3] - CorMat[2,3]*CorMat[1,2] )/(1 - CorMat[1,2]^2) ## Z
  Bet[2] = (sdt[J]/sdt[2])*(CorMat[2,3] - CorMat[1,3]*CorMat[1,2] )/(1 - CorMat[1,2]^2) ## for X
  
  Bet  
}

#@ Compute the inverse of the Correlation Matrix
InvMat3by3 <- function(Sig){
  InvSig <- matrix(0, nrow=3, ncol = 3)
  
  InvSig[1,1] = Sig[3,3]*Sig[2,2] - Sig[2,3]*Sig[2,3]
  InvSig[1,2] = Sig[1,3]*Sig[2,3] - Sig[3,3]*Sig[1,2]
  InvSig[1,3] = Sig[1,2]*Sig[2,3] - Sig[1,3]*Sig[2,2]
  
  InvSig[2,2] = Sig[3,3]*Sig[1,1] - Sig[1,3]*Sig[1,3]
  InvSig[2,3] = Sig[1,2]*Sig[1,3] - Sig[1,1]*Sig[2,3] ## 
  
  InvSig[3,3] = Sig[1,1]*Sig[2,2] - Sig[1,2]*Sig[1,2]
  
  
  DetVal <- Sig[1,1]*InvSig[1,1] + Sig[1,2]*InvSig[1,2] + Sig[1,3]*InvSig[1,3]
  DetVal
  
  InvSigF <- (1/DetVal)*(InvSig + t(InvSig) - diag(diag(InvSig)) )
  list(InverseCov = InvSigF, Partial_Correlation = (-1)*cov2cor(InvSigF))
}



names <-c("N","betx","betz","lbz","tau",
          "corSel",
          "corFull",
          "PcorSel",
          "PcorFull",
          
          "mb_1","mb_2",
          "Fmb_1","Fmb_2",
          "mab_1","mab_2", "mab_3",
          "Fmab_1","Fmab_2","Fmab_3",
          #texting
          "mab_intp","mab_alx","mab_alz","mab_pb1","mab_pb2","mab_pb3", 
          "mb_intp","mb_alx","mb_pb1","mb_pb2",
          "Fmab_intp","Fmab_alx","Fmab_alz","Fmab_pb1","Fmab_pb2","Fmab_pb3", 
          "Fmb_intp","Fmb_alx","Fmb_pb1","Fmb_pb2")


SimMisconcept12 <- function(N){
  
  bet.z <- .5     #runif(1,0,1)
  bet.x <- .45      #runif(1,0,1)
  lb.z <- .4
  al.x <- 0.     #runif(1,0,1)
  al.z <- .6     #runif(1,0,1)
  
  #--- simulate
  Z <- rnorm(n=N)
  #--- Simulate X 
  #varx.z <- 1 - lb.z^2
  varx.z <- 1 - var(Z)*lb.z^2
  X = lb.z*Z + rnorm(n=N, mean = 0, sd = sqrt(varx.z))
  #X <- X/sd(X) 
  #--- Simulate Y2  
  #vary2.z = 1 - al.z^2
  vary2.z = 1 - (var(Z)*al.z^2 + (al.x^2)*var(X) + 2*(al.x)*al.z*cov(X, Z))
  Y2 <- al.z*Z + rnorm(n=N, mean = 0, sd = sqrt(vary2.z))
  #Y2 <- Y2/sd(Y2)
  #-- Simulate Y1
  #vary1.xz <- 1 - (bet.z^2 +  (bet.x^2)*(lb.z^2) + 2*bet.z*bet.x*lb.z + (bet.x^2))
  vary1.xz <- 1 - (var(Z)*bet.z^2 +  var(X)*(bet.x^2) + 2*bet.z*bet.x*cov(X, Z))
  Y1 <-  bet.z*Z + bet.x*X + rnorm(n=N, mean = 0, sd= sqrt(vary1.xz))
  #Y1 <- Y1/sd(Y1)
  # CutVec <- seq(.1,2,length.out =300)
  # CoVec <- CorXZ(CutVec, Y1, X, Z)
  # Idp <- which(CoVec > 0)
  # IdN <- which(CoVec < 0)
  # 
  # FInd cut-off pount a for selecting on Y1
  OpRes <- optim(1.6, CorrelValZX, al.z=al.z, bet.z = bet.z, bet.x = bet.x, lb.z = lb.z)
  #taucut0 = 2.366733
  taucut0 = OpRes$par
  id <- which(Y1 > taucut0 )#taucut0*sd(Y1) + mean(Y1))
  #All_Models= list()
  mAB <- (lm(Y2[id] ~ X[id] + Z[id]))
  mB <- (lm(Y2[id] ~ X[id]))
  mA <- (lm(Y2[id] ~ Z[id]))
  
  # Res <- c( c((cbind(summary(mAB)$coefficients[,1], confint.lm(mAB)))), #vector of length 9 (mean - 2.5% - 97.5%)
  #           c((cbind(summary(mA)$coefficients[,1], confint.lm(mA)))), #vector of length 6 (mean - 2.5% - 97.5%)
  #           c((cbind(summary(mB)$coefficients[,1], confint.lm(mB)))), #vector of length 6 (mean - 2.5% - 97.5%)
  #           cor(FatA, Y), # crude correlations
  #           cor(FatB, Y), # crude correlations
  #           pcor(data.frame(Y = Y, FatB = FatB, FatA = FatA))$estimate[1,-1] # partial correlations vector of length 2
  #           )
  Parm = c(0.0, al.x, al.z)
  Res <- c( c(cor(X[id], Z[id]),taucut0, c(summary(mAB)$coefficients[,1], (confint.lm(mAB)[,2] > Parm)*(confint.lm(mAB)[,1] < Parm)) ), #vector of length 6 (mean - cover in the 95%(1) or mot (0) )
            c(summary(mA)$coefficients[,1], (confint.lm(mA)[,2] > Parm[-3])*(confint.lm(mA)[,1] < Parm[-3])), #(mean - cover in the 95%(1) or mot (0) )
            c(summary(mB)$coefficients[,1], (confint.lm(mB)[,2] > Parm[-2])*(confint.lm(mB)[,1] < Parm[-2])), #(mean - cover in the 95%(1) or mot (0) )
            
            cor(X[id], Y2[id]), # crude correlations
            cor(Z[id], Y2[id]), # crude correlations
            pcor(data.frame(Y2 = Y2[id], X = X[id], Z = Z[id]))$estimate[1,-1], # partial correlations vector of length 2
            pcor(data.frame(Y2 = Y2, X = X, Z = Z))$estimate[1,-1] # partial correlations vector of length 2
            #pcor(data.frame(Y = Y, FatA = FatA, FatB = FatB, X = X))$estimate[1,-1]
  )
  
  (Res)
}

#@ CovMat :  Compute the theoritical correlation  matrix (yield a covariance matrix for the vector Z, X, Y2, Y1)
CovMat <- function(al.z, bet.z, bet.x, lb.z, al.x){
  
  Sig <- matrix(0, nrow=4,ncol=4)
  Sig[1,2] = Sig[2,1] =  lb.z  # cor(X,Z)
  Sig[1,3] = Sig[3,1] =  al.z + al.x*lb.z # cor(Y2, Z)
  Sig[1,4] = Sig[4,1] = bet.z + bet.x*lb.z  #cor(Y1, Z)
  
  Sig[2,3] = Sig[3,2] = al.x + al.z*lb.z  #cor(Y2, X)
  Sig[2,4] = Sig[4,2] = bet.x + bet.z*lb.z  #cor(Y1, X)
  
  Sig[4,3] = Sig[3,4] = bet.x*(al.x + al.z*lb.z) + bet.z*(al.z + lb.z*al.x)   #(al.x*lb.z + al.z)*(bet.x*lb.z + bet.z)  #cor(Y2, Y1)
  diag(Sig) <- 1
  
  Sig
  
}

#@ EstimCov :  Compute empirical mean and covariance matrix between Z, X Y2
EstimCov<- function(Z,X,Y2){
  
  meanVal = c(mean(Z), mean(X), mean(Y2))
  
  dat <- data.frame(Z, X, Y2)
  Covm <- cov(dat)
  list(Mean = meanVal, Cov = Covm)
  
}

#@ CorrelValZX4 :  Compute theoritical means and covariance matrix between Z, X Y2, Y1 (given value of a)
CorrelValZX4 <- function(al.z, bet.z, bet.x, lb.z, a){
  
  Sig <- CovMat(al.z, bet.z, bet.x, lb.z, al.x)
  m1 = MTS::msqrt(Sig[1:3,1:3])
  
  
  #---------------
  Al = c(1/sqrt( 1 - (t(Sig[1:3,4])%*%solve(Sig[1:3,1:3])%*%(Sig[1:3,4]) ) ))*solve(Sig[1:3,1:3])%*%Sig[1:3,4]
  Del = (1/sqrt(c(1 + t(Al)%*%Sig[1:3,1:3]%*%Al)))*Sig[1:3,1:3]%*%Al
  
  # (1/sqrt( 1 - (t(Sig[1:3,4])%*%solve(Sig[1:3,1:3])%*%(Sig[1:3,4]) ) ))*(1.6)
  #  1.6*sqrt(c(1 + t(Al)%*%Sig[1:2,1:2]%*%Al))
  
  #delta = c(msqrt(1 + t(Lb1)%*%Sig[1:2,1:2]%*%(Lb1))$invsqrt)*Sig[1:2,1:2]%*%Lb1
  mu = c((Del*(dnorm(-a)/pnorm(-a))), dnorm(a)/pnorm(-a))
  
  # Sig[1:2,1:2] - diag(1/c(mu))%*%(tcrossprod(mu))%*%diag(1/c(mu)) #m1$invsqrt%*%(tcrossprod(mu))%*%m1$invsqrt
  SigF <- matrix(0, nrow=4, ncol=4)
  SigF[1:3,1:3] =  Sig[1:3,1:3] + (-(dnorm(-a)/pnorm(-a))^2 - (-a*dnorm(-a)/pnorm(-a)))*(tcrossprod(Del)) #m1$invsqrt%*%(tcrossprod(mu))%*%m1$invsqrt
  SigF[,4] = SigF[4,] = ( 1 + (-(dnorm(-a)/pnorm(-a))^2 - (-a*dnorm(-a)/pnorm(-a))) )*Sig[,4]
  
  
  Res = list(Mu = mu, Sig = SigF, CorSig = cov2cor(SigF)) 
  Res
  
}

#@ CorrelValObj :  objective dunction to obtain a
CorrelValObj <- function(al.z, bet.z, bet.x, lb.z, a){
  
  Sig <- CovMat(al.z, bet.z, bet.x, lb.z, al.x)
  m1 = MTS::msqrt(Sig[1:3,1:3])
  
  
  #---------------
  Al = c(1/sqrt( 1 - (t(Sig[1:3,4])%*%solve(Sig[1:3,1:3])%*%(Sig[1:3,4]) ) ))*solve(Sig[1:3,1:3])%*%Sig[1:3,4]
  Del = (1/sqrt(c(1 + t(Al)%*%Sig[1:3,1:3]%*%Al)))*Sig[1:3,1:3]%*%Al
  
  # (1/sqrt( 1 - (t(Sig[1:3,4])%*%solve(Sig[1:3,1:3])%*%(Sig[1:3,4]) ) ))*(1.6)
  #  1.6*sqrt(c(1 + t(Al)%*%Sig[1:2,1:2]%*%Al))
  
  #delta = c(msqrt(1 + t(Lb1)%*%Sig[1:2,1:2]%*%(Lb1))$invsqrt)*Sig[1:2,1:2]%*%Lb1
  mu = (Del*(dnorm(-a)/pnorm(-a)))
  
  # Sig[1:2,1:2] - diag(1/c(mu))%*%(tcrossprod(mu))%*%diag(1/c(mu)) #m1$invsqrt%*%(tcrossprod(mu))%*%m1$invsqrt
  SigF =  Sig[1:3,1:3] + (-(dnorm(-a)/pnorm(-a))^2 - (-a*dnorm(-a)/pnorm(-a)))*(tcrossprod(Del)) #m1$invsqrt%*%(tcrossprod(mu))%*%m1$invsqrt
  
  Res = list(Mu = mu, Sig = SigF, CorSig = cov2cor(SigF)) 
  (SigF[1,2])^2
  
}
CorrelValObj(al.z, bet.z, bet.x, lb.z, 1.6)


CorrelValObjV2 <- function(al.z, bet.z, bet.x, lb.z, a){
  
  (lb.z/((bet.z+lb.z*bet.x)*(bet.x+lb.z*bet.z)) + (-(dnorm(-a)/pnorm(-a))^2 + (a*dnorm(-a)/pnorm(-a))))^2
  
}

Obj = optimize(CorrelValObjV2, c(-7, 7), al.z=al.z, bet.z = bet.z, bet.x = bet.x, lb.z = lb.z)
Obj
CorrelValObjV2(al.z, bet.z, bet.x, lb.z, 1.6)
X0 = seq(-3,3,length.out=1000)

plot(X0, CorrelValObjV2(al.z, bet.z, bet.x, lb.z, X0))



optim(1.6, CorrelValObj, al.z=al.z, bet.z = bet.z, bet.x = bet.x, lb.z = lb.z)

Obj = optimize(CorrelValObj, c(-7, 7), al.z=al.z, bet.z = bet.z, bet.x = bet.x, lb.z = lb.z)
Obj

#--- SImulate data
#--- simulate
N = 50000
Z <- rnorm(n=N)
Z <- Z/sd(Z)
#--- Simulate X 
#----------------------------------------------------------------------
varx.z <- 1 - var(Z)*lb.z^2
X = lb.z*Z + rnorm(n=N, mean = 0, sd = sqrt(varx.z))
X <- X/sd(X) 
#--- Simulate Y2  
#----------------------------------------------------------------------
#vary2.z = 1 - (al.z^2 +  (al.x^2)*(lb.z^2) + 2*al.z*al.x*lb.z + (al.x^2))
vary2.z = 1 - (var(Z)*al.z^2 + (al.x^2)*var(X) + 2*(al.x)*al.z*cov(X, Z))
Y2 <- al.z*Z + al.x*X + rnorm(n=N, mean = 0, sd = sqrt(vary2.z))
Y2 <- Y2/sd(Y2)
#-- Simulate Y1
#----------------------------------------------------------------------
#vary1.xz <- 1 - (bet.z^2 +  (bet.x^2)*(lb.z^2) + 2*bet.z*bet.x*lb.z + (bet.x^2))
vary1.xz <- 1 - (var(Z)*bet.z^2 +  var(X)*(bet.x^2) + 2*bet.z*bet.x*cov(X, Z))
Y1 <-  bet.z*Z + bet.x*X + rnorm(n=N, mean = 0, sd = sqrt(vary1.xz))
Y1 <- Y1/sd(Y1)
id <- which(Y1 > Obj$minimum) ## id to be selected
length(id)

SigTe <- CovMat(al.z, bet.z, bet.x, lb.z, al.x) ## Covariance matrix of 4 times 4
SigTe
solve(SigTe[1:3,1:3])
InvMat3by3(SigTe[1:3,1:3])
dat=data.frame(Z=Z, X=X, Y2=Y2, Y1=Y1)
cov(dat)
PartialSlope(SigTe[1:3,1:3])

al.x*bet.x + al.x*lb.z*bet.z+al.z*bet.z + al.z*lb.z*bet.x
bet.x*(al.x + al.z*lb.z) + bet.z*(al.z + lb.z*al.x)


EstimCov(Z, X, Y2)                     ## Estimate covariance based on the data
EstimCov(Z[id], X[id], Y2[id])         ## Estimated Mean and Covriance after selection
dat=data.frame(Z=Z[id], X=X[id], Y2=Y2[id], Y1=Y1[id])
cov(dat)
cov2cor(cov(dat))
SigTea <- CorrelValZX4(al.z, bet.z, bet.x, lb.z, Obj$minimum)
SigTea
PartialSlope(SigTea$Sig[1:3,1:3])

summary(lm(Y2[id] ~  X[id]))
summary(lm(Y2[id] ~  Z[id]))

summary(lm(Y2[id] ~ X[id] + Z[id]))
summary(lm(Y2 ~ X+Z))
summary(lm(Y2 ~ X))

summary(lm(Y1[id] ~  X[id]))
summary(lm(Y1[id] ~  Z[id]))
summary(lm(Y1[id] ~ X[id] + Z[id]))
summary(lm(Y1 ~ X+Z))
summary(lm(Y1 ~ X))



a = Obj$minimum
1 + (-(dnorm(-a)/pnorm(-a))^2 - (-a*dnorm(-a)/pnorm(-a)))
(1 + (-(dnorm(-a)/pnorm(-a))^2 - (-a*dnorm(-a)/pnorm(-a))))*(cov(dat)[,4])

(bet.x + lb.z*bet.z)*(bet.z + bet.x*lb.z)

a <- seq(-10,10, length.out = 5000)
a= 2.36494
V0 <- ((dnorm(-a)/pnorm(-a))^2 + (-a*dnorm(-a)/pnorm(-a)))
plot(a, V0)

lb.z0 <- seq(0,1,length= 1000)
V1 <- (bet.x + lb.z0*bet.z)*(bet.z + bet.x*lb.z0)
plot(lb.z0, V1)
abline(a=0,b=1,lwd=3,col="red")

(bet.x + lb.z*bet.z)*(bet.z + bet.x*lb.z)

cor(X,Z)
cor(Z, Y2)
cor(Z, Y1)

cor(X,Y1)
lb.z*(bet.z + bet.x*lb.z)

cor(X, Y2)
cor(X[id], Y2[id])

cor(Y2, Y1)


Sig <-  CovMat(al.z, bet.z, bet.x, lb.z, al.x)[1:3,1:3] ## Check the partial correlation not conditionning on Y1
Sig <- CorrelValZX4(al.z, bet.z, bet.x, lb.z, Obj$minimum)$CorSig[1:3,1:3] ## Check the partial correlation not conditionning on Y1


#--- Compute the inverse of the Correlation Matrix

InvMat3by3 <- function(Sig){
  InvSig <- matrix(0, nrow=3, ncol = 3)
  
  InvSig[1,1] = Sig[3,3]*Sig[2,2] - Sig[2,3]*Sig[2,3]
  InvSig[1,2] = Sig[1,3]*Sig[2,3] - Sig[3,3]*Sig[1,2]
  InvSig[1,3] = Sig[1,2]*Sig[2,3] - Sig[1,3]*Sig[2,2]
  
  InvSig[2,2] = Sig[3,3]*Sig[1,1] - Sig[1,3]*Sig[1,3]
  InvSig[2,3] = Sig[1,2]*Sig[1,3] - Sig[1,1]*Sig[2,3] ## 
  
  InvSig[3,3] = Sig[1,1]*Sig[2,2] - Sig[1,2]*Sig[1,2]
  
  
  DetVal <- Sig[1,1]*InvSig[1,1] + Sig[1,2]*InvSig[1,2] + Sig[1,3]*InvSig[1,3]
  DetVal
  
  InvSigF <- (1/DetVal)*(InvSig + t(InvSig) - diag(diag(InvSig)) );
  InvSigF
}


