# Perform ML and Bayesian Modeling for the A5055 Data
library(nlme)
source(paste(PATH, 'function/mnlmmcm.fn2.r',sep=""))
source(paste(PATH, 'function/mnlmm.na.fn2.r',sep=""))
source(paste(PATH, 'function/Bay.mnlmmcm.MCMC.fn2.r',sep=""))
source(paste(PATH, 'function/Bay.mnlmm.na.fn2.r',sep=""))

# Read Data: A5055
HIV = read.table(paste(PATH, 'data/A5055data.txt', sep=''), header=T)
r = 2
cen.lrna = as.numeric(HIV$rna <= 50)
lrna = HIV$logrna
lrna[which(HIV$logrna <= log10(50))] = log10(50)

reHIV = data.frame(cbind(HIV$Subject, HIV$day, HIV$day/7, HIV$logrna, lrna, HIV$rna, HIV$cd4/HIV$cd8, log10(HIV$cd4), HIV$arm-1, cen.lrna))
colnames(reHIV) = c('Subject','Day','Time','lgcopy','lrna','rna','cd4.cd8','lcd4','arm','cen.lrna')
Data = groupedData(lrna + cd4.cd8 ~ Time|Subject, data=reHIV)
cluster = as.numeric(levels(Data$Subject))[Data$Subject]
N = length(unique(cluster))
Ni = matrix(0,N,1)
for (j in 1:N) Ni[j]=sum(cluster==j)

Subject = Time = Resp = Resp.raw = Censor = Var = NULL
for(i in 1: N){
 Subject = c(Subject, rep(reHIV[which(reHIV$Subject == i),]$Subject, r))
 Time = c(Time, rep(reHIV[which(reHIV$Subject == i),]$Time, r))
 Resp = c(Resp, c(reHIV[which(reHIV$Subject == i),]$lrna, reHIV[which(reHIV$Subject == i),]$cd4.cd8))
 Resp.raw = c(Resp.raw, c(Data[which(Data$Subject == i),]$lgcopy, Data[which(Data$Subject == i),]$cd4.cd8))
 Censor = c(Censor, c(reHIV[which(reHIV$Subject == i),]$cen.lrna, rep(0, Ni[i])))
 Var = c(Var, rep(1:r, each=Ni[i]))
}
Censor[is.na(Resp)] = NA
vec.reHIV = data.frame(cbind(Subject, Time, Var, Resp, Resp.raw, Censor, Mis=as.numeric(is.na(Censor))))
vecData = groupedData(Resp~Time|Subject, data=vec.reHIV)

# Initial values:
beta = c(6, 2, -0.05, 0.1, 14, -2, -1.7, 1.2, 0.3, -0.01, 0.001, 0.001)
Sigma = var(cbind(Data$lrna, Data$cd4.cd8), na.rm=T)
DD = 0.5*diag(6)
phi = 0.5
ga = 1
init.par = list(beta = beta, DD=DD, Sigma=Sigma, phi=phi, ga=ga)

### MNLMM ### 
# ML model fitting 
M0 = MNLMMna.ECM(init.par, tol=1e-5, typeCi='UNC', max.iter = 2000, per=100)
M1 = MNLMMna.ECM(init.par, tol=1e-5, typeCi='AR1', max.iter = 2000, per=100)
M2 = MNLMMna.ECM(init.par, tol=1e-5, typeCi='DEC', max.iter = 2000, per=100)

# Bayesian Model fitting
hyM0 = hyper.par(M0$para, M0$SD$FI, M0$b)
repeat{ 
IBM0 = try(IBF.MLMMna(chain=5, ITER=10000, hyper=hyM0, typeCi='UNC', nlag=10, per=100), silent=F)
          if(class(IBM0) != "try-error") break;
}

hyM1 = hyper.par(M1$para, M1$SD$FI, M1$b)
repeat{
IBM1 = try(IBF.MLMMna(chain=5, ITER=10000, hyper=hyM1, typeCi='AR1', nlag=10, per=100), silent=F)
          if(class(IBM1) != "try-error") break;
}

hyM2 = hyper.par(M2$para, M2$SD$FI, M2$b)
repeat{
IBM2 = try(IBF.MLMMna(chain=5, ITER=10000, hyper=hyM2, typeCi='DEC', nlag=10, per=100), silent=F)
          if(class(IBM2) != "try-error") break;
}

### MNLMM-CM ###
repeat{
CM0 = try(MNLMMcm.ECM(init.par, tol=1e-5, typeCi='UNC', max.iter = 2000, per=100), silent=F)
          if(class(CM0) != "try-error") break;
}
     
init.par = list(beta = CM0$para$beta, DD=CM0$para$DD, Sigma=CM0$para$Sigma, phi=0.75, ga=ga)
repeat{
CM1 = try(MNLMMcm.ECM(init.par, tol=1e-5, typeCi='AR1', max.iter = 2000, per=100), silent=F)
          if(class(CM1) != "try-error") break;
}
repeat{     
CM2 = try(MNLMMcm.ECM(init.par, tol=1e-5, typeCi='DEC', max.iter = 2000, per=100), silent=F)
          if(class(CM2) != "try-error") break;
}

# Bayesian Model fitting
hyCM0 = hyper.par(CM0$para, CM0$SD$FI, CM0$para$b)
repeat{
IBCM0 = try(IBF.MNLMMCM(chain=5, ITER=10000, hyper=hyCM0, typeCi='UNC', nlag=10, per=100, method='pseudo'), silent=F)
          if(class(IBCM0) != "try-error") break;
}

hyCM1 = hyper.par(CM1$para, CM1$SD$FI, CM1$para$b)
repeat{
IBCM1 = try(IBF.MNLMMCM(chain=5, ITER=10000, hyper=hyCM1, typeCi='AR1', nlag=10, per=100, method='pseudo'), silent=F)
          if(class(IBCM1) != "try-error") break;
}

hyCM2 = hyper.par(CM2$para, CM2$SD$FI, CM2$para$b)
repeat{
IBCM2 = try(IBF.MNLMMCM(chain=5, ITER=10000, hyper=hyCM2, typeCi='DEC', nlag=10, per=100, method='pseudo'), silent=F)
          if(class(IBCM2) != "try-error") break;
}
PATH = paste(getwd(),"/Data_and_Code_METR23/",sep="")
save.image(paste(PATH, 'data/A5055FitResult.RData',sep=""))
