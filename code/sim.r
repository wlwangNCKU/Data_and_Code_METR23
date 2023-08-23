source(paste(PATH, 'function/Bay.mnlmmcm.MCMC.fn.r',sep=""))
source(paste(PATH, 'function/Bay.mnlmm.na.fn.r',sep=""))
source(paste(PATH, 'function/mnlmm.na.fn.r',sep=""))
source(paste(PATH, 'function/mnlmmcm.fn.r',sep=""))

###mse function
mse = function(a, b) mean((a-b)^2,na.rm=T)

summary.Table = function(Loop, Case, vecData, mleCM, bayCM, mleM, bayM)
{
  Burn.in = c(Loop, Case, bayCM$burn.in, bayM$burn.in) 
  CPU = c(Loop, Case, mleCM$run.sec, bayCM$CPUtime, mleCM$run.sec, bayM$CPUtime) 
  
  ##Table : MLE for Missing + Censor
  LCI.cm<-mleCM$SD$SD$out[1,] - 1.96*mleCM$SD$SD$out[2,]
  UCI.cm<-mleCM$SD$SD$out[1,] + 1.96*mleCM$SD$SD$out[2,]
  Len.cm<-UCI.cm-LCI.cm
  mleCP.cm<-as.numeric(LCI.cm<=par.true & UCI.cm>=par.true)
  mle.cm<-t(rbind(Case,par.true,mleCM$SD$SD$out,LCI.cm,UCI.cm,Len.cm,t(mleCP.cm)))
  colnames(mle.cm)<-c("Case","par.true","Est","SE","LCI","UCI","Len","CP")
  
  ##Table : Bayesian for Missing + Censor
  bayCP.cm<-as.numeric(bayCM$pos.inf$theta.out[3,]<=par.true &
                         bayCM$pos.inf$theta.out[4,]>=par.true)
  bayLen.cm<-bayCM$pos.inf$theta.out[4,]-bayCM$pos.inf$theta.out[3,]
  bay.cm <- cbind(Case,par.true,t(bayCM$pos.inf$theta.out[-5,]),bayLen.cm,bayCP.cm)
  
  ##Table : MLE for Missing
  LCI.m<-mleM$SD$out[1,]-1.96*mleM$SD$out[2,]
  UCI.m<-mleM$SD$out[1,]+1.96*mleM$SD$out[2,]
  Len.m<-UCI.m-LCI.m
  mleCP.m<-as.numeric(LCI.m<=par.true & UCI.m>=par.true)
  mle.m<-t(rbind(Case,par.true,mleM$SD$out,LCI.m,UCI.m,Len.m,t(mleCP.m)))
  colnames(mle.m)<-c("Case","par.true","Est","SE","LCI","UCI","Len","CP")
  
  ###Table : Bayesian for Missing
  bay.m <- cbind(Case,par.true,t(bayM$est[-5,]))
  
  return(list(burn.in=Burn.in, CPU=CPU, mle.cm=mle.cm, bay.cm=bay.cm, mle.m=mle.m, bay.m=bay.m))
}

summary.yfit = function(Loop, Case, mleCM, bayCM, mleM, bayM)
{
  ##yfit : Missing + Censor
  mleCM.yfit <- mleCM$yfit
  bayCM.yfit <- bayCM$pos.inf$yfit.out
  
  ###yfit : Missing
  mleM.yfit <- mleM$yfit
  bayM.yfit <- bayM$pos.summ$yfit.out
  
  return(list(mleCM.yfit=mleCM.yfit, bayCM.yfit=bayCM.yfit, mleM.yfit=mleM.yfit, bayM.yfit=bayM.yfit))
}

r = 2; q =2 
Ti = seq(10, 100, 10)
si = length(Ti)
Xi = matrix(c(rep(1, si), 1: si), ncol=2)
Zi = matrix(Xi[, 1], ncol=1)
rXi = kronecker(diag(r), Xi)
rZi = kronecker(diag(r), Zi)
beta = c(5, 17, 7, 2, 0.05)
DD = matrix(c(1, 0.25, 0.25, 1), q)
Sigma = diag(r)
rho = 0.75
Sigma[1, 2] = Sigma[2, 1] = rho * sqrt(prod(diag(Sigma[1:2, 1:2])))
phi = 1e-6 ; ga = 1  #UNC 
par = list(beta=beta, DD=DD, Sigma=Sigma, phi=phi, ga=ga)
par.true = c(beta, DD[vech.posi(q)], Sigma[vech.posi(r)], phi, ga)
init.par = list(beta=beta+0.1*runif(1), DD=DD, Sigma=Sigma, phi=1e-6, ga=1)

Loop <- 1
while(Loop <= 100)
{
  cat(paste(c(rep('=', 25), rep(' ', 5), 'The ', Loop, 'th replication: N = ', N, rep(' ', 5), rep('=', 25)), sep = '', collapse = ''), '\n')
  gen.data = rmnlmm(Ti, rXi, rZi, beta, DD, Sigma, phi, ga, N, r)                # complete data
  
  ############################ Case 1. MAR + Censor 10% #############################
  data.na = rmnlmm.na(gen.data$Ymat, type='MAR')
  sim.data = rmnlmm.cen(data.na$Ymat, data.na$Ymat.na, brna=gen.data$brna, cen.rate=0.1)
  Data = sim.data$Data
  vecData = sim.data$vecData
  
  ### censor + missing (MAR + Censor 10%)  
  cat(paste(c(rep('=', 25), rep(' ', 5), 'Repication = ', Loop, '; MAR + Censor 10%', rep(' ', 5), rep('=', 25)), sep = '', collapse = ''), '\n')
  mleCM1 = MNLMMcm.ECM(init.par, tol=1e-3, typeCi='UNC', max.iter = 500, per=100)
  hyperCM1 = hyper.par(mleCM1$para, mleCM1$SD$FI)
  repeat{
    bayCM1 = try(MCMC.MNLMMCM(chain=2, ITER=2500, hyperCM1, typeCi='UNC', nlag=10, per=100, method='pseudo'),silent = T)
    if(class(bayCM1) != "try-error") break
  }
  
  ### only missing (MAR)
  cat(paste(c(rep('=', 25), rep(' ', 5), 'Repication = ', Loop, '; MAR', rep(' ', 5), rep('=', 25)), sep = '', collapse = ''), '\n')
  mle1 = MNLMMna.ECM(init.par, tol=1e-3, typeCi='UNC', max.iter = 500, per=100)
  hyper1 = hyper.par(mle1$para, mle1$SD)
  repeat{
    bay1 = try(MCMC.MNLMMna(chain=2, ITER=2500, hyper1, typeCi='UNC', nlag=10, per=100),silent = T)
    if(class(bay1) != "try-error") break
  }
  
  sum.tab1= summary.Table(Loop=Loop, Case=1, vecData, mleCM1, bayCM1, mle1, bay1)
  yfit.c1 = summary.yfit(Loop=Loop, Case=1, mleCM1, bayCM1, mle1, bay1)
  y.c1 = vecData$y
  mse.c1 = c(mse(mleCM1$yfit, vecData$y), mse(bayCM1$pos.inf$yfit.out[1,], vecData$y), mse(mle1$yfit, vecData$y), mse(bay1$pos.summ$yfit.out[1,], vecData$y))
  
  ########################## Case 2. MCAR 10% + Censor 10% ##########################
  data.na = rmnlmm.na(gen.data$Ymat, type='MCAR', na.rate=0.1)
  sim.data = rmnlmm.cen(data.na$Ymat, data.na$Ymat.na, brna=gen.data$brna, cen.rate=0.1)
  Data = sim.data$Data
  vecData = sim.data$vecDat
  
  ### censor + missing (MCAR + Censor 10%)  
  cat(paste(c(rep('=', 25), rep(' ', 5), 'Repication = ', Loop, '; MCAR 10% + Censor 10%', rep(' ', 5), rep('=', 25)), sep = '', collapse = ''), '\n')
  mleCM2 = MNLMMcm.ECM(init.par, tol=1e-3, typeCi='UNC', max.iter = 500, per=100)
  hyperCM2 = hyper.par(mleCM2$para, mleCM2$SD$FI)
  repeat{
    bayCM2 = try(MCMC.MNLMMCM(chain=2, ITER=2500, hyperCM2, typeCi='UNC', nlag=10, per=100, method='pseudo'),silent = T)
    if(class(bayCM2) != "try-error") break
  }
  
  ### only missing (MCAR 10%)
  cat(paste(c(rep('=', 25), rep(' ', 5), 'Repication = ', Loop, '; MCAR 10%', rep(' ', 5), rep('=', 25)), sep = '', collapse = ''), '\n')
  mle2 = MNLMMna.ECM(init.par, tol=1e-3, typeCi='UNC', max.iter = 500, per=100)
  hyper2 = hyper.par(mle2$para, mle2$SD)
  repeat{
    bay2 = try(MCMC.MNLMMna(chain=2, ITER=2500, hyper2, typeCi='UNC', nlag=10, per=100),silent = T)
    if(class(bay2) != "try-error") break
  }
  
  sum.tab2 = summary.Table(Loop=Loop, Case=2, vecData, mleCM2, bayCM2, mle2, bay2)
  yfit.c2 = summary.yfit(Loop=Loop, Case=2, mleCM2, bayCM2, mle2, bay2)
  y.c2 = vecData$y
  mse.c2 = c(mse(mleCM2$yfit, vecData$y), mse(bayCM2$pos.inf$yfit.out[1,], vecData$y), mse(mle2$yfit, vecData$y), mse(bay2$pos.summ$yfit.out[1,], vecData$y))
  
  ############################# Case 3. MAR + Censor 30% ############################
  data.na = rmnlmm.na(gen.data$Ymat, type='MAR') 
  sim.data = rmnlmm.cen(data.na$Ymat, data.na$Ymat.na, brna=gen.data$brna, cen.rate=0.3)
  Data = sim.data$Data
  vecData = sim.data$vecDat
  
  ### censor + missing (MAR + Censor 30%)  
  cat(paste(c(rep('=', 25), rep(' ', 5), 'Repication = ', Loop, '; MCAR 10% + Censor 30%', rep(' ', 5), rep('=', 25)), sep = '', collapse = ''), '\n')
  mleCM3 = MNLMMcm.ECM(init.par, tol=1e-3, typeCi='UNC', max.iter = 500, per=100)
  hyperCM3 = hyper.par(mleCM3$para, mleCM3$SD$FI)
  repeat{
    bayCM3 = try(MCMC.MNLMMCM(chain=2, ITER=2500, hyperCM3, typeCi='UNC', nlag=10, per=100, method='pseudo'),silent = T)
    if(class(bayCM3) != "try-error") break
  }
  
  ### only missing (MAR 10%)
  cat(paste(c(rep('=', 25), rep(' ', 5), 'Repication = ', Loop, '; MCAR 10%', rep(' ', 5), rep('=', 25)), sep = '', collapse = ''), '\n')
  mle3 = MNLMMna.ECM(init.par, tol=1e-3, typeCi='UNC', max.iter = 500, per=100)
  hyper3 = hyper.par(mle3$para, mle3$SD)
  repeat{
    bay3 = try(MCMC.MNLMMna(chain=2, ITER=2500, hyper3, typeCi='UNC', nlag=10, per=100),silent = T)
    if(class(bay3) != "try-error") break
  }
  
  sum.tab3 = summary.Table(Loop=Loop, Case=3, vecData, mleCM3, bayCM3, mle3, bay3)
  yfit.c3 = summary.yfit(Loop=Loop, Case=3, mleCM3, bayCM3, mle3, bay3)
  y.c3 = vecData$y
  mse.c3 = c(mse(mleCM3$yfit, vecData$y), mse(bayCM3$pos.inf$yfit.out[1,], vecData$y), mse(mle3$yfit, vecData$y), mse(bay3$pos.summ$yfit.out[1,], vecData$y))
  
  ########################## Case 4. MCAR 10% + Censor 30% ##########################
  data.na = rmnlmm.na(gen.data$Ymat, type='MCAR', na.rate=0.1)
  sim.data = rmnlmm.cen(data.na$Ymat, data.na$Ymat.na, brna=gen.data$brna, cen.rate=0.3)
  Data = sim.data$Data
  vecData = sim.data$vecDat
  
  ### censor + missing (MCAR + Censor 30%)  
  cat(paste(c(rep('=', 25), rep(' ', 5), 'Repication = ', Loop, '; MCAR 10% + Censor 10%', rep(' ', 5), rep('=', 25)), sep = '', collapse = ''), '\n')
  mleCM4 = MNLMMcm.ECM(init.par, tol=1e-3, typeCi='UNC', max.iter = 500, per=100)
  hyperCM4 = hyper.par(mleCM4$para, mleCM4$SD$FI)
  repeat{
    bayCM4 = try(MCMC.MNLMMCM(chain=2, ITER=2500, hyperCM4, typeCi='UNC', nlag=10, per=100, method='pseudo'),silent = T)
    if(class(bayCM4) != "try-error") break
  }
  
  ### only missing (MCAR 10%)
  cat(paste(c(rep('=', 25), rep(' ', 5), 'Repication = ', Loop, '; MCAR 10%', rep(' ', 5), rep('=', 25)), sep = '', collapse = ''), '\n')
  mle4 = MNLMMna.ECM(init.par, tol=1e-3, typeCi='UNC', max.iter = 500, per=100)
  hyper4 = hyper.par(mle4$para, mle4$SD)
  repeat{
    bay4 = try(MCMC.MNLMMna(chain=2, ITER=2500, hyper4, typeCi='UNC', nlag=10, per=100),silent = T)
    if(class(bay4) != "try-error") break
  }
  
  sum.tab4 = summary.Table(Loop=Loop, Case=4, vecData, mleCM4, bayCM4, mle4, bay4)
  yfit.c4 = summary.yfit(Loop=Loop, Case=4, mleCM4, bayCM4, mle4, bay4)
  y.c4 = vecData$y
  mse.c4 = c(mse(mleCM4$yfit, vecData$y), mse(bayCM4$pos.inf$yfit.out[1,], vecData$y), mse(mle4$yfit, vecData$y), mse(bay4$pos.summ$yfit.out[1,], vecData$y))
  
  ############################ Case 5. MAR + Censor 50% #############################
  data.na = rmnlmm.na(gen.data$Ymat, type='MAR')
  sim.data = rmnlmm.cen(data.na$Ymat, data.na$Ymat.na, brna=gen.data$brna, cen.rate=0.5)
  Data = sim.data$Data
  vecData = sim.data$vecData
  
  ### censor + missing (MAR + Censor 50%)  
  cat(paste(c(rep('=', 25), rep(' ', 5), 'Repication = ', Loop, '; MAR + Censor 50%', rep(' ', 5), rep('=', 25)), sep = '', collapse = ''), '\n')
  mleCM5 = MNLMMcm.ECM(init.par, tol=1e-3, typeCi='UNC', max.iter = 500, per=100)
  hyperCM5 = hyper.par(mleCM5$para, mleCM5$SD$FI)
  repeat{
    bayCM5 = try(MCMC.MNLMMCM(chain=2, ITER=2500, hyperCM5, typeCi='UNC', nlag=10, per=100, method='pseudo'),silent = T)
    if(class(bayCM5) != "try-error") break
  }
  
  ### only missing (MAR)
  cat(paste(c(rep('=', 25), rep(' ', 5), 'Repication = ', Loop, '; MAR', rep(' ', 5), rep('=', 25)), sep = '', collapse = ''), '\n')
  mle5 = MNLMMna.ECM(init.par, tol=1e-3, typeCi='UNC', max.iter = 500, per=100)
  hyper5 = hyper.par(mle5$para, mle5$SD)
  repeat{
    bay5 = try(MCMC.MNLMMna(chain=2, ITER=2500, hyper5, typeCi='UNC', nlag=10, per=100),silent = T)
    if(class(bay5) != "try-error") break
  }
  
  sum.tab5= summary.Table(Loop=Loop, Case=5, vecData, mleCM5, bayCM5, mle5, bay5)
  yfit.c5 = summary.yfit(Loop=Loop, Case=5, mleCM5, bayCM5, mle5, bay5)
  y.c5 = vecData$y
  mse.c5 = c(mse(mleCM5$yfit, vecData$y), mse(bayCM5$pos.inf$yfit.out[1,], vecData$y), mse(mle5$yfit, vecData$y), mse(bay5$pos.summ$yfit.out[1,], vecData$y))
  
  ########################## Case 6. MCAR 10% + Censor 50% ##########################
  data.na = rmnlmm.na(gen.data$Ymat, type='MCAR', na.rate=0.1)
  sim.data = rmnlmm.cen(data.na$Ymat, data.na$Ymat.na, brna=gen.data$brna, cen.rate=0.5)
  Data = sim.data$Data
  vecData = sim.data$vecDat
  
  ### censor + missing (MCAR + Censor 50%)  
  cat(paste(c(rep('=', 25), rep(' ', 5), 'Repication = ', Loop, '; MCAR 10% + Censor 50%', rep(' ', 5), rep('=', 25)), sep = '', collapse = ''), '\n')
  mleCM6 = MNLMMcm.ECM(init.par, tol=1e-3, typeCi='UNC', max.iter = 500, per=100)
  hyperCM6 = hyper.par(mleCM6$para, mleCM6$SD$FI)
  repeat{
    bayCM6 = try(MCMC.MNLMMCM(chain=2, ITER=2500, hyperCM6, typeCi='UNC', nlag=10, per=100, method='pseudo'),silent = T)
    if(class(bayCM6) != "try-error") break
  }
  
  ### only missing (MCAR 10%)
  cat(paste(c(rep('=', 25), rep(' ', 5), 'Repication = ', Loop, '; MCAR 10%', rep(' ', 5), rep('=', 25)), sep = '', collapse = ''), '\n')
  mle6 = MNLMMna.ECM(init.par, tol=1e-3, typeCi='UNC', max.iter = 500, per=100)
  hyper6 = hyper.par(mle6$para, mle6$SD)
  repeat{
    bay6 = try(MCMC.MNLMMna(chain=2, ITER=2500, hyper6, typeCi='UNC', nlag=10, per=100),silent = T)
    if(class(bay6) != "try-error") break
  }
  
  sum.tab6 = summary.Table(Loop=Loop, Case=6, vecData, mleCM6, bayCM6, mle6, bay6)
  yfit.c6 = summary.yfit(Loop=Loop, Case=6, mleCM6, bayCM6, mle6, bay6)
  y.c6 = vecData$y
  mse.c6 = c(mse(mleCM6$yfit, vecData$y), mse(bayCM6$pos.inf$yfit.out[1,], vecData$y), mse(mle6$yfit, vecData$y), mse(bay6$pos.summ$yfit.out[1,], vecData$y))

  ############################### Print Outputs ###############################
  ### mle.cm ###
  mleCM.table = rbind(sum.tab1$mle.cm, sum.tab2$mle.cm, sum.tab3$mle.cm, sum.tab4$mle.cm, sum.tab5$mle.cm, sum.tab6$mle.cm)
  
  ### bay.cm ###
  bayCM.table = rbind(sum.tab1$bay.cm, sum.tab2$bay.cm, sum.tab3$bay.cm, sum.tab4$bay.cm, sum.tab5$bay.cm, sum.tab6$bay.cm)
  
  ### mle.m ###
  mleM.table = rbind(sum.tab1$mle.m, sum.tab2$mle.m, sum.tab3$mle.m, sum.tab4$mle.m, sum.tab5$mle.m, sum.tab6$mle.m)
  
  ### bay.m ###
  bayM.table = rbind(sum.tab1$bay.m, sum.tab2$bay.m, sum.tab3$bay.m, sum.tab4$bay.m, sum.tab5$bay.m, sum.tab6$bay.m)
  
  ### mse.yfit ###
  mse.yfit = data.frame(Loop, c('MAR.C10','MCAR10.C10','MAR.C30','MCAR10.C30','MAR.C50','MCAR10.C50'), 
                   rbind(mse.c1, mse.c2, mse.c3, mse.c4, mse.c5, mse.c6))
  
  ###Export files
  write.table(mse.yfit, paste(PATH1, 'mseyfit.txt',sep=""), append=T, row.names = F, col.names = F)
  write.table(data.frame(Loop, c("b1","b2","b3","b4","b5","d11","d21","d22","sigma11","sigma21","sigma22","phi","gamma"), 
                    mleCM.table), paste(PATH1, 'table(mleCM).txt',sep=""), append=T, row.names = F, col.names = F)
  write.table(data.frame(Loop, c("b1","b2","b3","b4","b5","d11","d21","d22","sigma11","sigma21","sigma22","phi","gamma"), 
                    bayCM.table), paste(PATH1, 'table(bayCM).txt',sep=""), append=T, row.names = F, col.names = F)
  write.table(data.frame(Loop, c("b1","b2","b3","b4","b5","d11","d21","d22","sigma11","sigma21","sigma22","phi","gamma"), 
                    mleM.table), paste(PATH1, 'table(mleM).txt',sep=""), append=T, row.names = F, col.names = F)
  write.table(data.frame(Loop, c("b1","b2","b3","b4","b5","d11","d21","d22","sigma11","sigma21","sigma22","phi","gamma"), 
                    bayM.table), paste(PATH1, 'table(bayM).txt',sep=""), append=T, row.names = F, col.names = F)
  
  Loop=Loop+1
}
