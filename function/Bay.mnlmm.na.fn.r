library(nlme)
library(mvtnorm)
library(MCMCpack)

# Approximate observed log-likelihood
Tay.aloglik.na = function(Ytilde, Xtilde, TrZ, beta, DD, Sigma, phi, ga, cumsum.oi, cumsum.ni)
{
  N = length(unique(Data$Subject))
  Y.na = vecData$Resp
  na.ind = which(is.na(Y.na))
  n = length(Y.na)
  q = ncol(DD)
  Ytilde.o = Ytilde[-na.ind]
  Xtilde.o = Xtilde[-na.ind, ]
  no = length(Ytilde.o)
  Y.cent = Ytilde.o - Xtilde.o %*% beta
  TLambda = matrix(0, ncol=n, nrow=n)
  for(i in 1: N)
  {
     if(i == 1) idx = 1: cumsum.oi[i]
     else idx = (cumsum.oi[i-1]+1): cumsum.oi[i]
     rZi = TrZ[idx, ((i-1)*q+1): (i*q)]
     TLambda[idx, idx] = rZi %*% DD %*% t(rZi)+ kronecker(Sigma, DEC(phi, ga, Data$Time[Data$Subject == i]))
  }
  TLambda.oo = TLambda[-na.ind, -na.ind]
  sum.delta.i = log.det.Lam.inv = 0
  for(i in 1: N){
    if(i == 1) idx.o = 1: cumsum.ni[i]
    else idx.o = (cumsum.ni[i-1]+1): cumsum.ni[i]
    Lam.oo.inv = solve(TLambda.oo[idx.o, idx.o])
    log.det.Lam.inv = log.det.Lam.inv + log(det(Lam.oo.inv))
    sum.delta.i = sum.delta.i + t(Y.cent[idx.o]) %*% Lam.oo.inv %*% Y.cent[idx.o]
  }
  n.loglik= -0.5*log(2*pi)*no + 0.5*log.det.Lam.inv - 0.5*sum.delta.i
  return(n.loglik)
}

# M-H algorithm for phi and ga
log.pos.Phi = function(Phi, Y.cent, typeCi, cumsum.oi, cumsum.ni, r)
{
   N=length(unique(Data$Subject))
   r = nrow(Sigma)
   N=length(unique(Data$Subject))
   Y.na = vecData$Resp
   na.ind = which(is.na(Y.na))
   n = length(Y.na)
   Yo.cent = Y.cent[-na.ind]
   TSC.inv = matrix(0, n, n)
   TSC.inv[1: cumsum.oi[1], 1: cumsum.oi[1]] = solve(kronecker(Sigma, DEC(phi, ga, Data$Time[Data$Subject == 1])))
   for(i in 2: N) TSC.inv[(cumsum.oi[i-1]+1): cumsum.oi[i], (cumsum.oi[i-1]+1): cumsum.oi[i]] = solve(kronecker(Sigma, DEC(phi, ga, Data$Time[Data$Subject == i])))
   TSC.oo.inv = TSC.inv[-na.ind, -na.ind]
   if(typeCi == 'AR1'){ phi = Phi; ga = 1
   } else{
    phi = Phi[1]; ga = Phi[2]
   }
   log.det.Ri = log(det(as.matrix(TSC.oo.inv[1: cumsum.ni[1], 1: cumsum.ni[1]])))
   sum.g = t(Yo.cent[1:cumsum.ni[1]]) %*% TSC.oo.inv[1: cumsum.ni[1], 1: cumsum.ni[1]] %*% Yo.cent[1:cumsum.ni[1]]
   for(i in 2: N){
     log.det.Ri = log.det.Ri + log(det(as.matrix(TSC.oo.inv[(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]])))
     sum.g = sum.g + t(Yo.cent[(cumsum.ni[i-1]+1): cumsum.ni[i]]) %*% TSC.oo.inv[(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]] %*% Yo.cent[(cumsum.ni[i-1]+1): cumsum.ni[i]]
   }
   log.pos.phi = 0.5 * (log.det.Ri - sum.g) - 2*log(1+ga)
   return(log.pos.phi)
}

MH.Phi = function(phi, ga, Y.cent, typeCi, cumsum.oi, cumsum.ni, r, Phi.mle)
{
   if(typeCi == 'AR1'){
     Phi.old = phi
     Phi.star.old = log(Phi.old/(1-Phi.old))
     sd.phi.star = 0.01*sqrt(2.4*(1/(Phi.mle[1]*(1-Phi.mle[1])))^2)
#     sd.phi.star = 0.1 * (2.4/sqrt(2))^2 * sd.Phi
     Phi.star.new = rnorm(1, Phi.star.old, sd.phi.star)
     Phi.new = exp(Phi.star.new) / (1 + exp(Phi.star.new))
   }
   if(typeCi == 'DEC'){
     Phi.old = c(phi, ga)
     Phi.star.old = c(log(Phi.old[1]/(1-Phi.old[1])), log(Phi.old[2]))
#     cov.Phi.star = ((2.4/sqrt(2))^2 * diag(c(1/(Phi.mle[1]*(1-Phi.mle[1])), 1/Phi.mle[2]))) %*% sd.Phi %*% diag(c(1/(Phi.mle[1]*(1-Phi.mle[1])), 1/Phi.mle[2]))
#     cov.Phi.star = 0.1*(2.4/sqrt(2))^2 * sd.Phi
     cov.Phi.star = matrix(c((1/(Phi.mle[1]*(1-Phi.mle[1])))^2, 0, 0, 1/Phi.mle[2]^2), 2,2)
     repeat{
     Phi.star.new = rmvnorm(1, Phi.star.old, cov.Phi.star)
     if(Phi.star.new[2]>-2 & Phi.star.new[2]< 0.7) break
     }
     Phi.new = c(exp(Phi.star.new[1]) / (1 + exp(Phi.star.new[1])), exp(Phi.star.new[2]))
   }
   log.pos.old = log.pos.Phi(Phi.old, Y.cent, typeCi, cumsum.oi, cumsum.ni, r)
   log.pos.new = log.pos.Phi(Phi.new, Y.cent, typeCi, cumsum.oi, cumsum.ni, r)
   log.accept = log.pos.new - log.pos.old
   U = runif(1)
   if(log(U) > log.accept) Phi.new = Phi.old
   return(Phi.new)
}

# MCMC
MCMC.MNLMMna = function(chain, ITER, hyper, typeCi=c('UNC','AR1','DEC'), nlag=10, per=100)
{
   cat(paste(c(rep('-', 25), rep(' ', 5), 'Running ... Bayesian MCMC for MNLMM', rep(' ', 5), rep('-', 30)), sep = '', collapse = ''), '\n')
   begin=proc.time()[1]
   typeCi = typeCi[1]
# data
   N = length(unique(vecData$Subject))
   Y.na = vecData$Resp
   na.ind = which(is.na(Y.na))
   n = length(Y.na)
   p = length(init.par$beta); q = nrow(init.par$DD); r = nrow(init.par$Sigma)
   Y = Y.na
   Y[na.ind] = 9999
   A = diag(p)
   B = matrix(c(1,rep(0, 7), 1, 0), ncol=q)
   TrZ = matrix(0, ncol=N*q, nrow=n)
   SigCi = SigCi.inv = as.list(N)
   TLam = matrix(0, nrow=n, ncol=n)
   ti = matrix(NA, 2, N)
   for(i in 1: N){
    ti[1, i] = length(na.omit(Data$ycm1[Data$Subject == i]))
    ti[2, i] = length(na.omit(Data$ycm2[Data$Subject == i]))
   }
   ni = colSums(ti)
   cumsum.ni = cumsum(ni)
   si = numeric(N)
   for(i in 1: N) si[i] = length(Data$Ti[Data$Subject == i])
   oi = si*r
   cumsum.oi = cumsum(oi)
   na.num = c(0,cumsum((oi-ni)[which((oi-ni)!=0)]))
   obs.sub = which(c(oi == colSums(ti)))
   mis.sub = which(c(oi != colSums(ti)))
   num.obs.sub = length(obs.sub)
   num.mis.sub = length(mis.sub)
   mis.ind = Yim = as.list(num.mis.sub)
   for(i in 1: num.mis.sub){
      if(mis.sub[i] == 1) mis.ind[[i]] = is.na(Y.na[1: cumsum.oi[mis.sub[i]]])
      else mis.ind[[i]] = is.na(Y.na[(cumsum.oi[mis.sub[i]-1]+1): cumsum.oi[mis.sub[i]]])
   }
   TSoo = matrix(0, n, n)
# hyper parameters:
   beta0 = hyper$beta0
   J0 = hyper$J0
   J0.inv = solve(hyper$J0)
   d0 = hyper$d0
   G0 = hyper$G0
   s0 = hyper$s0
   H0 = hyper$H0
   Phi.mle = hyper$pg.mle
   m1 = q*(q+1)/2; m2 = r*(r+1)/2
   m = p + m1 + m2 + 2
   Theta = array(NA, dim=c(ITER, (m+1), chain))
   TB = array(NA, dim=c(ITER, N*q, chain))
   TYimp = array(NA, dim=c(ITER, length(na.ind), chain))
   TYfit = array(NA, dim=c(ITER, n, chain))
   ytil = as.list(N)
   loglik.iter = matrix(NA, nrow=chain, ncol=ITER)
   acc.phi = acc.ga = matrix(NA, nrow=chain, ncol=ITER)
   acc.b = array(NA, dim=c(N, ITER, chain))
   vech.DD = vech.posi((q))
   vech.Sig = vech.posi(r)
   Init = NULL

   for(k in 1: chain){
   cat(paste(c(rep('=', 25), rep(' ', 5), 'N = ', N, ' : The ', k,'th chain', rep(' ', 5), rep('=', 30)), sep = '', collapse = ''), '\n')
# random initial values:
   beta = c(rmvnorm(1, beta0, J0))
   DD = riwish(d0, G0)
   Sigma = riwish(s0, H0) 
   phi = runif(1, 0, 1)
   ga = runif(1, 0, 1)
   b = rmvnorm(N, rep(0, q), DD)
   if(r==1 & q==1){ Init = rbind(Init, c(beta, DD, Sigma, phi, ga))
   } else Init = rbind(Init, c(beta, DD[vech.DD], Sigma[vech.Sig], phi, ga))
    
   MU = Xtilde = NULL
   for(i in 1: N){
     eta.i = A %*% beta + B %*% b[i, ]
     mu.i = mu.fn(eta.i, Data$Time[Data$Subject == i])
     MU = c(MU, mu.i$mu1, mu.i$mu2)
     dmu.i = dmu(eta.i, Data$Time[Data$Subject == i])
     Xtilde = rbind(Xtilde, dmu.i$dmu1 %*% A, dmu.i$dmu2 %*% A)
     if(i == 1) idx = 1: cumsum.oi[i]
     else idx = (cumsum.oi[i-1]+1): cumsum.oi[i]
     TrZ[idx, ((i-1)*q+1): (i*q)] = rbind(dmu.i$dmu1 %*% B, dmu.i$dmu2 %*% B)
     SigCi[[i]] = kronecker(Sigma, DEC(phi, ga, Data$Time[Data$Subject == i]))
     SigCi.inv[[i]] = solve(SigCi[[i]])
     TLam[idx, idx] = TrZ[idx, ((i-1)*q+1): (i*q)] %*% DD %*% t(TrZ[idx, ((i-1)*q+1): (i*q)])+ SigCi[[i]]
   }
   Ytilde = Y - MU + Xtilde %*% beta + TrZ %*% as.vector(t(b))
   loglik.iter[k, 1] = Tay.aloglik.na(Ytilde, Xtilde, TrZ, beta, DD, Sigma, phi, ga, cumsum.oi, cumsum.ni)
   cat(rep('=', 25), 'MNLMM.na with ', typeCi, ' errors; ',  '%; missing = ', na.num/n*100, '%', rep('=', 25), sep = '', '\n')
   cat(paste(rep('-', 20), sep = '', collapse = ''), 'MCMC: Chain = ', k, '\t iter = 0, loglik = ', loglik.iter[k, 1], paste(rep('-', 20), sep = '', collapse = ''), '\n')
   cat('Initial values: beta = ',  round(beta, 2), ',\t D = ', round(DD[vech.DD], 2), '\t Sigma = ', round(Sigma[vech.Sig], 3), ',\t (phi,ga) = ', round(c(phi,ga), 2), '\n')

   Theta[1,,k] = c(beta, DD[vech.DD], Sigma[vech.Sig], phi, ga, loglik.iter[k, 1])
   TB[1,,k] = as.vector(t(b))
   TYfit[1, , k] = MU
   for(i in 1: num.mis.sub){
     if(mis.sub[i] == 1) idx = 1: cumsum.oi[i]
     if(mis.sub[i] != 1) idx = (cumsum.oi[mis.sub[i]-1]+1): cumsum.oi[mis.sub[i]]
     TYimp[1,(na.num[i]+1):na.num[(i+1)],k] = mean(Y.na[idx], na.rm=T)
   }
   for(iter in 2: ITER)
   {
# imputation for ym
# Ym & b (pseudo data)
     D.inv = solve(DD)
     cent = Ytilde - Xtilde %*% beta
     for(i in 1: num.obs.sub){
       if(obs.sub[i] == 1) idx = 1: cumsum.oi[i]
       if(obs.sub[i] != 1) idx = (cumsum.oi[obs.sub[i]-1]+1): cumsum.oi[obs.sub[i]]
       ytil[[obs.sub[i]]] = Ytilde[idx]
       TSoo[idx, idx] = solve(TLam[idx, idx])
       Sb = solve(t(TrZ[idx, ((obs.sub[i]-1)*q+1): (obs.sub[i]*q)]) %*% SigCi.inv[[obs.sub[i]]] %*% TrZ[idx, ((obs.sub[i]-1)*q+1): (obs.sub[i]*q)] + D.inv)
       mu.b = Sb %*% t(TrZ[idx, ((obs.sub[i]-1)*q+1): (obs.sub[i]*q)]) %*% SigCi.inv[[obs.sub[i]]] %*% cent[idx]
       b[obs.sub[i],] = rmvnorm(1, mu.b, Sb)
     }
     for(i in 1: num.mis.sub){
       Oi = matrix(diag(oi[mis.sub[i]])[-which(mis.ind[[i]]==T), ], ncol=oi[mis.sub[i]])
       Mi = matrix(diag(oi[mis.sub[i]])[which(mis.ind[[i]]==T), ], ncol=oi[mis.sub[i]])
       if(mis.sub[i] == 1) idx = 1: cumsum.oi[i]
       if(mis.sub[i] != 1) idx = (cumsum.oi[mis.sub[i]-1]+1): cumsum.oi[mis.sub[i]]
       Soo = t(Oi) %*% solve(TLam[idx, idx][-which(mis.ind[[i]]==T), -which(mis.ind[[i]]==T)]) %*% Oi
       mu2.1 = Mi %*% (Xtilde[idx, ]%*%beta + TLam[idx,idx]%*% Soo %*% cent[idx])
       S22.1 = Mi %*% TLam[idx, idx] %*% (diag(oi[mis.sub[i]])-Soo %*% TLam[idx,idx]) %*% t(Mi)
       ym.til = c(rmvnorm(1, mu2.1, S22.1))
       ytil[[mis.sub[i]]] = c(t(Oi) %*% as.matrix(Ytilde[idx][-which(mis.ind[[i]]==T)]) + t(Mi) %*% ym.til)
       TSoo[idx, idx] = Soo
       Joo = t(Oi) %*% solve(SigCi[[mis.sub[i]]][-which(mis.ind[[i]]==T), -which(mis.ind[[i]]==T)]) %*% Oi
       Sb = solve(t(TrZ[idx, ((mis.sub[i]-1)*q+1): (mis.sub[i]*q)]) %*% Joo %*% TrZ[idx, ((mis.sub[i]-1)*q+1): (mis.sub[i]*q)] + D.inv)
       mu.b = Sb %*% t(TrZ[idx, ((mis.sub[i]-1)*q+1): (mis.sub[i]*q)]) %*% Joo %*% cent[idx]
       b[mis.sub[i], ] = rmvnorm(1, mu.b, Sb)
    }
    Y.til.hat = NULL
    for(i in 1: N) Y.til.hat = c(Y.til.hat, ytil[[i]])
    Y.hat = Y.til.hat + MU - Xtilde %*% beta - TrZ %*% as.vector(t(b))
    TYimp[iter,,k] = Y.hat[na.ind]
    TYfit[iter, , k] = MU

# b
     TB[iter,,k] = as.vector(t(b))
# beta (pseudo data)
     Sig.beta = solve(t(Xtilde) %*% TSoo %*% Xtilde + J0.inv)
     mu.beta = Sig.beta %*% (t(Xtilde) %*% TSoo %*% Ytilde + J0.inv %*% beta0)
     beta = c(rmvnorm(1, mu.beta, Sig.beta))
# DD
     Mb = matrix(b, ncol = N)
     DD =  riwish(N+d0, t(b)%*%b+G0)
# Sigma
     MU = NULL
     for(i in 1: N){
       eta.i = A %*% beta + B %*% b[i, ]
       mu.i = mu.fn(eta.i, Data$Time[Data$Subject == i])
       MU = c(MU, mu.i$mu1, mu.i$mu2)
     }
     Y.cent = Y.hat - MU
     e.c = NULL
     for(i in 2: N) e.c = c(e.c, rep(0, n), Y.cent[(cumsum.oi[i-1]+1): cumsum.oi[i]])
     e.c = matrix(c(Y.cent[1: cumsum.oi[1]], e.c), ncol = N)
     TE = e.c %*% t(e.c)
     Psi = diag(r)
     for(j in 1: r)for(l in 1: r)
     {
       Ce=0
       Ce=sum(solve(DEC(phi, ga, Data$Time[Data$Subject == 1]))*TE[1: cumsum.oi[1], 1: cumsum.oi[1]][((j-1)*si[1]+1): (j*si[1]), ((l-1)*si[1]+1): (l*si[1])])
       for(i in 2: N) Ce = Ce + sum(solve(DEC(phi, ga, Data$Time[Data$Subject == i]))*TE[(cumsum.oi[i-1]+1): cumsum.oi[i], (cumsum.oi[i-1]+1): cumsum.oi[i]][((j-1)*si[i]+1): (j*si[i]), ((l-1)*si[i]+1): (l*si[i])])
       Psi[j, l] = Ce
     }
     Sigma = riwish(sum(si)+s0, Psi+H0)
# phi and ga
     if(typeCi == 'UNC') phi = 1e-6; ga = 1
     if(typeCi == 'AR1'){
        ga = 1
        phi = MH.Phi(phi, ga, Y.cent, typeCi, cumsum.oi, cumsum.ni, r, Phi.mle)
     }
     if(typeCi == 'DEC'){
        Phi = MH.Phi(phi, ga, Y.cent, typeCi, cumsum.oi, cumsum.ni, r, Phi.mle)
        phi = Phi[1]; ga = Phi[2]
     }
     pos.samp = c(beta, DD[vech.DD], Sigma[vech.Sig], phi, ga)
# approximate log-likelihood
     MU = Xtilde = NULL
     for(i in 1: N){
       eta.i = A %*% beta + B %*% b[i, ]
       mu.i = mu.fn(eta.i, Data$Time[Data$Subject == i])
       MU = c(MU, mu.i$mu1, mu.i$mu2)
       dmu.i = dmu(eta.i, Data$Time[Data$Subject == i])
       Xtilde = rbind(Xtilde, dmu.i$dmu1 %*% A, dmu.i$dmu2 %*% A)
       if(i == 1) idx = 1: cumsum.oi[i]
       else idx = (cumsum.oi[i-1]+1): cumsum.oi[i]
       TrZ[idx, ((i-1)*q+1): (i*q)] = rbind(dmu.i$dmu1 %*% B, dmu.i$dmu2 %*% B)
       SigCi[[i]] = kronecker(Sigma, DEC(phi, ga, Data$Time[Data$Subject == i]))
       SigCi.inv[[i]] = solve(SigCi[[i]])
       TLam[idx, idx] = TrZ[idx, ((i-1)*q+1): (i*q)] %*% DD %*% t(TrZ[idx, ((i-1)*q+1): (i*q)])+ SigCi[[i]]
     }
     Ytilde = Y - MU + Xtilde %*% beta + TrZ %*% as.vector(t(b))
     loglik.iter[k, iter] = Tay.aloglik.na(Ytilde, Xtilde, TrZ, beta, DD, Sigma, phi, ga, cumsum.oi, cumsum.ni)
     if(iter%%per == 0)  cat('chain = ', k, '; iter = ', iter, '; loglik = ', loglik.iter[k, iter], ',\t pos samples = ', round(pos.samp, 3), sep = ' ', '\n')
     Theta[iter, , k] = c(pos.samp, loglik.iter[k, iter])

#     if(PATH != NULL){
#      write(c(pos.samp, loglik.iter[k, iter]), paste(PATH,'theta.txt',sep=""), ncol=length(init)+1, append=T)
#      write(as.vector(t(b)), paste(PATH,'b.txt',sep=""), ncol=(q*N), append=T)
#     }
   }}
   end = proc.time()[1]
   run.sec = end - begin
   cat('It took', run.sec, 'seconds.\n')
# posterior results
   CHECK = MPSRF.multi(Theta[, 1: (m-2), ], max.bins=50, start=0)   
   burn.in = CHECK$conv.step 
   if(burn.in == 'NA') burn.in = ITER/2
   cho = seq(0, ITER-burn.in, nlag)[-1]
   size = length(cho)

   theta.hat = apply(Theta[-c(1:burn.in),-(m+1),][cho,,], 2, mean)
   pos.sd = apply(Theta[-c(1:burn.in),-(m+1),][cho,,], 2, sd)
   pos.CI = apply(Theta[-c(1:burn.in),-(m+1),][cho,,], 2, quantile, prob=c(0.025, 0.975, 0.5))
   cp01 = as.numeric(pos.CI[1,]<=par.true & pos.CI[2,]>=par.true)
   len = pos.CI[2,] - pos.CI[1,]
   theta.out = rbind(theta.hat, pos.sd, pos.CI, len, cp01)
   rownames(theta.out) = c('Mean','SD','2.5%','97.5%','Median','CI.len','CP01')
   b.hat = matrix(apply(TB[-c(1:burn.in),,][cho,,], 2, mean), nrow=N, byrow=T)
   yimp = apply(TYimp[-c(1:burn.in),,][cho,,], 2, mean)
 #  pos.summ = list(theta.hat=theta.hat, b.hat=b.hat, yimp=yimp)

# DIC
   D.bar = -2 * mean(loglik.iter[,-c(1:burn.in)][,cho])
   beta.hat = theta.hat[1:p]
   Sig.hat = diag(r); D.hat = diag(q)
   inv.vech.DD = cbind(vech.DD[,2], vech.DD[,1])
   inv.vech.Sig = cbind(vech.Sig[,2], vech.Sig[,1])
   D.hat[vech.DD] = D.hat[inv.vech.DD] = theta.hat[-c(1:p)][1:m1]
   Sig.hat[vech.Sig] = Sig.hat[inv.vech.Sig] = theta.hat[-c(1:(p+m1))][1:m2]
   phi.hat = rev(theta.hat)[2]
   ga.hat = rev(theta.hat)[1]
   MU = Xtilde = NULL
   for(i in 1: N){
     eta.i = A %*% beta.hat + B %*% b.hat[i, ]
     mu.i = mu.fn(eta.i, Data$Time[Data$Subject == i])
     MU = c(MU, mu.i$mu1, mu.i$mu2)
     dmu.i = dmu(eta.i, Data$Time[Data$Subject == i])
     Xtilde = rbind(Xtilde, dmu.i$dmu1 %*% A, dmu.i$dmu2 %*% A)
     if(i == 1) idx = 1: cumsum.oi[i]
     else idx = (cumsum.oi[i-1]+1): cumsum.oi[i]
     TrZ[idx, ((i-1)*q+1): (i*q)] = rbind(dmu.i$dmu1 %*% B, dmu.i$dmu2 %*% B)
   }
   Ytilde = Y - MU + Xtilde %*% beta.hat + TrZ %*% as.vector(t(b.hat))
   D.theta.bar = -2 * Tay.aloglik.na(Ytilde, Xtilde, TrZ, beta.hat, D.hat, Sig.hat, phi.hat, ga.hat, cumsum.oi, cumsum.ni)
   DIC = 2*D.bar - D.theta.bar
   para.pos = list(theta.hat=theta.hat, beta=beta.hat, Sigma=Sig.hat, DD=D.hat, Phi=c(phi.hat, ga.hat), MU)
   if(typeCi == 'UNC') ss = 0
   if(typeCi == 'AR1') ss = 1
   if(typeCi == 'DEC') ss = 2
   m1 = q * (q+1) / 2; m2 = r * (r+1) / 2
   num.par = p + m1 + m2 + ss
   EAIC = D.bar + 2* num.par
   EBIC = D.bar + num.par * log(N)
# LPML
   new.ITER = (ITER - burn.in)/nlag
   theta.new = Theta[-c(1:burn.in),,][cho,,]
   TB.new = TB[-c(1: burn.in),,][cho,,]
   Sig.l = diag(r); DD.l = diag(q)
   ind.lik = array(0, dim=c(new.ITER, N, chain))
   for(k in 1: chain){
   for(iter in 1: new.ITER){
   MU = Xtilde = NULL
   bb = matrix(TB.new[iter, ,k], nrow=N, byrow=T)
   for(i in 1: N){
     eta.i = A %*% theta.new[iter, 1:p, k] + B %*% bb[i, ]
     mu.i = mu.fn(eta.i, Data$Time[Data$Subject == i])
     MU = c(MU, mu.i$mu1, mu.i$mu2)
     dmu.i = dmu(eta.i, Data$Time[Data$Subject == i])
     Xtilde = rbind(Xtilde, dmu.i$dmu1 %*% A, dmu.i$dmu2 %*% A)
     if(i == 1) idx = 1: cumsum.oi[i]
     else idx = (cumsum.oi[i-1]+1): cumsum.oi[i]
     TrZ[idx, ((i-1)*q+1): (i*q)] = rbind(dmu.i$dmu1 %*% B, dmu.i$dmu2 %*% B)
   }
   Ytilde = Y - MU + Xtilde %*% theta.new[iter, 1:p, k] + TrZ %*% as.vector(t(bb))
   Y.cent = Ytilde - Xtilde %*% theta.new[iter, 1:p, k]
   beta.l = theta.new[iter, 1:p, k]
   DD.l[vech.DD] = DD.l[inv.vech.DD] = theta.new[iter, -c(1:p), k][1:m1]
   Sig.l[vech.Sig] = Sig.l[inv.vech.Sig] = theta.new[iter, -c(1:(p+m1)), k][1:m2]
   phi.l = rev(theta.new[iter, ,k])[3]
   ga.l = rev(theta.new[iter,,k])[2]
   ind.lik[iter,,k] = Tay.alik.ind.na(Y.cent, TrZ, beta.l, DD.l, Sig.l, phi.l, ga.l, ni, cumsum.oi, cumsum.ni)
   }}
   CPO = 1 / apply(1/ind.lik, 2, mean)
   LPML = sum(log(CPO))
   KL = -log(CPO) + apply(ind.lik, 2, mean)
   model.inf = list(m=num.par, DIC=DIC, EAIC=EAIC, EBIC=EBIC, CPO=CPO, LPML=LPML, ind.lnL=ind.lik, KL=KL)

# ym imputation
   Ym.hat = apply(TYimp[-c(1:burn.in),,][cho,,], 2, mean)
   Ym.sd = apply(TYimp[-c(1:burn.in),,][cho,,], 2, sd)
   Ym.CI = apply(TYimp[-c(1:burn.in),,][cho,,], 2, quantile, prob=c(0.025, 0.975, 0.5))
   Ym.out = rbind(Ym.hat, Ym.sd, Ym.CI)
# fitted responses 
   Yfit = apply(TYfit[-c(1:burn.in),,][cho,,], 2, mean)
   Yfit.sd = apply(TYfit[-c(1:burn.in),,][cho,,], 2, sd)
   Yfit.CI = apply(TYfit[-c(1:burn.in),,][cho,,], 2, quantile, prob=c(0.025, 0.975, 0.5))
   yfit.out = rbind(Yfit, Yfit.sd, Yfit.CI)
   rownames(Ym.out) = rownames(yfit.out) = c('Mean','SD','2.5%','97.5%','Median')  
   pos.summ = list(theta.out=theta.out, Ym.out=Ym.out, yfit.out=yfit.out, b.hat=b.hat)
   return(list(run.sec = run.sec, Theta = Theta, model.inf = model.inf, TB = TB, TYimp=TYimp, para.pos = para.pos, est = theta.out, pos.summ=pos.summ, alnL = loglik.iter, Init=Init, burn.in=burn.in, size=size*chain, CPUtime=end-begin))
}

Tay.alik.ind.na = function(Y.cent, TrZ, beta, DD, Sigma, phi, ga, ni, cumsum.oi, cumsum.ni)
{
   N=length(unique(Data$Subject))
   ind.lik = numeric(N)
   Y.na = vecData$Resp
   na.ind = which(is.na(Y.na))
   n = length(Y.na)
   q = ncol(DD)
   Yo.cent = Y.cent[-na.ind]
   TLambda = matrix(0, ncol=n, nrow=n)
   for(i in 1: N){
     if(i == 1) idx = 1: cumsum.oi[1]
     else idx = (cumsum.oi[i-1]+1): cumsum.oi[i]
     rZi = TrZ[idx, ((i-1)*q+1): (i*q)]
     TLambda[idx, idx] = rZi %*% DD %*% t(rZi)+ kronecker(Sigma, DEC(phi, ga, Data$Time[Data$Subject == i]))
   }
   TLambda.oo = TLambda[-na.ind, -na.ind]
   for(i in 1: N){
#    cat('i=', i, '\n')
    if(i == 1){ idx.o = 1: cumsum.ni[i]
    } else idx.o = (cumsum.ni[i-1]+1): cumsum.ni[i]
    Lam.oo.inv = solve(TLambda.oo[idx.o, idx.o])
    ind.lik[i] = exp(0.5* (log(det(Lam.oo.inv)) - log(2*pi)*ni[i] - t(Yo.cent[idx.o]) %*% Lam.oo.inv %*% Yo.cent[idx.o]))
   }
   return(ind.lik)
}
