library(nlme)
library(mvtnorm)
library(MCMCpack)

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
     sd.phi.star = 0.1* sqrt((1/(Phi.mle[1]*(1-Phi.mle[1])))^2)
     Phi.star.new = rnorm(1, Phi.star.old, sd.phi.star)
     Phi.new = exp(Phi.star.new) / (1 + exp(Phi.star.new))
   }
   if(typeCi == 'DEC'){
     Phi.old = c(phi, ga)
     Phi.star.old = c(log(Phi.old[1]/(1-Phi.old[1])), log(Phi.old[2]))
     cov.Phi.star = 0.005*matrix(c((1/(Phi.mle[1]*(1-Phi.mle[1])))^2, 0, 0, 1/Phi.mle[2]^2), 2,2)
     repeat{
     Phi.star.new = rmvnorm(1, Phi.star.old, round(cov.Phi.star,6))
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

log.pos.den.n = function(b, Ytilde, beta, DD, Sigma, phi, ga, beta0, J0, d0, G0, s0, H0, Xtilde, TrZ, cumsum.oi)
{
  q = ncol(DD)
  SigCi = kronecker(Sigma, DEC(phi, ga, Data$Time[Data$Subject == 1]))
  y.b = log(dmvnorm(Ytilde[1: cumsum.oi[1]], mean = Xtilde[1: cumsum.oi[1], ] %*% beta + TrZ[1: cumsum.oi[1], ((i-1)*q+1): (i*q)] %*% b[1: q], sigma = SigCi) * dmvnorm(b[1: q], mean = rep(0, q), sigma = DD))
  for(i in 2: N){
      SigCi = kronecker(Sigma, DEC(phi, ga, Data$Ti[Data$Subject == i]))
      y.b = y.b + log(dmvnorm(Ytilde[(cumsum.oi[i-1]+1): cumsum.oi[i]], mean = Xtilde[(cumsum.oi[i-1]+1): cumsum.oi[i], ] %*% beta + 
            TrZ[(cumsum.oi[i-1]+1): cumsum.oi[i], ((i-1)*q+1): (i*q)] %*% b[((i-1)*q+1): (i*q)], sigma = SigCi) * dmvnorm(b[((i-1)*q+1): (i*q)], mean = rep(0, q), sigma = DD))
  }
  log.poster.n = y.b + log(dmvnorm(c(beta), beta0, J0) * diwish(DD, d0, G0) * diwish(Sigma, s0, H0)) - 2*log(1+ga)
  return(log.poster.n)
}

IBF.b.ym = function(ITER, par.mle, hyper, Xtilde, TrZ)
{
   begin=proc.time()[1]
   J = 2 * ITER
   beta = par.mle$beta; p = length(beta)
   DD = par.mle$DD; q= ncol(DD)
   D.inv = solve(DD)
   Sigma = par.mle$Sigma 
   phi = par.mle$phi
   ga = par.mle$ga
   b = hyper$mle.b
   beta0 = hyper$beta0 
   J0 = hyper$J0 
   d0 = hyper$d0 
   G0 = hyper$G0 
   s0 = hyper$s0 
   H0 = hyper$H0

   N = length(unique(Data$Subject))
   Y.na = vecData$Resp
   na.ind = which(is.na(Y.na))
   n = length(Y.na)
   Y = Y.na
   Y[na.ind] = 9999
   A = matrix(0, 6, p)
   B = diag(6)
   TrZ = matrix(0, ncol=N*q, nrow=n)
   SigCi = SigCi.inv = as.list(N)
   TLam = matrix(0, nrow=n, ncol=n)
   TB = matrix(NA, nrow = J, ncol = N*q)
   Ym = matrix(NA, nrow = J, ncol = length(na.ind))
   log.pos.den = rep(0, J)
   ti = matrix(NA, 2, N)
   for(i in 1: N){
    ti[1, i] = length(na.omit(Data$lrna[Data$Subject == i]))
    ti[2, i] = length(na.omit(Data$cd4.cd8[Data$Subject == i]))
   }
   ni = colSums(ti)
   cumsum.ni = cumsum(ni)
   si = numeric(N)
   for(i in 1: N) si[i] = length(Data$Time[Data$Subject == i])
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

   MU = Xtilde = NULL
   for(i in 1: N){
     A[1, 1:2] = A[2, 3:4] = A[5, 9:10] = A[6, 11:12] = c(1, Data$arm[which(Data$Subject==i)][1])
     A[3, 5:6] = A[4, 7:8] = c(1, Data$lcd4[which(Data$Subject==i)][1])
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
   cent = Ytilde - Xtilde %*% beta
   bi = ym = as.list(N)
 # Ym & b #
   for(j in 1: J){
     for(i in 1: num.obs.sub){
       if(obs.sub[i] == 1) idx = 1: cumsum.oi[i]
       if(obs.sub[i] != 1) idx = (cumsum.oi[obs.sub[i]-1]+1): cumsum.oi[obs.sub[i]]
       Sb = solve(t(TrZ[idx, ((obs.sub[i]-1)*q+1): (obs.sub[i]*q)]) %*% SigCi.inv[[obs.sub[i]]] %*% TrZ[idx, ((obs.sub[i]-1)*q+1): (obs.sub[i]*q)] + D.inv)
       mu.b = Sb %*% t(TrZ[idx, ((obs.sub[i]-1)*q+1): (obs.sub[i]*q)]) %*% SigCi.inv[[obs.sub[i]]] %*% cent[idx]
       bi[[obs.sub[i]]] = rmvnorm(1, mu.b, round(Sb, 3))
     }
     for(i in 1: num.mis.sub){
       Oi = matrix(diag(oi[mis.sub[i]])[-which(mis.ind[[i]]==T), ], ncol=oi[mis.sub[i]])
       Mi = matrix(diag(oi[mis.sub[i]])[which(mis.ind[[i]]==T), ], ncol=oi[mis.sub[i]])
       if(mis.sub[i] == 1) idx = 1: cumsum.oi[i]
       if(mis.sub[i] != 1) idx = (cumsum.oi[mis.sub[i]-1]+1): cumsum.oi[mis.sub[i]]
       Soo = t(Oi) %*% solve(TLam[idx, idx][-which(mis.ind[[i]]==T), -which(mis.ind[[i]]==T)]) %*% Oi
       Joo = t(Oi) %*% solve(SigCi[[mis.sub[i]]][-which(mis.ind[[i]]==T), -which(mis.ind[[i]]==T)]) %*% Oi
       if(mis.sub[i] == 1) idx = 1: cumsum.oi[i]
       if(mis.sub[i] != 1) idx = (cumsum.oi[mis.sub[i]-1]+1): cumsum.oi[mis.sub[i]]
       mu2.1 = Mi %*% (Xtilde[idx, ]%*%beta + TLam[idx,idx]%*% Soo %*% cent[idx])
       S22.1 = Mi %*% TLam[idx, idx] %*% (diag(oi[mis.sub[i]])-Soo %*% TLam[idx,idx]) %*% t(Mi)
       ym[[i]] = rmvnorm(1, mu2.1, S22.1)
       Sb = solve(t(TrZ[idx, ((mis.sub[i]-1)*q+1): (mis.sub[i]*q)]) %*% Joo %*% TrZ[idx, ((mis.sub[i]-1)*q+1): (mis.sub[i]*q)] + D.inv)
       mu.b = Sb %*% t(TrZ[idx, ((mis.sub[i]-1)*q+1): (mis.sub[i]*q)]) %*% Joo %*% cent[idx]
       bi[[mis.sub[i]]] = rmvnorm(1, mu.b, round(Sb, 3))
     }
     
     bb = ymm = NULL
     for(i in 1: N) bb = c(bb, bi[[i]])
     for(i in 1: num.mis.sub) ymm = c(ymm, ym[[i]])          
     TB[j, ] = bb                   
     Ym[j, ] = ymm
     Ytilde[na.ind] = ymm
     log.pos.den[j] = log.pos.den.n(bb, Ytilde, beta, DD, Sigma, phi, ga, beta0, J0, d0, G0, s0, H0, Xtilde, TrZ, cumsum.oi)
     if(j%%1000 == 0) cat('Imputation: iter = ', j, 'log.pos.den = ', log.pos.den[j], '\n')
   }
   Const = - max(abs(log.pos.den))
   weight = (log.pos.den + Const)/sum((log.pos.den + Const))
   ibf.samp = sample(1:J, size = ITER, replace = F, prob = weight)
   b.IBF = TB[ibf.samp, ]
   ym.IBF = Ym[ibf.samp, ]
   end = proc.time()[1]
   cat('The imputation step took', end - begin, 'seconds.\n')
   return(list(b.IBF=b.IBF, ym.IBF=ym.IBF))
}

IBF.MLMMna = function(chain, ITER, hyper, typeCi=c('UNC','AR1','DEC'), nlag=10, per=100)
{
   begin=proc.time()[1]
# hyper parameters:
   beta0 = hyper$beta0; p = length(beta0) 
   J0 = hyper$J0
   J0.inv = solve(hyper$J0) 
   d0 = hyper$d0 
   G0 = hyper$G0; q = ncol(G0) 
   s0 = hyper$s0 
   H0 = hyper$H0
   FIpg = hyper$FIpg
   Phi.mle = hyper$pg.mle
   par.mle = hyper$par.mle
   m1 = q * (q+1) / 2; m2 = r * (r+1) / 2
   m = p + m1 + m2 + 2

   N = length(unique(Data$Subject))
   Y.na = vecData$Resp
   na.ind = which(is.na(Y.na))
   n = length(Y.na)
   num.na = length(na.ind)
   Y = Y.na
   Y[na.ind] = 9999
   A = matrix(0, 6, p)
   B = diag(6)
   TrZ = matrix(0, ncol=N*q, nrow=n)
   SigCi = SigCi.inv = as.list(N)
   TLam = matrix(0, nrow=n, ncol=n)
   ti = matrix(NA, 2, N)
   for(i in 1: N){
    ti[1, i] = length(na.omit(Data$lrna[Data$Subject == i]))
    ti[2, i] = length(na.omit(Data$cd4.cd8[Data$Subject == i]))
   }
   ni = colSums(ti)
   cumsum.ni = cumsum(ni)
   si = numeric(N)
   for(i in 1: N) si[i] = length(Data$Time[Data$Subject == i])
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

# save posterior samples:
   Theta = array(NA, dim=c(ITER, (m+1), chain))
   TB = array(NA, dim=c(ITER, N*q, chain))
   TYimp = array(NA, dim=c(ITER, num.na, chain))
   TYfit = array(NA, dim=c(ITER, n, chain))
   ytil = as.list(N)
   loglik.iter = matrix(NA, nrow=chain, ncol=ITER)
   acc.phi = acc.ga = matrix(NA, nrow=chain, ncol=ITER)
   acc.b = array(NA, dim=c(N, ITER, chain))
   vech.DD = vech.posi((q))
   vech.Sig = vech.posi(r)
   Init = NULL

   for(k in 1: chain){
# initial values: (generate from the priors)
   beta = hyper$par.mle$beta
   DD = hyper$par.mle$DD
   Sigma = hyper$par.mle$Sigma 
   phi = hyper$par.mle$phi
   ga = hyper$par.mle$ga
   b = hyper$mle.b
   if(r==1 & q==1){ Init = rbind(Init, c(beta, DD, Sigma, phi, ga))
   } else Init = rbind(Init, c(beta, DD[vech.DD], Sigma[vech.Sig], phi, ga))
    
   MU = Xtilde = NULL
   for(i in 1: N){
     A[1, 1:2] = A[2, 3:4] = A[5, 9:10] = A[6, 11:12] = c(1, Data$arm[which(Data$Subject==i)][1])
     A[3, 5:6] = A[4, 7:8] = c(1, Data$lcd4[which(Data$Subject==i)][1])
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
   cat(paste(rep('-', 20), sep = '', collapse = ''), 'Chain = ', k, '\t iter = 0, loglik = ', loglik.iter[k, 1], paste(rep('-', 20), sep = '', collapse = ''), '\n')
   cat('Initial values: beta = ',  round(beta, 2), ',\t D = ', round(DD[vech.DD], 2), '\t Sigma = ', round(Sigma[vech.Sig], 3), ',\t (phi,ga) = ', round(c(phi,ga), 2), '\n')
# Imputation Step
 # Ym & b #
   par.mle = hyper$par.mle
   ibf.samp = IBF.b.ym(ITER, par.mle, hyper, Xtilde, TrZ)
   TB[,, k] = ibf.samp$b.IBF
   TYimp[,, k] = ibf.samp$ym.IBF
# Posterior Step via Hetropolis within Gibbs sampler
   for(iter in 1: ITER)
   {
     Ytilde[na.ind] = TYimp[iter, , k] 
     cent = Ytilde - Xtilde %*% beta
     for(i in 1: num.obs.sub){
       if(obs.sub[i] == 1) idx = 1: cumsum.oi[i]
       if(obs.sub[i] != 1) idx = (cumsum.oi[obs.sub[i]-1]+1): cumsum.oi[obs.sub[i]]
       TSoo[idx, idx] = solve(TLam[idx, idx])
     }
     for(i in 1: num.mis.sub){
       if(mis.sub[i] == 1) idx = 1: cumsum.oi[i]
       if(mis.sub[i] != 1) idx = (cumsum.oi[mis.sub[i]-1]+1): cumsum.oi[mis.sub[i]]
       Oi = diag(oi[mis.sub[i]])[-which(mis.ind[[i]]==T), ]
       TSoo[idx, idx] = t(Oi) %*% solve(TLam[idx, idx][-which(mis.ind[[i]]==T), -which(mis.ind[[i]]==T)]) %*% Oi
     }
 # beta #   
     Sig.beta = solve(t(Xtilde) %*% TSoo %*% Xtilde + J0.inv)
     mu.beta = Sig.beta %*% (t(Xtilde) %*% TSoo %*% Ytilde + J0.inv %*% beta0)
     beta = c(rmvnorm(1, mu.beta, round(Sig.beta, 3)))
 # D #
     Mb = matrix(TB[iter,,k], ncol = N)
     DD = riwish(N+d0, Mb%*%t(Mb)+G0)
 # Sigma #
     MU = NULL
     for(i in 1: N){
       eta.i = A %*% beta + B %*% Mb[, i]
       mu.i = mu.fn(eta.i, Data$Time[Data$Subject == i])
       MU = c(MU, mu.i$mu1, mu.i$mu2)
     }
     Y.hat = Ytilde + MU - Xtilde %*% beta - TrZ %*% as.vector(t(TB[iter,,k]))
     Y.cent = Y.hat - MU
     e = Y.cent - TrZ %*% TB[iter,, k]
     e.c = NULL
     for(i in 2: N) e.c = c(e.c, rep(0, n), e[(cumsum.oi[i-1]+1): cumsum.oi[i]])
     e.c = matrix(c(e[1: cumsum.oi[1]], e.c), ncol = N)
     TE = e.c %*% t(e.c) 
     Psi = diag(r)
     if(r == 1){
       Ce = 0
       Ce = sum(solve(DEC(phi, ga, Data$Time[Data$Subject == 1])) * TE[1: cumsum.oi[1], ][1: oi[1], 1: oi[1]])
       for(i in 2: N) Ce = Ce + sum(solve(DEC(phi, ga, Data$Time[Data$Subject == i]))*TE[(cumsum.oi[i-1]+1): cumsum.oi[i], (cumsum.oi[i-1]+1): cumsum.oi[i]][1: oi[i], 1: oi[i]])
       Psi = Ce
     }
    else{
     for(j in 1: r)for(l in 1: r){
       Ce = 0
       Ce=sum(solve(DEC(phi, ga, Data$Time[Data$Subject == 1]))*TE[1: cumsum.oi[1], 1: cumsum.oi[1]][((j-1)*si[1]+1): (j*si[1]), ((l-1)*si[1]+1): (l*si[1])])
       for(i in 2: N) Ce = Ce + sum(solve(DEC(phi, ga, Data$Time[Data$Subject == i]))*TE[(cumsum.oi[i-1]+1): cumsum.oi[i], (cumsum.oi[i-1]+1): cumsum.oi[i]][((j-1)*si[i]+1): (j*si[i]), ((l-1)*si[i]+1): (l*si[i])])
       Psi[j, l] = Ce
     }}
     Sigma = riwish(sum(oi)+s0, Psi+H0)
 # phi & ga #
     if(typeCi == 'UNC'){ phi=1e-6; ga = 1}
     if(typeCi == 'AR1'){
        ga = 1
        phi = MH.Phi(phi, ga, Y.cent, typeCi, cumsum.oi, cumsum.ni, r, Phi.mle)
      }
     if(typeCi == 'DEC'){   
        par.DEC = MH.Phi(phi, ga, Y.cent, typeCi, cumsum.oi, cumsum.ni, r, Phi.mle)
        phi = par.DEC[1]; ga = par.DEC[2]
     }
     if(r==1 & q ==1) pos.samp = c(beta, DD, Sigma, phi, ga)
     else pos.samp = c(beta, DD[vech.DD], Sigma[vech.Sig], phi, ga)
     MU = Xtilde = NULL
     for(i in 1: N){
       A[1, 1:2] = A[2, 3:4] = A[5, 9:10] = A[6, 11:12] = c(1, Data$arm[which(Data$Subject==i)][1])
       A[3, 5:6] = A[4, 7:8] = c(1, Data$lcd4[which(Data$Subject==i)][1])
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
     if(iter%%per == 0) cat('iter = ', iter, ',\t loglik = ', loglik.iter[k, iter], ',\t pos.samples = ', round(pos.samp, 3), '\n')
     Theta[iter, , k] = c(pos.samp, loglik.iter[k, iter])
     TYfit[iter, , k] = MU
   }}
   end = proc.time()[1]
   run.sec = end - begin
   cat('The posterior step took', run.sec, 'seconds.\n')
   if(typeCi == 'UNC') num.par = m - 2
   if(typeCi == 'AR1') num.par = m - 1
   if(typeCi == 'DEC') num.par = m
   CHECK = MPSRF.multi(Theta[, 1: p, ], max.bins=50, start=0)   
   burn.in = ifelse(CHECK$conv.step=='NA', ITER/2, CHECK$conv.step) 
   cho = seq(0, ITER-burn.in, nlag)[-1]
   size = length(cho)
   
# posterior inference for parameters:
   inv.vech.DD = cbind(vech.DD[,2], vech.DD[,1])
   inv.vech.Sig = cbind(vech.Sig[,2], vech.Sig[,1])
   r1 = nrow(inv.vech.DD); r2 = nrow(inv.vech.Sig)
# theta.hat
   theta.hat = apply(Theta[-c(1:burn.in), -(m+1),][cho,,], 2, mean)
   beta.hat = theta.hat[1: p]
   D.hat = matrix(NA, q, q)
   Sig.hat = matrix(NA, r, r)
   D.hat[vech.DD]=D.hat[inv.vech.DD] = theta.hat[(p+1): (p+r1)]
   Sig.hat[vech.Sig]=Sig.hat[inv.vech.Sig] = theta.hat[(p+r1+1):(p+r1+r2)]
   Phi.hat = theta.hat[-(1:(p+r1+r2))]
   pos.sd = apply(Theta[-c(1:burn.in),-(m+1),][cho,,], 2, sd)
   pos.CI = apply(Theta[-c(1:burn.in),-(m+1),][cho,,], 2, quantile, prob=c(0.025, 0.975, 0.5))
   len = pos.CI[2,] - pos.CI[1,]
   theta.out = rbind(theta.hat, pos.sd, pos.CI, len)
   rownames(theta.out) = c('Mean','SD','2.5%','97.5%','Median','CI.len')
# ym imputation
   Ym.hat = apply(TYimp[-c(1:burn.in),,][cho,,], 2, mean)
   Ym.sd = apply(TYimp[-c(1:burn.in),,][cho,,], 2, sd)
   Ym.CI = apply(TYimp[-c(1:burn.in),,][cho,,], 2, quantile, prob=c(0.025, 0.975, 0.5))
   theta.out = rbind(theta.hat, pos.sd, pos.CI)
   Ym.out = rbind(Ym.hat, Ym.sd, Ym.CI)
# fitted responses
   Yfit = apply(TYfit[-c(1:burn.in),,][cho,,], 2, mean)
   Yfit.sd = apply(TYfit[-c(1:burn.in),,][cho,,], 2, sd)
   Yfit.CI = apply(TYfit[-c(1:burn.in),,][cho,,], 2, quantile, prob=c(0.025, 0.975, 0.5))
   yfit.out = rbind(Yfit, Yfit.sd, Yfit.CI)
   rownames(theta.out) = rownames(Ym.out) = rownames(yfit.out) = c('Mean','SD','2.5%','97.5%','Median')
# random effects estimation
   b.hat = matrix(apply(TB[-c(1:burn.in),,][cho,,], 2, mean), nrow=N, byrow=T)
   pos.summ = list(theta.out=theta.out, Ym.out=Ym.out, yfit.out=yfit.out, b.hat=b.hat)

# model selection criteria:
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
     A[1, 1:2] = A[2, 3:4] = A[5, 9:10] = A[6, 11:12] = c(1, Data$arm[which(Data$Subject==i)][1])
     A[3, 5:6] = A[4, 7:8] = c(1, Data$lcd4[which(Data$Subject==i)][1])
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
   para.pos = list(theta.hat=theta.hat, beta=beta.hat, Sigma=Sig.hat, DD=D.hat, Phi=c(phi.hat, ga.hat))
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
     A[1, 1:2] = A[2, 3:4] = A[5, 9:10] = A[6, 11:12] = c(1, Data$arm[which(Data$Subject==i)][1])
     A[3, 5:6] = A[4, 7:8] = c(1, Data$lcd4[which(Data$Subject==i)][1])
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
   CPO = 1 / apply(1/exp(ind.lik), 2, mean)
   LPML = sum(log(CPO))
   model.inf = list(m=m1, DIC=DIC, EAIC=EAIC, EBIC=EBIC, CPO=CPO, LPML=LPML)
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
