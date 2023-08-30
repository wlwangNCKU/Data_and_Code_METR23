library(nlme)
library(mvtnorm)
# DEC structure
DEC = function(phi, ga, Ti) phi^((abs(outer(Ti, Ti, '-')))^ga)
DEC.dot.phi = function(phi, ga, Ti)
{
  tem = abs(outer(Ti, Ti,'-'))^ga
  dot.DEC = tem*(phi^(tem-1))
  return(dot.DEC)
}
DEC.dot.ga = function(phi, ga, Ti)
{
  tem = abs(outer(Ti, Ti,'-'))
  dot.ga.DEC = log(phi^tem)*DEC(phi, ga, Ti)
  return(dot.ga.DEC)
}

# ML Estimation
mu.fn = function(eta, x)
{
  mu1 = eta[1] / (1 + exp((eta[2]-x)/eta[3]))
  mu2 = eta[4] + eta[5]*x
  return(list(mu1 = mu1, mu2 = mu2, mu = c(mu1, mu2)))
}

dmu = function(eta, x)
{
  p = length(eta)
  si = length(x)
  dmu1 = dmu2 = matrix(0, ncol=p, nrow=si)
  a = exp((eta[2]-x)/eta[3])
  dmu1[,1] = 1/(1+a)
  dmu1[,2] = -(eta[1]*a) / (eta[3]*(1+a)^2)
  dmu1[,3] = (eta[1]*(eta[2]-x)*a) / (eta[3]^2*(1+a)^2)
  dmu2[,4] = 1
  dmu2[,5] = x
  return(list(dmu1 = dmu1, dmu2 = dmu2))
}

# vech:
vech.posi=function(dim) cbind(rep(1:dim, 1:dim), unlist(mapply(':', 1, 1:dim)))

# ML Estimation
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

# phiga part in Q-function:
phi.ga.Qna=function(phiga, Sigma, Ome, typeCi='DEC', cumsum.oi)
{
   if(typeCi == 'AR1'){ phi = phiga; ga = 1}
   if(typeCi == 'DEC'){ phi = phiga[1]; ga = phiga[2]}
   Y.na = vecData$Resp
   n=length(Y.na)
   SigCi = as.list(N)
   log.sum = 0
   for(i in 1: N){
      if(i == 1) idx = 1: cumsum.oi[1]
      else idx = (cumsum.oi[i-1]+1): cumsum.oi[i]
      SigCi[[i]] = kronecker(Sigma, DEC(phi, ga, Data$Time[Data$Subject == i]))
      log.sum = log.sum + r*log(det(DEC(phi, ga, Data$Time[Data$Subject == i]))) + sum(solve(SigCi[[i]]) * Ome[[i]])
   }
   return(log.sum)
}
 
MNLMMna.ECM = function(init.par, tol=1e-6, typeCi=c('UNC','AR1','DEC'), max.iter = 10000, per=100)
{
   cat(paste(c(rep('-', 25), rep(' ', 5), 'Running ... Pseudo ECM for MNLMM', rep(' ', 5), rep('-', 30)), sep = '', collapse = ''), '\n')
   begin = proc.time()[1]
# missing information
   N=length(unique(Data$Subject))
   Y.na = vecData$Resp
   na.ind = which(is.na(Y.na))
   n = length(Y.na)
   ti = matrix(NA, 2, N)
   for(i in 1: N){
    ti[1, i] = length(na.omit(Data$ycm1[Data$Subject == i]))
    ti[2, i] = length(na.omit(Data$ycm2[Data$Subject == i]))
   }
   ni = colSums(ti)
   cumsum.ni = cumsum(ni)
   oi = numeric(N)
   for(i in 1: N) oi[i] = length(Data$Ti[Data$Subject == i])
   cumsum.oi = cumsum(oi*r)
   obs.sub = which(c(oi*r == colSums(ti)))
   mis.sub = which(c(oi*r != colSums(ti)))
   num.obs.sub = length(obs.sub)
   num.mis.sub = length(mis.sub)
   mis.ind = Yim = as.list(num.mis.sub)
   for(i in 1: num.mis.sub){
      if(mis.sub[i] == 1) mis.ind[[i]] = is.na(Y.na[1: cumsum.oi[mis.sub[i]]])
      else  mis.ind[[i]] = is.na(Y.na[(cumsum.oi[mis.sub[i]-1]+1): cumsum.oi[mis.sub[i]]])
   }
# initial values of parameters
   beta = init.par$beta
   DD = init.par$DD
   D.inv = solve(DD)
   Sigma = init.par$Sigma
   phi = init.par$phi
   ga = init.par$ga
   p = length(beta)
   q = ncol(DD)
   typeCi = typeCi[1]

   b = rmvnorm(N, rep(0, q), DD)
# basic setting
   m1 = q*(q+1)/2
   m2 = r*(r+1)/2
   vechD = vech.posi(q)
   vechS = vech.posi(r)
   iter = 0
   A = diag(p)
   B = matrix(c(1,rep(0, 7), 1, 0), ncol=q)
   TrZ = matrix(0, ncol=N*q, nrow=n)
   Lambda = SigCi = Vb = Ome = Ytil.hat = bi = as.list(N)
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
   }
   Ytilde = Y.na - MU + Xtilde %*% beta + TrZ %*% as.vector(t(b))
   old.loglik = Tay.aloglik.na(Ytilde, Xtilde, TrZ, beta, DD, Sigma, phi, ga, cumsum.oi, cumsum.ni)
   theta.old = c(beta, DD[vechD], Sigma[vechS], phi, ga)
   iter.lnL = old.loglik
   cat('iter = ', iter, ',\t obs.logli = ', old.loglik, sep = '', '\n')
   repeat{
     iter = iter + 1
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
       Lambda[[i]] = TrZ[idx, ((i-1)*q+1): (i*q)] %*% DD %*% t(TrZ[idx, ((i-1)*q+1): (i*q)])+ SigCi[[i]]
     }
     Ytilde = Y.na - MU + Xtilde %*% beta + TrZ %*% as.vector(t(b))
     Y.cent = Ytilde - Xtilde %*% beta
# EM-step:
     if(num.obs.sub != 0){
     for(i in 1: num.obs.sub){
       if(obs.sub[i] == 1) idx = 1: cumsum.oi[i]
       if(obs.sub[i] != 1) idx = (cumsum.oi[obs.sub[i]-1]+1): cumsum.oi[obs.sub[i]]
           Ytil.hat[[obs.sub[i]]] = Ytilde[idx]
           bi[[obs.sub[i]]] = DD %*% t(TrZ[idx, ((obs.sub[i]-1)*q+1): (obs.sub[i]*q)]) %*% solve(Lambda[[obs.sub[i]]]) %*% Y.cent[idx]
           Vb[[obs.sub[i]]] = solve(D.inv + t(TrZ[idx, ((obs.sub[i]-1)*q+1): (obs.sub[i]*q)]) %*% solve(SigCi[[obs.sub[i]]]) %*% TrZ[idx, ((obs.sub[i]-1)*q+1): (obs.sub[i]*q)])
     }}
     if(num.mis.sub != 0){
     for(i in 1: num.mis.sub){
       if(mis.sub[i] == 1) idx = 1: cumsum.oi[i]
       if(mis.sub[i] != 1) idx = (cumsum.oi[mis.sub[i]-1]+1): cumsum.oi[mis.sub[i]]
       Xim = matrix(Xtilde[idx, ][which(mis.ind[[i]]==T), ], ncol=p)
       Xio = matrix(Xtilde[idx, ][-which(mis.ind[[i]]==T), ], ncol=p)
       Zio = matrix(TrZ[idx, ((mis.sub[i]-1)*q+1): (mis.sub[i]*q)][-which(mis.ind[[i]]==T), ], ncol=q)
       Lam.oo.inv = solve(Lambda[[mis.sub[i]]][-which(mis.ind[[i]]==T), -which(mis.ind[[i]]==T)])
       Yim[[i]] = c(Xim %*% beta + Lambda[[mis.sub[i]]][which(mis.ind[[i]]==T), -which(mis.ind[[i]]==T)] %*% Lam.oo.inv %*% Y.cent[idx][-which(mis.ind[[i]]==T)])
       Oi = matrix(diag(oi[mis.sub[i]]*r)[-which(mis.ind[[i]]==T), ], ncol=oi[mis.sub[i]]*r)
       Mi = matrix(diag(oi[mis.sub[i]]*r)[which(mis.ind[[i]]==T), ], ncol=oi[mis.sub[i]]*r)
       Ytil.hat[[mis.sub[i]]] = t(Oi) %*% Ytilde[idx][-which(mis.ind[[i]]==T)] + t(Mi) %*% Yim[[i]]
       bi[[mis.sub[i]]] = DD %*% t(Zio) %*% Lam.oo.inv %*% Y.cent[idx][-which(mis.ind[[i]]==T)]
       Vb[[mis.sub[i]]] = solve(D.inv + t(Zio) %*% solve(SigCi[[mis.sub[i]]][-which(mis.ind[[i]]==T), -which(mis.ind[[i]]==T)]) %*% Zio)
     }}
     Y.hat = bb = NULL
     for(i in 1: N){
       Y.hat = c(Y.hat, Ytil.hat[[i]])
       bb = c(bb, bi[[i]])
     }
  # beta #
     term1 = t(Xtilde[1: cumsum.oi[1], ]) %*% solve(Lambda[[1]]) %*% Xtilde[1: cumsum.oi[1], ]
     term2 = t(Xtilde[1: cumsum.oi[1], ]) %*% solve(Lambda[[1]]) %*% Y.hat[1: cumsum.oi[1]]
     for(i in 2: N){
       term1 = term1 + t(Xtilde[(cumsum.oi[i-1]+1): cumsum.oi[i], ]) %*% solve(Lambda[[i]]) %*% Xtilde[(cumsum.oi[i-1]+1): cumsum.oi[i], ]
       term2 = term2 + t(Xtilde[(cumsum.oi[i-1]+1): cumsum.oi[i], ]) %*% solve(Lambda[[i]]) %*% Y.hat[(cumsum.oi[i-1]+1): cumsum.oi[i]]
     }
     beta = solve(term1) %*% term2
  # b #
     b = t(matrix(bb, ncol=N))
     sum.Vb = Vb[[1]]
     for(i in 2: N) sum.Vb = sum.Vb + Vb[[i]]
     Psi = t(b) %*% b + sum.Vb
  # DD #
     DD = Psi / N
     D.inv = solve(DD)
  # Sigma #
     e = Y.hat - Xtilde %*% beta - TrZ %*% bb
     if(num.obs.sub != 0){
     for(i in 1: num.obs.sub){
       if(obs.sub[i] == 1) Ome[[obs.sub[i]]] = e[1: cumsum.oi[i]] %*% t(e[1: cumsum.oi[i]]) + (SigCi[[i]] - SigCi[[i]] %*% solve(Lambda[[i]]) %*% SigCi[[i]])
       if(obs.sub[i] != 1) Ome[[obs.sub[i]]] = e[(cumsum.oi[obs.sub[i]-1]+1): cumsum.oi[obs.sub[i]]] %*% t(e[(cumsum.oi[obs.sub[i]-1]+1): cumsum.oi[obs.sub[i]]]) + (SigCi[[obs.sub[i]]] - SigCi[[obs.sub[i]]] %*% solve(Lambda[[obs.sub[i]]]) %*% SigCi[[obs.sub[i]]])
     }}
     if(num.mis.sub != 0){
     for(i in 1: num.mis.sub){
        Oi = matrix(diag(oi[mis.sub[i]]*r)[-which(mis.ind[[i]]==T), ], ncol=oi[mis.sub[i]]*r)
        if(mis.sub[i] == 1) Ome[[mis.sub[i]]] = e[1: cumsum.oi[i]] %*% t(e[1: cumsum.oi[i]]) + (SigCi[[i]] - SigCi[[i]] %*% t(Oi) %*% solve(Oi %*% Lambda[[i]] %*% t(Oi)) %*% Oi %*% SigCi[[i]])
        if(mis.sub[i] != 1) Ome[[mis.sub[i]]] = e[(cumsum.oi[mis.sub[i]-1]+1): cumsum.oi[mis.sub[i]]] %*% t(e[(cumsum.oi[mis.sub[i]-1]+1): cumsum.oi[mis.sub[i]]]) + (SigCi[[mis.sub[i]]] - SigCi[[mis.sub[i]]] %*% t(Oi) %*% solve(Oi %*% Lambda[[mis.sub[i]]] %*% t(Oi)) %*% Oi %*% SigCi[[mis.sub[i]]])
     }}
     if(r==1)
     {
       Ce=0
       Ce=sum(solve(DEC(phi, ga, Data$Time[Data$Subject == 1])) * Ome[[1]])
       for(i in 2: N) Ce = Ce + sum(solve(DEC(phi, ga, Data$Time[Data$Subject == i])) * Ome[[i]])
       Sigma = Ce / n
     } else{
     for(j in 1: r)for(l in 1: r)
     {
       Ce=0
       Ce=sum(solve(DEC(phi, ga, Data$Time[Data$Subject == 1])) * Ome[[1]][((j-1)*oi[1]+1): (j*oi[1]), ((l-1)*oi[1]+1): (l*oi[1])])
       for(i in 2: N) Ce = Ce + sum(solve(DEC(phi, ga, Data$Time[Data$Subject == i])) * Ome[[i]][((j-1)*oi[i]+1): (j*oi[i]), ((l-1)*oi[i]+1): (l*oi[i])])
       Sigma[j, l] = Ce /(n / r)
     }}
  # phi, ga #
     if(typeCi == 'UNC'){ phi=1e-6; ga = 1}
     if(typeCi == 'AR1'){
         ga = 1
         phi = optim(par = phi, fn = phi.ga.Qna, method = "L-BFGS-B", lower = 1e-6, upper = 1-1e-6, Sigma=Sigma, Ome=Ome, typeCi='AR1', cumsum.oi=cumsum.oi)$par
     }
     if(typeCi == 'DEC'){
        par.DEC = optim(par = c(phi, ga), fn = phi.ga.Qna, method = "L-BFGS-B", lower = c(1e-6, 1e-6), upper = c(1-1e-6, Inf), Sigma=Sigma, Ome=Ome, typeCi='DEC', cumsum.oi=cumsum.oi)$par
        phi = par.DEC[1]; ga = par.DEC[2]
     }
     new.loglik = Tay.aloglik.na(Ytilde, Xtilde, TrZ, beta, DD, Sigma, phi, ga, cumsum.oi, cumsum.ni)
     iter.lnL = c(iter.lnL, new.loglik)
     diff.lnL = (new.loglik-old.loglik)
     theta.new = c(beta, DD[vechD], Sigma[vechS], phi, ga)
     diff = max(abs((theta.new-theta.old)/theta.old))
     if(iter%%per == 0) cat('iter = ', iter, ',\t obs.loglik = ', new.loglik, ',\t lnL.diff = ', diff.lnL, ',\t theta.diff = ', diff, ',\t beta=', beta, ',\t D = ', DD[vechD], ',\t Sigma =', Sigma[vechS], ',\t phi=', phi, ',\t ga=', ga, sep = ' ', '\n')
     if(abs(diff.lnL) < tol || diff < tol || iter >= max.iter) break
     old.loglik = new.loglik; theta.old = theta.new
   }
# Summarize results:
   cat('iter = ', iter, ', \t obs.loglik = ', new.loglik, sep = '', '\n')
   end = proc.time()[1]
   cat('beta =', beta, '\n')
   cat('D =\n'); print(DD)
   cat('Sigma =\n'); print(Sigma)
   cat('phi =', phi, ',\t ga=', ga, '\n')
   if(typeCi == 'UNC') ss = 0
   if(typeCi == 'AR1') ss = 1
   if(typeCi == 'DEC') ss = 2
   num.par = p + m1 + m2 + ss
   aic = -2 * new.loglik + 2* num.par
   bic = -2 * new.loglik + num.par * log(N)
   cat('loglik = ', new.loglik, ',\t AIC = ', aic, '.\t BIC = ', bic, '\n')
   model.inf = list(loglik=new.loglik, iter.lnL = iter.lnL, aic=aic, bic=bic, num.par=num.par)
   para = list(beta = beta, DD = DD, Sigma = Sigma, phi = phi, ga = ga)
# Missing Data Imputation
   pred = IMP.FIT.mnlmm.na(para, b, N, n, p, q, A, B, obs.sub, mis.sub, num.obs.sub, num.mis.sub, mis.ind, cumsum.oi, oi)
   EST = c(beta, DD[vechD], Sigma[vechS], phi, ga)
   SD = FI.mnlmm.na(para, pred$b, typeCi, ti, oi)
   ehat = Y.na - pred$Yhat
   Ehat = matrix(ehat[1:cumsum.oi[1]], ncol=r)
   for(i in 2: N) Ehat = rbind(Ehat, matrix(ehat[(cumsum.oi[i-1]+1): cumsum.oi[i]], ncol=r))
   run.sec = end - begin
   cat("It took", run.sec, "seconds.\n")
   return(list(run.sec = run.sec, iter = iter, model.inf = model.inf, para = para, est = EST, SD = SD, b = pred$b, yfit=MU, Yhat = pred$Yhat, Yimp = pred$Yimp, yimp = pred$Yimp[na.ind], Ehat = Ehat, Del.e = pred$Del.e/ni, Del.b = pred$Del.b/q))
}

# Pseudo Fisher Information Matrix
FI.mnlmm.na = function(par, b, typeCi, ti, oi)
{
 # para.est:
   beta = par$beta
   DD = par$DD
   Sigma = par$Sigma
   phi = par$phi
   ga = par$ga
   N=length(unique(Data$Subject))
   Y.na = vecData$Resp
   na.ind = which(is.na(Y.na))
   n = length(Y.na)
   ni = colSums(ti)
   cumsum.ni = cumsum(ni)
   cumsum.oi = cumsum(oi*r)
   obs.sub = which(c(oi*r == colSums(ti)))
   mis.sub = which(c(oi*r != colSums(ti)))
   num.obs.sub = length(obs.sub)
   num.mis.sub = length(mis.sub)
   mis.ind = as.list(num.mis.sub)
   for(i in 1: num.mis.sub){
      if(mis.sub[i] == 1) mis.ind[[i]] = is.na(Y.na[1: cumsum.oi[mis.sub[i]]])
      else  mis.ind[[i]] = is.na(Y.na[(cumsum.oi[mis.sub[i]-1]+1): cumsum.oi[mis.sub[i]]])
   }
   if(phi == 1e-6 & ga == 1) typeCi = 'UNC'
   else if(ga == 1) typeCi = 'AR1'
   else typeCi = 'DEC'
   p = length(beta)
   q = ncol(DD)
   m1 = q*(q+1)/2
   m2 = r*(r+1)/2
   vechD = vech.posi(q)
   vechS = vech.posi(r)
   g = m1+m2+2
   A = diag(p)
   B = matrix(c(1,rep(0, 7), 1, 0), ncol=q)
   TrZ = matrix(0, ncol=N*q, nrow=n)
   Lambda = as.list(N)
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
     Lambda[[i]] = TrZ[idx, ((i-1)*q+1): (i*q)] %*% DD %*% t(TrZ[idx, ((i-1)*q+1): (i*q)])+ kronecker(Sigma, DEC(phi, ga, Data$Time[Data$Subject == i]))
   }
   FI = matrix(0, p+g, p+g)
 # I_{beta}:
   if(num.obs.sub != 0){
   for(i in 1: num.obs.sub){
      if(obs.sub[i] == 1) idx = 1: cumsum.oi[i]
      if(obs.sub[i] != 1) idx = (cumsum.oi[obs.sub[i]-1]+1): cumsum.oi[obs.sub[i]]
      FI[1: p, 1: p] = FI[1: p, 1: p] + t(Xtilde[idx, ]) %*% solve(Lambda[[obs.sub[i]]]) %*% Xtilde[idx, ]
   }}
   if(num.mis.sub != 0){
   for(i in 1: num.mis.sub){
      if(mis.sub[i] == 1) Xio = matrix(Xtilde[1: cumsum.oi[i], ][-which(mis.ind[[i]]==T), ], ncol=p)
      if(mis.sub[i] != 1) Xio = matrix(Xtilde[(cumsum.oi[mis.sub[i]-1]+1): cumsum.oi[mis.sub[i]], ][-which(mis.ind[[i]]==T), ], ncol=p)
        FI[1: p, 1: p] = FI[1: p, 1: p] + t(Xio) %*% solve(Lambda[[mis.sub[i]]][-which(mis.ind[[i]]==T), -which(mis.ind[[i]]==T)]) %*% Xio
   }}
# I_{omega}:
   dot.L = as.list(matrix(0, N, g))
   for(l in 1: m1)
   {
     dot.DD = matrix(0, q, q)
     dot.DD[matrix(vechD[l, ], 1)] = dot.DD[matrix(rev(vechD[l, ]), 1)] = 1
     dot.L[[(l-1)*N+1]] = TrZ[1: cumsum.oi[1], ((i-1)*q+1): (i*q)] %*% dot.DD %*% t(TrZ[1: cumsum.oi[1], ((i-1)*q+1): (i*q)])
     for(i in 2: N) dot.L[[(l-1)*N+i]] = TrZ[(cumsum.oi[i-1]+1): cumsum.oi[i], ((i-1)*q+1): (i*q)] %*% dot.DD %*% t(TrZ[(cumsum.oi[i-1]+1): cumsum.oi[i], ((i-1)*q+1): (i*q)])
   }
   for(l in 1: m2)
   {
     dot.Sig = matrix(0, r, r)
     dot.Sig[matrix(vechS[l, ], 1)] = dot.Sig[matrix(rev(vechS[l, ]), 1)] = 1
     for(i in 1: N) dot.L[[m1*N+(l-1)*N+i]] = kronecker(dot.Sig, DEC(phi, ga, Data$Time[Data$Subject == i]))
   }
   for(i in 1: N) dot.L[[(m1+m2)*N+i]] = kronecker(Sigma, DEC.dot.phi(phi, ga, Data$Time[Data$Subject == i]))
   for(i in 1: N) dot.L[[(g-1)*N+i]] = kronecker(Sigma, DEC.dot.ga(phi, ga, Data$Time[Data$Subject == i]))
   Linv.dotL = as.list(numeric(N*g))
   for(s in 1: g)for(l in 1: s){
     if(num.obs.sub != 0){
     for(i in 1: num.obs.sub){
       Linv.dotL[[(s-1)*N+obs.sub[i]]] = solve(Lambda[[obs.sub[i]]]) %*% dot.L[[(s-1)*N+obs.sub[i]]]
       FI[p+s,p+l] = FI[p+s,p+l] + sum(diag(Linv.dotL[[(s-1)*N+obs.sub[i]]] %*% Linv.dotL[[(l-1)*N+obs.sub[i]]]))/2
     }}
     if(num.mis.sub != 0){
     for(i in 1: num.mis.sub){
       Oi = matrix(diag(oi[mis.sub[i]]*r)[-which(mis.ind[[i]]==T), ], ncol=oi[mis.sub[i]]*r)
       Linv.dotL[[(s-1)*N+mis.sub[i]]] = solve(Lambda[[mis.sub[i]]][-which(mis.ind[[i]]==T), -which(mis.ind[[i]]==T)]) %*% Oi %*% dot.L[[(s-1)*N+mis.sub[i]]] %*% t(Oi)
       FI[p+s,p+l] = FI[p+s,p+l] + sum(diag(Linv.dotL[[(s-1)*N+mis.sub[i]]] %*% Linv.dotL[[(l-1)*N+mis.sub[i]]]))/2
   }}}
   for(l in (p+1):(p+g-1))for(s in (l+1):(p+g)) FI[l, s] = FI[s, l]
   if(typeCi == 'UNC') SE = c(sqrt(diag(solve(FI[1:(p+g-2), 1: (p+g-2)]))), rep(0,2))
   if(typeCi == 'AR1') SE = c(sqrt(diag(solve(FI[1:(p+g-1), 1: (p+g-1)]))), 0)
   if(typeCi == 'DEC') SE = sqrt(diag(solve(FI)))
   if(r == 1) EST=c(beta, DD, Sigma, phi, ga)
   else EST = c(beta, DD[vechD], Sigma[vechS], phi, ga)
   out = rbind(EST, SE)
   return(list(out=out, FI=round(FI, 4), Ibeta = FI[1:p, 1:p], FI.omega = FI[-(1:p), -(1:p)]))
}

# Imputation and Fitted values:
IMP.FIT.mnlmm.na = function(par, b, N, n, p, q, A, B, obs.sub, mis.sub, num.obs.sub, num.mis.sub, mis.ind, cumsum.oi, oi)
{
   Y.na = vecData$Resp
   beta = par$beta
   DD = par$DD
   D.inv = solve(DD)
   Sigma = par$Sigma
   phi = par$phi
   ga = par$ga
   Yim = as.list(num.mis.sub)
   TrZ = matrix(0, ncol=N*q, nrow=n)
   Lambda = SigCi = Ytil.hat = bi = as.list(N)
   Del.e = Del.b = numeric(N)
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
     Lambda[[i]] = TrZ[idx, ((i-1)*q+1): (i*q)] %*% DD %*% t(TrZ[idx, ((i-1)*q+1): (i*q)])+ SigCi[[i]]
   }
   Ytilde = Y.na - MU + Xtilde %*% beta + TrZ %*% as.vector(t(b))
   Y.cent = Ytilde - Xtilde %*% beta
   ycent.na = Y.na - MU
# EM-step:
   if(num.obs.sub != 0){
   for(i in 1: num.obs.sub){
     if(obs.sub[i] == 1) idx = 1: cumsum.oi[i]
     if(obs.sub[i] != 1) idx = (cumsum.oi[obs.sub[i]-1]+1): cumsum.oi[obs.sub[i]]
         Ytil.hat[[obs.sub[i]]] = Ytilde[idx]
         bi[[obs.sub[i]]] = DD %*% t(TrZ[idx, ((obs.sub[i]-1)*q+1): (obs.sub[i]*q)]) %*% solve(Lambda[[obs.sub[i]]]) %*% Y.cent[idx]
         Del.e[obs.sub[i]] = t(ycent.na[idx]) %*% solve(SigCi[[obs.sub[i]]]) %*% ycent.na[idx]
         Del.b[obs.sub[i]] = t(b[obs.sub[i], ]) %*% D.inv %*% b[obs.sub[i], ]
   }}
   if(num.mis.sub != 0){
   for(i in 1: num.mis.sub){
     if(mis.sub[i] == 1) idx = 1: cumsum.oi[i]
     if(mis.sub[i] != 1) idx = (cumsum.oi[mis.sub[i]-1]+1): cumsum.oi[mis.sub[i]]
     Xim = matrix(Xtilde[idx, ][which(mis.ind[[i]]==T), ], ncol=p)
     Xio = matrix(Xtilde[idx, ][-which(mis.ind[[i]]==T), ], ncol=p)
     Zio = matrix(TrZ[idx, ((mis.sub[i]-1)*q+1): (mis.sub[i]*q)][-which(mis.ind[[i]]==T), ], ncol=q)
     Lam.oo.inv = solve(Lambda[[mis.sub[i]]][-which(mis.ind[[i]]==T), -which(mis.ind[[i]]==T)])
     Yim[[i]] = c(Xim %*% beta + Lambda[[mis.sub[i]]][which(mis.ind[[i]]==T), -which(mis.ind[[i]]==T)] %*% Lam.oo.inv %*% Y.cent[idx][-which(mis.ind[[i]]==T)])
     Oi = matrix(diag(oi[mis.sub[i]]*r)[-which(mis.ind[[i]]==T), ], ncol=oi[mis.sub[i]]*r)
     Mi = matrix(diag(oi[mis.sub[i]]*r)[which(mis.ind[[i]]==T), ], ncol=oi[mis.sub[i]]*r)
     Ytil.hat[[mis.sub[i]]] = t(Oi) %*% Ytilde[idx][-which(mis.ind[[i]]==T)] + t(Mi) %*% Yim[[i]]
     bi[[mis.sub[i]]] = DD %*% t(Zio) %*% Lam.oo.inv %*% Y.cent[idx][-which(mis.ind[[i]]==T)]
     Del.e[mis.sub[i]] = t(ycent.na[idx][-which(mis.ind[[i]]==T)]) %*% solve(SigCi[[mis.sub[i]]][-which(mis.ind[[i]]==T), -which(mis.ind[[i]]==T)]) %*% ycent.na[idx][-which(mis.ind[[i]]==T)]
     Del.b[mis.sub[i]] = t(b[mis.sub[i], ]) %*% D.inv %*% b[mis.sub[i], ]
   }}
   Y.hat = bb = NULL
   for(i in 1: N){
     Y.hat = c(Y.hat, Ytil.hat[[i]])
     bb = c(bb, bi[[i]])
   }
   b = t(matrix(bb, ncol=N))
   Yimp = Y.hat + MU - Xtilde %*% beta - TrZ %*% bb
   return(list(b=b, Yhat=MU, Yimp=Yimp, Del.e=Del.e, Del.b=Del.b))
}
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         