library(nlme)
library(mvtnorm)
library(MCMCpack)
library(tmvtnorm)

### Diagnosis convergence: (max.bins) Maximum number of bins ###
rdc=function(x,n) x[ceiling(n/2+1):n]
cent.sum=function(X) (t(X)-colMeans(X))%*%t(t(X)-colMeans(X))

MPSRF.multi=function(M, max.bins=50, start=0)
{
   N=apply(M, 3, nrow)[1]
   m=length(apply(M, 3, nrow))
   bin.width=floor((N-start)/max.bins)
   npar=apply(M, 3, ncol)[1]
   Rphat=numeric(0)
   for(i in 1: max.bins)
   {
     sub.M=apply(apply(M, c(3, 2), rdc, n=start+i*bin.width), c(3, 2), cbind)
     n=apply(sub.M, 3, nrow)[1]
     W=matrix(rowSums(apply(sub.M, 3, cent.sum)), npar, npar)/(m*(n-1))
     theta.bdd=rowSums(apply(sub.M, 3, colSums))/(m*n)
     B.over.n=(apply(sub.M, 3, colMeans)-theta.bdd)%*%t(apply(sub.M, 3, colMeans)-theta.bdd)/(m-1)
     lambda1=Re(eigen(solve(W)%*%B.over.n)$values[1])
     Rphat[i]=(n-1)/n+((m+1)/m)*lambda1
   }
   sum.id = numeric(max.bins)
   for(k in 1: max.bins) sum.id[k] = sum((Rphat <= 1.1)[k:max.bins])
   a = min(which(sum.id == (max.bins-(1:max.bins)+1)))
   if(a=='Inf'){ conv.step = 'NA'
   } else conv.step = seq(start, nrow(M), bin.width)[-1][a]
   return(list(Rphat=Rphat, conv.step=conv.step))
}

# vech:
vech.posi=function(dim) cbind(rep(1:dim, 1:dim), unlist(mapply(':', 1, 1:dim)))

DEC = function(phi, ga, Ti) phi^((abs(outer(Ti, Ti, '-')))^ga)

### Specify the values of hyperparameters ###
hyper.par = function(par, FI)
{
  beta0 = par$beta 
  J0 = 9 * solve(FI$Ibeta)  
  G0 = par$DD
  H0 = par$Sigma
  q = ncol(G0)
  r = ncol(H0)
  d0 = q +2
  s0 = r + 2
  pg.mle = c(par$phi, par$ga)
  return(list(par.mle = par, beta0 = beta0, J0 = J0, d0 = d0, G0 = G0, s0 = s0, H0 = H0, pg.mle = pg.mle))
}

# ML Estimation
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

### MCMC for MNLMM-CM ###
MCMC.MNLMMCM = function(chain, ITER, hyper, typeCi=c('UNC','AR1','DEC'), nlag=10, per=100, method=c('pseudo','raw'))
{
   begin=proc.time()[1]
# data
   N = length(unique(vecData$Subject))
   y.na = vecData$Resp
   n = length(y.na)
   na.ind = which(is.na(y.na))
   y = y.na
   y[na.ind] = 999

   Data.o = vecData[which(vecData$Censor==0),]
   Data.c = vecData[which(vecData$Censor==1),]
   Data.m = vecData[na.ind,]
   Data.p = vecData[-na.ind, ]
   cen.ind = which(Data.p$Censor == 1)
   yp = y[-na.ind]
   np = length(yp)
   yo = yp[-cen.ind]
   no = length(yo)

   si = oi = pii = numeric(N)
   for(i in 1: N) si[i] = length(unique(Data$Time[Data$Subject == i]))
   ni = si * r
   cumsum.ni = cumsum(ni)
   for(i in 1: N) oi[i] = length(Data.o$Time[Data.o$Subject == i])
   cumsum.oi = cumsum(oi)
   for(i in 1: N) pii[i] = length(Data.p$Time[Data.p$Subject == i])
   cumsum.pi = cumsum(pii)
   cen.subj = unique(Data.c$Subject)
   cen.subj = as.numeric(levels(cen.subj))[cen.subj]
   na.subj = unique(Data.m$Subject)
   na.subj = as.numeric(levels(na.subj))[na.subj]
   obs.subj = c(1:N)[-sort(unique(c(cen.subj, na.subj)))]
   Nc = length(cen.subj)
   ci = numeric(Nc)
   for(i in 1: Nc) ci[i] = length(Data.c$Time[Data.c$Subject == cen.subj[i]])
   cumsum.nc = cumsum(ci)
   num.cen = length(cen.ind)
   num.na = length(na.ind)

# hyper parameters:
   beta0 = hyper$beta0 
   J0 = hyper$J0 
   J0.inv = solve(J0)
   d0 = hyper$d0 
   G0 = hyper$G0 
   s0 = hyper$s0 
   H0 = hyper$H0
   pg.mle = hyper$pg.mle
   p = length(beta0)
   q = ncol(G0)
   r = nrow(H0)
   m = p + q*(q+1)/2 + r*(r+1)/2 + 2
   typeCi = typeCi[1]
   method = method[1]
   vech.DD = vech.posi((q))
   vech.Sig = vech.posi(r)

   A = diag(p)
   B = matrix(c(1,rep(0, 7), 1, 0), ncol=q)
   TrZ = matrix(0, ncol=N*q, nrow=n)
   TP = diag(n)[-na.ind, ]
   TM = diag(n)[na.ind, ]
   TO = diag(np)[-cen.ind, ]
   TC = diag(np)[cen.ind, ]
   if(num.cen == 1) TC = t(TC)

# save posterior samples:
   Theta = array(NA, dim=c(ITER, (m+1), chain))
   TB = array(NA, dim=c(ITER, N*q, chain))
   TYimp = array(NA, dim=c(ITER, num.na, chain))
   TYcen = array(NA, dim=c(ITER, num.cen, chain))
   TYfit = array(NA, dim=c(ITER, n, chain))
   Init = NULL

   for(cha in 1: chain){
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
     if(i == 1) idx = 1: cumsum.ni[i]
     else idx = (cumsum.ni[i-1]+1): cumsum.ni[i]
     TrZ[idx, ((i-1)*q+1): (i*q)] = rbind(dmu.i$dmu1 %*% B, dmu.i$dmu2 %*% B)
   }
   Ytilde = y - MU + Xtilde %*% beta + TrZ %*% as.vector(t(b))
   Xp.tilde = Xtilde[-na.ind, ]
   Xm.tilde = Xtilde[na.ind, ]
   Xc.tilde = Xp.tilde[cen.ind, ]
   Xo.tilde = Xp.tilde[-cen.ind, ]
   Yp.tilde = Ytilde[-na.ind]
   Yo.tilde = Yp.tilde[-cen.ind]
   Utilde = Yp.tilde[cen.ind]

 # observed log-likelihood:
   TLam = TOme = matrix(0, ncol=n, nrow=n)
   TOme[1: cumsum.ni[1], 1: cumsum.ni[1]] = kronecker(Sigma, DEC(phi,ga, Data$Time[Data$Subject == 1]))
   Zi = matrix(TrZ[1: cumsum.ni[1], 1: q], ncol=q)
   TLam[1: cumsum.ni[1], 1: cumsum.ni[1]] = as.matrix(TrZ[1: cumsum.ni[1], 1:q]) %*% DD %*% t(TrZ[1: cumsum.ni[1], 1:q]) + TOme[1: cumsum.ni[1], 1: cumsum.ni[1]]
   for(i in 2: N){
    Zi = matrix(TrZ[(cumsum.ni[i-1]+1): cumsum.ni[i], ((i-1)*q+1): (i*q)], ncol=q)
    TOme[(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]] = kronecker(Sigma, DEC(phi,ga, Data$Time[Data$Subject == i]))
    TLam[(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]] = Zi %*% DD %*% t(Zi) + TOme[(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]]    
   }
   TLam.pp = TLam[-na.ind, -na.ind]
   TLam.mp = TLam[na.ind, -na.ind]
   TLam.mm = TLam[na.ind, na.ind]
   TLam.oo = TLam.pp[-cen.ind, -cen.ind]
   TLam.co = TLam.pp[cen.ind, -cen.ind]
   TLam.cc = TLam.pp[cen.ind, cen.ind]
   if(num.cen == 1){
     TLam.co = t(TLam.co); TLam.cc = as.matrix(TLam.cc)
   }
   if(num.na == 1) TLam.mp = t(TLam.mp) 
   TLam.oo.inv = solve(TLam.oo)
   yo.cent.tilde = Yo.tilde - Xo.tilde %*% beta
   mu.co = Xc.tilde %*% beta + TLam.co %*% TLam.oo.inv %*% yo.cent.tilde
   Sig.cc.o = TLam.cc - TLam.co %*% TLam.oo.inv %*% t(TLam.co)
   if(ci[1]==1){ log.cdf = log(pnorm(Utilde[1:cumsum.nc[1]], mean=mu.co[1:cumsum.nc[1]], sd=sqrt(Sig.cc.o[1:cumsum.nc[1], 1:cumsum.nc[1]])))
   } else log.cdf = log(pmvnorm(lower=rep(-Inf, ci[1]), upper=Utilde[1:cumsum.nc[1]], mean=c(mu.co[1:cumsum.nc[1]]), sigma=Sig.cc.o[1:cumsum.nc[1], 1:cumsum.nc[1]]))[1]
   if(Nc != 1){
   for(i in 2: Nc){
   if(ci[i]==1){ log.cdf = log.cdf + log(pnorm(Utilde[(cumsum.nc[i-1]+1): cumsum.nc[i]], mean=mu.co[(cumsum.nc[i-1]+1): cumsum.nc[i]], sd=sqrt(Sig.cc.o[(cumsum.nc[i-1]+1): cumsum.nc[i], (cumsum.nc[i-1]+1): cumsum.nc[i]])))
   } else log.cdf = log.cdf + log(pmvnorm(lower=rep(-Inf, ci[i]), upper=Utilde[(cumsum.nc[i-1]+1): cumsum.nc[i]], mean=c(mu.co[(cumsum.nc[i-1]+1): cumsum.nc[i]]), sigma=Sig.cc.o[(cumsum.nc[i-1]+1): cumsum.nc[i], (cumsum.nc[i-1]+1): cumsum.nc[i]]))[1]
   }}
   log.det.Sig.inv = log(det(TLam.oo.inv[1:cumsum.oi[1], 1:cumsum.oi[1]]))
   for(i in 2: N) log.det.Sig.inv = log.det.Sig.inv + log(det(as.matrix(TLam.oo.inv[(cumsum.oi[i-1]+1): cumsum.oi[i], (cumsum.oi[i-1]+1): cumsum.oi[i]])))
   loglik = log.cdf -0.5*log(2*pi)*no + 0.5*log.det.Sig.inv - 0.5*t(yo.cent.tilde) %*% TLam.oo.inv %*% yo.cent.tilde
   cat(rep('=', 25), 'MNLMM with ', typeCi, ' errors; ', 'censoring = ', num.cen/n*100, '%; missing = ', num.na/n*100, '%', rep('=', 25), sep = '', '\n')
   cat(paste(rep('-', 20), sep = '', collapse = ''), 'MCMC with ', method, ' method: Chain = ', cha, '\t iter = 0, loglik = ', loglik, paste(rep('-', 20), sep = '', collapse = ''), '\n')
   cat('Initial values: beta = ',  round(beta, 2), ',\t D = ', round(DD[vech.DD], 2), '\t Sigma = ', round(Sigma[vech.Sig], 3), ',\t (phi,ga) = ', round(c(phi,ga), 2), '\n')
   
   for(iter in 1: ITER)
   {
# Imputation Step:
 # Yc #
     yc.tilde.hat = c(rtmvnorm(1, mean=c(mu.co), sigma=Sig.cc.o, lower=rep(-Inf, num.cen), upper=Utilde, algorithm="gibbs", burn.in.samples=10))  # Should we choose different burn.in.sample size? 
 # Ym #
     yp.tilde.hat = t(TO) %*% Yo.tilde + t(TC) %*% yc.tilde.hat 
     yp.tilde.cent = yp.tilde.hat - Xp.tilde %*% beta
     TLam.pp.inv = solve(TLam.pp)
     mu.mp = Xm.tilde %*% beta + TLam.mp %*% TLam.pp.inv %*% yp.tilde.cent
     Sig.mm.p = TLam.mm - TLam.mp %*% TLam.pp.inv %*% t(TLam.mp)
     ym.tilde.hat = c(rmvnorm(1, mu.mp, Sig.mm.p))
     y.tilde.hat = t(TP) %*% yp.tilde.hat + t(TM) %*% ym.tilde.hat
 # b #
     TOme.pp.inv = solve(TOme[-na.ind, -na.ind])
     Sb = solve(t(TrZ[-na.ind, ]) %*% TOme.pp.inv %*% TrZ[-na.ind, ] + solve(kronecker(diag(N), DD)))
     mu.b = Sb %*% t(TrZ[-na.ind, ]) %*% TOme.pp.inv %*% yp.tilde.cent
     bb = c(rmvnorm(1, mu.b, Sb))
     b = matrix(bb, nrow=N, byrow=T)

# Posterior Step:
 # beta (pseudo data) #
     Sig.beta = solve(t(Xp.tilde) %*% TLam.pp.inv %*% Xp.tilde + J0.inv)
     mu.beta = Sig.beta %*% (t(Xp.tilde) %*% TLam.pp.inv %*% yp.tilde.hat + J0.inv %*% beta0)
     beta = c(rmvnorm(1, mu.beta, Sig.beta))
 # DD #
     DD = riwish(N+d0, t(b)%*%b+G0)
 # Sigma #
     if(method=='raw'){ 
     MU = NULL
     for(i in 1: N){
       eta.i = A %*% beta + B %*% b[i, ]
       mu.i = mu.fn(eta.i, Data$Time[Data$Subject == i])
       MU = c(MU, mu.i$mu1, mu.i$mu2)
     }
     Y.hat = y.tilde.hat + MU - Xtilde %*% beta - TrZ %*% as.vector(t(b))
     Y.cent = Y.hat - MU
     e.c = NULL
     for(i in 2: N) e.c = c(e.c, rep(0, n), Y.cent[(cumsum.ni[i-1]+1): cumsum.ni[i]])
     e.c = matrix(c(Y.cent[1: cumsum.ni[1]], e.c), ncol = N)
     }
     if(method=='pseudo'){
     Y.tilde.cent = y.tilde.hat - Xtilde %*% beta - TrZ %*%as.vector(t(b))
     e.c = NULL
     for(i in 2: N) e.c = c(e.c, rep(0, n), Y.tilde.cent[(cumsum.ni[i-1]+1): cumsum.ni[i]])
     e.c = matrix(c(Y.tilde.cent[1: cumsum.ni[1]], e.c), ncol = N)
     }
     TE = e.c %*% t(e.c)
     Psi = diag(r)
     if(r == 1){
       Ce = 0
       Ce = sum(solve(DEC(phi, ga, unique(Data$Time[Data$Subject == 1]))) * TE[1: cumsum.ni[1], ][1: si[1], 1: si[1]])
       for(i in 2: N) Ce = Ce + sum(solve(DEC(phi, ga, unique(Data$Time[Data$Subject == i])))*TE[(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]][1: si[i], 1: si[i]])
       Psi = Ce
     } else{
     for(j in 1: r)for(l in 1: r)
     {
       Ce=0
       Ce=sum(solve(DEC(phi, ga, Data$Time[Data$Subject == 1])) * TE[1: cumsum.ni[1], 1: cumsum.ni[1]][((j-1)*si[1]+1): (j*si[1]), ((l-1)*si[1]+1): (l*si[1])])
       for(i in 2: N) Ce = Ce + sum(solve(DEC(phi, ga, Data$Time[Data$Subject == i]))*TE[(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]][((j-1)*si[i]+1): (j*si[i]), ((l-1)*si[i]+1): (l*si[i])])
       Psi[j, l] = Ce
     }}
     Sigma = riwish(sum(si)+s0, Psi+H0)
 # phi and ga #
     if(typeCi == 'UNC') phi = 1e-6; ga = 1
     if(typeCi == 'AR1'){
        ga = 1
        if(method=='raw') yp.cent = Y.cent[-na.ind]
        if(method=='pseudo') yp.cent = yp.tilde.cent
        phi = MH.phi.nonli(phi, DD, Sigma, yp.cent, pg.mle, TrZ, method=method)
     }
     if(typeCi == 'DEC'){
        if(method=='raw') yp.cent = Y.cent[-na.ind]
        if(method=='pseudo') yp.cent = yp.tilde.cent
        Phi = MH.DEC.nonli(phi, ga, DD, Sigma, yp.cent, pg.mle, TrZ, method=method)
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
       if(i == 1) idx = 1: cumsum.ni[i]
       else idx = (cumsum.ni[i-1]+1): cumsum.ni[i]
       TrZ[idx, ((i-1)*q+1): (i*q)] = rbind(dmu.i$dmu1 %*% B, dmu.i$dmu2 %*% B)
     }
     Ytilde = y - MU + Xtilde %*% beta + TrZ %*% as.vector(t(b))
     Xp.tilde = Xtilde[-na.ind, ]
     Xm.tilde = Xtilde[na.ind, ]
     Xc.tilde = Xp.tilde[cen.ind, ]
     Xo.tilde = Xp.tilde[-cen.ind, ]
     Yp.tilde = Ytilde[-na.ind]
     Yo.tilde = Yp.tilde[-cen.ind]
     Utilde = Yp.tilde[cen.ind]
     Y.hat = y.tilde.hat + MU - Xtilde %*% beta - TrZ %*% as.vector(t(b))

 # observed log-likelihood:
     TLam = TOme = matrix(0, ncol=n, nrow=n)
     TOme[1: cumsum.ni[1], 1: cumsum.ni[1]] = kronecker(Sigma, DEC(phi,ga, Data$Time[Data$Subject == 1]))
     Zi = matrix(TrZ[1: cumsum.ni[1], 1: q], ncol=q)
     TLam[1: cumsum.ni[1], 1: cumsum.ni[1]] = as.matrix(TrZ[1: cumsum.ni[1], 1:q]) %*% DD %*% t(TrZ[1: cumsum.ni[1], 1:q]) + TOme[1: cumsum.ni[1], 1: cumsum.ni[1]]
     for(i in 2: N){
      Zi = matrix(TrZ[(cumsum.ni[i-1]+1): cumsum.ni[i], ((i-1)*q+1): (i*q)], ncol=q)
      TOme[(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]] = kronecker(Sigma, DEC(phi,ga, Data$Time[Data$Subject == i]))
      TLam[(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]] = Zi %*% DD %*% t(Zi) + TOme[(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]]    
     }
     TLam.pp = TLam[-na.ind, -na.ind]
     TLam.mp = TLam[na.ind, -na.ind]
     TLam.mm = TLam[na.ind, na.ind]
     TLam.oo = TLam.pp[-cen.ind, -cen.ind]
     TLam.co = TLam.pp[cen.ind, -cen.ind]
     TLam.cc = TLam.pp[cen.ind, cen.ind]
     if(num.cen == 1){
       TLam.co = t(TLam.co); TLam.cc = as.matrix(TLam.cc)
     }
     if(num.na == 1) TLam.mp = t(TLam.mp) 
     TLam.oo.inv = solve(TLam.oo)
     yo.cent.tilde = Yo.tilde - Xo.tilde %*% beta
     mu.co = Xc.tilde %*% beta + TLam.co %*% TLam.oo.inv %*% yo.cent.tilde
     Sig.cc.o = TLam.cc - TLam.co %*% TLam.oo.inv %*% t(TLam.co)
     if(ci[1]==1){ log.cdf = log(pnorm(Utilde[1:cumsum.nc[1]], mean=mu.co[1:cumsum.nc[1]], sd=sqrt(Sig.cc.o[1:cumsum.nc[1], 1:cumsum.nc[1]])))
     } else log.cdf = log(pmvnorm(lower=rep(-Inf, ci[1]), upper=Utilde[1:cumsum.nc[1]], mean=c(mu.co[1:cumsum.nc[1]]), sigma=Sig.cc.o[1:cumsum.nc[1], 1:cumsum.nc[1]]))[1]
     if(Nc != 1){
     for(i in 2: Nc){
     if(ci[i]==1){ log.cdf = log.cdf + log(pnorm(Utilde[(cumsum.nc[i-1]+1): cumsum.nc[i]], mean=mu.co[(cumsum.nc[i-1]+1): cumsum.nc[i]], sd=sqrt(Sig.cc.o[(cumsum.nc[i-1]+1): cumsum.nc[i], (cumsum.nc[i-1]+1): cumsum.nc[i]])))
     } else log.cdf = log.cdf + log(pmvnorm(lower=rep(-Inf, ci[i]), upper=Utilde[(cumsum.nc[i-1]+1): cumsum.nc[i]], mean=c(mu.co[(cumsum.nc[i-1]+1): cumsum.nc[i]]), sigma=Sig.cc.o[(cumsum.nc[i-1]+1): cumsum.nc[i], (cumsum.nc[i-1]+1): cumsum.nc[i]]))[1]
     }}
     log.det.Sig.inv = log(det(TLam.oo.inv[1:cumsum.oi[1], 1:cumsum.oi[1]]))
     for(i in 2: N) log.det.Sig.inv = log.det.Sig.inv + log(det(as.matrix(TLam.oo.inv[(cumsum.oi[i-1]+1): cumsum.oi[i], (cumsum.oi[i-1]+1): cumsum.oi[i]])))
     loglik = log.cdf -0.5*log(2*pi)*no + 0.5*log.det.Sig.inv - 0.5*t(yo.cent.tilde) %*% TLam.oo.inv %*% yo.cent.tilde

# save posterior samples:
     if(r==1 & q==1) pos.samp = c(beta, DD, Sigma, phi, ga)
     else pos.samp = c(beta, DD[vech.DD], Sigma[vech.Sig], phi, ga)
     if(iter%%per == 0) cat('iter = ', iter, ',\t loglik = ', loglik, ',\t pos.samples = ', round(pos.samp, 3), '\n')
     Theta[iter, , cha] = c(pos.samp, loglik)
     TB[iter, , cha] = bb
     TYimp[iter, , cha] = Y.hat[na.ind]
     TYcen[iter, , cha] = Y.hat[-na.ind][cen.ind]
     TYfit[iter, , cha] = MU
#     if(PATH != NULL){
#      write(c(pos.samp, loglik.iter[k, iter]), paste(PATH,'theta.txt',sep=""), ncol=length(init)+1, append=T)
#      write(as.vector(t(b)), paste(PATH,'b.txt',sep=""), ncol=(q*N), append=T)
#     }
   }}
   end = proc.time()[1]
   cat('It took', end - begin, 'seconds.\n')
   data.inf = list(obs.subj=obs.subj, cen.subj=cen.subj, na.subj=na.subj, cen.ind=which(vecData$Censor==1), na.ind=na.ind)

# summary (after convergence):
   if(typeCi == 'UNC') m1 = m - 2
   if(typeCi == 'AR1') m1 = m - 1
   if(typeCi == 'DEC') m1 = m
   CHECK = MPSRF.multi(Theta[, 1: m1, ], max.bins=50, start=0)   
   burn.in = CHECK$conv.step 
   if(burn.in == 'NA') burn.in = ITER/2
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
   phi.hat = Phi.hat[1] 
   ga.hat = Phi.hat[2]
   pos.sd = apply(Theta[-c(1:burn.in),-(m+1),][cho,,], 2, sd)
   pos.CI = apply(Theta[-c(1:burn.in),-(m+1),][cho,,], 2, quantile, prob=c(0.025, 0.975, 0.5))
# ym imputation
   Ym.hat = apply(TYimp[-c(1:burn.in),,][cho,,], 2, mean)
   Ym.sd = apply(TYimp[-c(1:burn.in),,][cho,,], 2, sd)
   Ym.CI = apply(TYimp[-c(1:burn.in),,][cho,,], 2, quantile, prob=c(0.025, 0.975, 0.5))
# yc imputation
   Yc.hat = apply(TYcen[-c(1:burn.in),,][cho,,], 2, mean)
   Yc.sd = apply(TYcen[-c(1:burn.in),,][cho,,], 2, sd)
   Yc.CI = apply(TYcen[-c(1:burn.in),,][cho,,], 2, quantile, prob=c(0.025, 0.975, 0.5))
   theta.out = rbind(theta.hat, pos.sd, pos.CI)
   Ym.out = rbind(Ym.hat, Ym.sd, Ym.CI)
   Yc.out = rbind(Yc.hat, Yc.sd, Yc.CI)
# fitted responses 
   Yfit = apply(TYfit[-c(1:burn.in),,][cho,,], 2, mean)
   Yfit.sd = apply(TYfit[-c(1:burn.in),,][cho,,], 2, sd)
   Yfit.CI = apply(TYfit[-c(1:burn.in),,][cho,,], 2, quantile, prob=c(0.025, 0.975, 0.5))
   yfit.out = rbind(Yfit, Yfit.sd, Yfit.CI)
   rownames(theta.out) = rownames(Ym.out) = rownames(Yc.out) = rownames(yfit.out) = c('Mean','SD','2.5%','97.5%','Median')  
# random effects estimation
   b.hat = matrix(apply(TB[-c(1:burn.in),,][cho,,], 2, mean), nrow=N, byrow=T)
   pos.summ = list(theta.out=theta.out, Ym.out=Ym.out, Yc.out=Yc.out, yfit.out=yfit.out, b.hat=b.hat)

# model selection criteria:
   D.bar = -2 * mean(Theta[-c(1:burn.in), (m+1),][cho,])
 # DIC #
   MU = Xtilde = NULL
   for(i in 1: N){
     eta.i = A %*% beta.hat + B %*% b.hat[i, ]
     mu.i = mu.fn(eta.i, Data$Time[Data$Subject == i])
     MU = c(MU, mu.i$mu1, mu.i$mu2)
     dmu.i = dmu(eta.i, Data$Time[Data$Subject == i])
     Xtilde = rbind(Xtilde, dmu.i$dmu1 %*% A, dmu.i$dmu2 %*% A)
     if(i == 1) idx = 1: cumsum.ni[i]
     else idx = (cumsum.ni[i-1]+1): cumsum.ni[i]
     TrZ[idx, ((i-1)*q+1): (i*q)] = rbind(dmu.i$dmu1 %*% B, dmu.i$dmu2 %*% B)
   }
   Ytilde = y - MU + Xtilde %*% beta.hat + TrZ %*% as.vector(t(b.hat))
   Xp.tilde = Xtilde[-na.ind, ]
#   Xm.tilde = Xtilde[na.ind, ]
   Xc.tilde = Xp.tilde[cen.ind, ]
   Xo.tilde = Xp.tilde[-cen.ind, ]
   Yp.tilde = Ytilde[-na.ind]
   Yo.tilde = Yp.tilde[-cen.ind]
   Utilde = Yp.tilde[cen.ind]

 # observed log-likelihood:
   TLam = TOme = matrix(0, ncol=n, nrow=n)
   TOme[1: cumsum.ni[1], 1: cumsum.ni[1]] = kronecker(Sig.hat, DEC(phi.hat, ga.hat, Data$Time[Data$Subject == 1]))
   Zi = matrix(TrZ[1: cumsum.ni[1], 1: q], ncol=q)
   TLam[1: cumsum.ni[1], 1: cumsum.ni[1]] = as.matrix(TrZ[1: cumsum.ni[1], 1:q]) %*% D.hat %*% t(TrZ[1: cumsum.ni[1], 1:q]) + TOme[1: cumsum.ni[1], 1: cumsum.ni[1]]
   for(i in 2: N){
    Zi = matrix(TrZ[(cumsum.ni[i-1]+1): cumsum.ni[i], ((i-1)*q+1): (i*q)], ncol=q)
    TOme[(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]] = kronecker(Sig.hat, DEC(phi.hat, ga.hat, Data$Time[Data$Subject == i]))
    TLam[(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]] = Zi %*% D.hat %*% t(Zi) + TOme[(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]]    
   }
   TLam.pp = TLam[-na.ind, -na.ind]
   TLam.mp = TLam[na.ind, -na.ind]
   TLam.mm = TLam[na.ind, na.ind]
   TLam.oo = TLam.pp[-cen.ind, -cen.ind]
   TLam.co = TLam.pp[cen.ind, -cen.ind]
   TLam.cc = TLam.pp[cen.ind, cen.ind]
   if(num.cen == 1){
     TLam.co = t(TLam.co); TLam.cc = as.matrix(TLam.cc)
   }
   if(num.na == 1) TLam.mp = t(TLam.mp) 
   TLam.oo.inv = solve(TLam.oo)
   yo.cent.tilde = Yo.tilde - Xo.tilde %*% beta.hat
   mu.co = Xc.tilde %*% beta.hat + TLam.co %*% TLam.oo.inv %*% yo.cent.tilde
   Sig.cc.o = TLam.cc - TLam.co %*% TLam.oo.inv %*% t(TLam.co)
   if(ci[1]==1){ log.cdf = log(pnorm(Utilde[1:cumsum.nc[1]], mean=mu.co[1:cumsum.nc[1]], sd=sqrt(Sig.cc.o[1:cumsum.nc[1], 1:cumsum.nc[1]])))
   } else log.cdf = log(pmvnorm(lower=rep(-Inf, ci[1]), upper=Utilde[1:cumsum.nc[1]], mean=c(mu.co[1:cumsum.nc[1]]), sigma=Sig.cc.o[1:cumsum.nc[1], 1:cumsum.nc[1]]))[1]
   if(Nc != 1){
   for(i in 2: Nc){
   if(ci[i]==1){ log.cdf = log.cdf + log(pnorm(Utilde[(cumsum.nc[i-1]+1): cumsum.nc[i]], mean=mu.co[(cumsum.nc[i-1]+1): cumsum.nc[i]], sd=sqrt(Sig.cc.o[(cumsum.nc[i-1]+1): cumsum.nc[i], (cumsum.nc[i-1]+1): cumsum.nc[i]])))
   } else log.cdf = log.cdf + log(pmvnorm(lower=rep(-Inf, ci[i]), upper=Utilde[(cumsum.nc[i-1]+1): cumsum.nc[i]], mean=c(mu.co[(cumsum.nc[i-1]+1): cumsum.nc[i]]), sigma=Sig.cc.o[(cumsum.nc[i-1]+1): cumsum.nc[i], (cumsum.nc[i-1]+1): cumsum.nc[i]]))[1]
   }}
   log.det.Sig.inv = log(det(TLam.oo.inv[1:cumsum.oi[1], 1:cumsum.oi[1]]))
   for(i in 2: N) log.det.Sig.inv = log.det.Sig.inv + log(det(as.matrix(TLam.oo.inv[(cumsum.oi[i-1]+1): cumsum.oi[i], (cumsum.oi[i-1]+1): cumsum.oi[i]])))
   loglik.hat = log.cdf -0.5*log(2*pi)*no + 0.5*log.det.Sig.inv - 0.5*t(yo.cent.tilde) %*% TLam.oo.inv %*% yo.cent.tilde

   D.theta.bar = -2 * loglik.hat
   DIC = 2*D.bar - D.theta.bar
   EAIC = D.bar + 2* m1
   EBIC = D.bar + m1 * log(N)

 # CPO & LPML #
   new.ITER = (ITER - burn.in)/nlag
   ind.lik = array(0, dim=c(new.ITER, N, chain))
   theta.new = Theta[-c(1:burn.in),-(m+1),][cho,,]
   DD.l = diag(q); Sig.l = diag(r)
   TB.new = TB[-c(1:burn.in),,][cho,,]
   for(k in 1: chain){
   for(iter in 1: new.ITER){
     beta.l = theta.new[iter, 1:p, k]
     DD.l[vech.DD] = DD.l[inv.vech.DD] = theta.new[iter, -c(1:p), k][1:r1]
     Sig.l[vech.Sig] = Sig.l[inv.vech.Sig] = theta.new[iter, -c(1:(p+r1)), k][1:r2]
     phi.l = rev(theta.new[iter,,k])[2] 
     ga.l = rev(theta.new[iter,,k])[1]
     b.l = matrix(TB.new[iter,,k], nrow=N, byrow=T) 
     ind.lik[iter,,k] = mnlmmcm.ind.loglik(beta.l, DD.l, Sig.l, phi.l, ga.l, b.l)
   }}
   CPO = 1 / apply(1/exp(ind.lik), 2, mean)
   LPML = sum(log(CPO))
   KL = -log(CPO) + apply(ind.lik, 2, mean)
   model.inf = list(m=m1, DIC=DIC, EAIC=EAIC, EBIC=EBIC, CPO=CPO, LPML=LPML, ind.lnL=ind.lik, KL=KL)
   return(list(model.inf=model.inf, pos.inf=pos.summ, Theta = Theta, TB = TB, TYimp=TYimp, TYcen=TYcen, Init=Init, data.inf=data.inf, burn.in=burn.in, size=size*chain, CPUtime=end-begin))
}

# individual log-likelihood for mlmm-cm
mnlmmcm.ind.loglik = function(beta, DD, Sig, phi, ga, b) 
{
   N = length(unique(vecData$Subject))
   y.na = vecData$Resp
   n = length(y.na)
   na.ind = which(is.na(y.na))
   y = y.na
   y[na.ind] = 999

   Data.o = vecData[which(vecData$Censor==0),]
   Data.c = vecData[which(vecData$Censor==1),]
   Data.p = vecData[-na.ind, ]
   cen.ind = which(Data.p$Censor == 1)
   yp = y[-na.ind]
   np = length(yp)
   yo = yp[-cen.ind]
   no = length(yo)

   si = oi = numeric(N)
   for(i in 1: N) si[i] = length(unique(Data$Time[Data$Subject == i]))
   ni = si * r
   cumsum.ni = cumsum(ni)
   for(i in 1: N) oi[i] = length(Data.o$Time[Data.o$Subject == i])
   cumsum.oi = cumsum(oi)
   cen.subj = unique(Data.c$Subject)
   cen.subj = as.numeric(levels(cen.subj))[cen.subj]
   Nc = length(cen.subj)
   ci = numeric(Nc)
   for(i in 1: Nc) ci[i] = length(Data.c$Time[Data.c$Subject == cen.subj[i]])
   cumsum.nc = cumsum(ci)
   num.cen = length(cen.ind)
   num.na = length(na.ind)

   p = length(beta)
   q = ncol(DD)
   A = diag(p)
   B = matrix(c(1,rep(0, 7), 1, 0), ncol=q)
   TrZ = matrix(0, ncol=N*q, nrow=n)
   MU = Xtilde = NULL
   for(i in 1: N){
     eta.i = A %*% beta + B %*% b[i, ]
     mu.i = mu.fn(eta.i, Data$Time[Data$Subject == i])
     MU = c(MU, mu.i$mu1, mu.i$mu2)
     dmu.i = dmu(eta.i, Data$Time[Data$Subject == i])
     Xtilde = rbind(Xtilde, dmu.i$dmu1 %*% A, dmu.i$dmu2 %*% A)
     if(i == 1) idx = 1: cumsum.ni[i]
     else idx = (cumsum.ni[i-1]+1): cumsum.ni[i]
     TrZ[idx, ((i-1)*q+1): (i*q)] = rbind(dmu.i$dmu1 %*% B, dmu.i$dmu2 %*% B)
   }
   Ytilde = y - MU + Xtilde %*% beta + TrZ %*% as.vector(t(b))
   Xp.tilde = Xtilde[-na.ind, ]
   Xm.tilde = Xtilde[na.ind, ]
   Xc.tilde = Xp.tilde[cen.ind, ]
   Xo.tilde = Xp.tilde[-cen.ind, ]
   Yp.tilde = Ytilde[-na.ind]
   Yo.tilde = Yp.tilde[-cen.ind]
   Utilde = Yp.tilde[cen.ind]

 # observed log-likelihood:
   TLam = TOme = matrix(0, ncol=n, nrow=n)
   TOme[1: cumsum.ni[1], 1: cumsum.ni[1]] = kronecker(Sigma, DEC(phi,ga, Data$Time[Data$Subject == 1]))
   Zi = matrix(TrZ[1: cumsum.ni[1], 1: q], ncol=q)
   TLam[1: cumsum.ni[1], 1: cumsum.ni[1]] = as.matrix(TrZ[1: cumsum.ni[1], 1:q]) %*% DD %*% t(TrZ[1: cumsum.ni[1], 1:q]) + TOme[1: cumsum.ni[1], 1: cumsum.ni[1]]
   for(i in 2: N){
    Zi = matrix(TrZ[(cumsum.ni[i-1]+1): cumsum.ni[i], ((i-1)*q+1): (i*q)], ncol=q)
    TOme[(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]] = kronecker(Sigma, DEC(phi,ga, Data$Time[Data$Subject == i]))
    TLam[(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]] = Zi %*% DD %*% t(Zi) + TOme[(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]]    
   }
   TLam.pp = TLam[-na.ind, -na.ind]
   TLam.mp = TLam[na.ind, -na.ind]
   TLam.mm = TLam[na.ind, na.ind]
   TLam.oo = TLam.pp[-cen.ind, -cen.ind]
   TLam.co = TLam.pp[cen.ind, -cen.ind]
   TLam.cc = TLam.pp[cen.ind, cen.ind]
   if(num.cen == 1){
     TLam.co = t(TLam.co); TLam.cc = as.matrix(TLam.cc)
   }
   if(num.na == 1) TLam.mp = t(TLam.mp) 
   TLam.oo.inv = solve(TLam.oo)
   yo.cent.tilde = Yo.tilde - Xo.tilde %*% beta
   mu.co = Xc.tilde %*% beta + TLam.co %*% TLam.oo.inv %*% yo.cent.tilde
   Sig.cc.o = TLam.cc - TLam.co %*% TLam.oo.inv %*% t(TLam.co)

   ind.lik = rep(0, N)
   if(ci[1]==1){
    ind.lik[cen.subj[1]] = log(pnorm(Utilde[1:cumsum.nc[1]], mean=mu.co[1:cumsum.nc[1]], sd=sqrt(Sig.cc.o[1:cumsum.nc[1], 1:cumsum.nc[1]])))
   } else ind.lik[cen.subj[1]] = log(pmvnorm(lower=rep(-Inf, ci[1]), upper=Utilde[1:cumsum.nc[1]], mean=c(mu.co[1:cumsum.nc[1]]), sigma=Sig.cc.o[1:cumsum.nc[1], 1:cumsum.nc[1]]))[1]
   if(Nc != 1){
   for(i in 2: Nc){
   if(ci[i]==1){ ind.lik[cen.subj[i]] = log(pnorm(Utilde[(cumsum.nc[i-1]+1): cumsum.nc[i]], mean=mu.co[(cumsum.nc[i-1]+1): cumsum.nc[i]], sd=sqrt(Sig.cc.o[(cumsum.nc[i-1]+1): cumsum.nc[i], (cumsum.nc[i-1]+1): cumsum.nc[i]])))
   } else ind.lik[cen.subj[i]] = log(pmvnorm(lower=rep(-Inf, ci[i]), upper=Utilde[(cumsum.nc[i-1]+1): cumsum.nc[i]], mean=c(mu.co[(cumsum.nc[i-1]+1): cumsum.nc[i]]), sigma=Sig.cc.o[(cumsum.nc[i-1]+1): cumsum.nc[i], (cumsum.nc[i-1]+1): cumsum.nc[i]]))[1]
   }}
   log.det.Sig.inv = log(det(TLam.oo.inv[1:cumsum.oi[1], 1:cumsum.oi[1]]))
   ind.lik[1] = ind.lik[1] - 0.5*log(2*pi)*oi[1] + 0.5*log.det.Sig.inv - 0.5*t(yo.cent.tilde[1:cumsum.oi[1]]) %*% TLam.oo.inv[1:cumsum.oi[1], 1:cumsum.oi[1]] %*% yo.cent.tilde[1:cumsum.oi[1]] 
   for(i in 2: N){
      log.det.Sig.inv = log(det(as.matrix(TLam.oo.inv[(cumsum.oi[i-1]+1): cumsum.oi[i], (cumsum.oi[i-1]+1): cumsum.oi[i]])))
      ind.lik[i] = ind.lik[i] - 0.5*log(2*pi)*oi[i] + 0.5*log.det.Sig.inv - 0.5*t(yo.cent.tilde[(cumsum.oi[i-1]+1): cumsum.oi[i]]) %*% TLam.oo.inv[(cumsum.oi[i-1]+1): cumsum.oi[i], (cumsum.oi[i-1]+1): cumsum.oi[i]] %*% yo.cent.tilde[(cumsum.oi[i-1]+1): cumsum.oi[i]]
   } 
   return(ind.lik)
}

# M-H algorithm for generating phi and ga
log.pseudo.pos.phiga = function(phi, ga, DD, Sigma, yp.cent, TrZ)
{
   y.na = vecData$Resp
   n = length(y.na)
   na.ind = which(is.na(y.na))
   Data.p = vecData[-na.ind, ]
   si = numeric(N)
   for(i in 1: N) si[i] = length(unique(Data$Time[Data$Subject == i]))
   ni = si * r
   cumsum.ni = cumsum(ni)
   pii = numeric(N)
   for(i in 1: N) pii[i] = length(unique(Data.p$Time[Data.p$Subject == i]))
   cumsum.pi = cumsum(pii)
   
   TLam = TOme = matrix(0, ncol=n, nrow=n)
   TOme[1: cumsum.ni[1], 1: cumsum.ni[1]] = kronecker(Sigma, DEC(phi, ga, Data$Time[Data$Subject == 1]))
   Zi = matrix(TrZ[1: cumsum.ni[1], 1: q], ncol=q)
   TLam[1: cumsum.ni[1], 1: cumsum.ni[1]] = as.matrix(TrZ[1: cumsum.ni[1], 1:q]) %*% DD %*% t(TrZ[1: cumsum.ni[1], 1:q]) + TOme[1: cumsum.ni[1], 1: cumsum.ni[1]]
   for(i in 2: N){
    Zi = matrix(TrZ[(cumsum.ni[i-1]+1): cumsum.ni[i], ((i-1)*q+1): (i*q)], ncol=q)
    TOme[(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]] = kronecker(Sigma, DEC(phi,ga, Data$Time[Data$Subject == i]))
    TLam[(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]] = Zi %*% DD %*% t(Zi) + TOme[(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]]    
   }
   TLam.pp = TLam[-na.ind, -na.ind]
   TLam.pp.inv = solve(TLam.pp)
   log.det.Lam.inv = log(det(TLam.pp.inv[1: cumsum.pi[1], 1: cumsum.pi[1]]))
   for(i in 2: N) log.det.Lam.inv = log.det.Lam.inv + log(det(TLam.pp.inv[(cumsum.pi[i-1]+1): cumsum.pi[i], (cumsum.pi[i-1]+1): cumsum.pi[i]]))  
   log.pos = 0.5*log.det.Lam.inv - 0.5*t(yp.cent) %*% TLam.pp.inv %*% yp.cent
   return(log.pos)
}

log.npos.phiga = function(phi, ga, DD, Sigma, yp.cent)
{
   y.na = vecData$Resp
   n = length(y.na)
   na.ind = which(is.na(y.na))
   Data.p = vecData[-na.ind, ]
   si = numeric(N)
   for(i in 1: N) si[i] = length(unique(Data$Time[Data$Subject == i]))
   ni = si * r
   cumsum.ni = cumsum(ni)
   pii = numeric(N)
   for(i in 1: N) pii[i] = length(unique(Data.p$Time[Data.p$Subject == i]))
   cumsum.pi = cumsum(pii)
   
   TOme = matrix(0, ncol=n, nrow=n)
   TOme[1: cumsum.ni[1], 1: cumsum.ni[1]] = kronecker(Sigma, DEC(phi, ga, Data$Time[Data$Subject == 1]))
   for(i in 2: N) TOme[(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]] = kronecker(Sigma, DEC(phi,ga, Data$Time[Data$Subject == i]))
   TOme.pp = TOme[-na.ind, -na.ind]
   TOme.pp.inv = solve(TOme.pp)
   log.det.Ome.inv = log(det(TOme.pp.inv[1: cumsum.pi[1], 1: cumsum.pi[1]]))
   for(i in 2: N) log.det.Ome.inv = log.det.Ome.inv + log(det(TOme.pp.inv[(cumsum.pi[i-1]+1): cumsum.pi[i], (cumsum.pi[i-1]+1): cumsum.pi[i]]))  
   log.pos = 0.5*log.det.Ome.inv - 0.5*t(yp.cent) %*% TOme.pp.inv %*% yp.cent
   return(log.pos)
}

MH.phi.nonli = function(phi, DD, Sigma, yp.cent, pg.mle, TrZ, method)
{
   phi.old = phi
   phi.star.old = log(phi.old) - log(1-phi.old)
   if(method=='pseudo') log.pos.old = log.pseudo.pos.phiga(phi.old, ga=1, DD, Sigma, yp.cent, TrZ)
   if(method=='raw') log.pos.old = log.npos.phiga(phi.old, ga=1, DD, Sigma, yp.cent)
   sd.phi.star = sqrt(2*(1/(pg.mle[1]*(1-pg.mle[1])))^2)
   phi.star.new = rnorm(1, phi.star.old, sd.phi.star)
   phi.new = exp(phi.star.new)/(1 + exp(phi.star.new))
   if(method=='pseudo') log.pos.new = log.pseudo.pos.phiga(phi.new, ga=1, DD, Sigma, yp.cent, TrZ)
   if(method=='raw') log.pos.new = log.npos.phiga(phi.new, ga=1, DD, Sigma, yp.cent) 
   log.accept = log.pos.new - log.pos.old
   if(log(runif(1)) > log.accept) phi.new = phi.old
   return(phi.new)
}

MH.DEC.nonli = function(phi, ga, DD, Sigma, yp.cent, pg.mle, TrZ, method)
{
   phi.ga.old = c(phi, ga)
   phi.ga.star.old = c(log(phi.ga.old[1])-log(1-phi.ga.old[1]), log(phi.ga.old[2]))
   if(method=='pseudo') log.pos.old = log.pseudo.pos.phiga(phi.ga.old[1], phi.ga.old[2], DD, Sigma, yp.cent, TrZ)
   if(method=='raw') log.pos.old = log.npos.phiga(phi.ga.old[1], phi.ga.old[2], DD, Sigma, yp.cent)
   cov.pg.star = matrix(c((1/(pg.mle[1]*(1-pg.mle[1])))^2, 0, 0, 1/pg.mle[2]^2), 2,2)
   repeat{
   phi.ga.star.new = as.vector(rmvnorm(1, phi.ga.star.old, cov.pg.star))
   if(phi.ga.star.new[2]>-1 & phi.ga.star.new[2]<0.5) break
   }
   phi.ga.new = c(exp(phi.ga.star.new[1])/(1 + exp(phi.ga.star.new[1])), exp(phi.ga.star.new[2]))
   if(method=='pseudo') log.pos.new = log.pseudo.pos.phiga(phi.ga.new[1], phi.ga.new[2], DD, Sigma, yp.cent, TrZ)
   if(method=='raw') log.pos.new = log.npos.phiga(phi.ga.new[1], phi.ga.new[2], DD, Sigma, yp.cent)
   log.accept = log.pos.new - log.pos.old
   if(log(runif(1)) > log.accept) phi.ga.new = phi.ga.old
   return(phi.ga.new)
}
