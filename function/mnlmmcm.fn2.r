library(nlme)
library(mvtnorm)
vech.posi=function(dim) cbind(rep(1:dim, 1:dim), unlist(mapply(':', 1, 1:dim)))
DEC = function(phi, ga, Ti) phi^((abs(outer(Ti, Ti, '-')))^ga)

mu.fn = function(eta, x)
{
   mu1 = log10(exp(eta[1]-eta[2]*x) + exp(eta[3]-eta[4]*x))
   mu2 = eta[5] + eta[6]*x
   return(list(mu1 = mu1, mu2 = mu2, mu = c(mu1, mu2)))
}

dmu = function(eta, x)
{
   p = length(eta)
   si = length(x)
   dmu1 = dmu2 = matrix(0, ncol=p, nrow=si)
   a1 = 1/((exp(eta[1]-eta[2]*x) + exp(eta[3]-eta[4]*x))*log(10))
   dmu1[,1] = a1 * exp(eta[1]-eta[2]*x)
   dmu1[,2] = - a1 * exp(eta[1]-eta[2]*x)*x
   dmu1[,3] = a1 * exp(eta[3]-eta[4]*x)
   dmu1[,4] = - a1 * exp(eta[3]-eta[4]*x)*x
   dmu2[,5] = 1
   dmu2[,6] = x
   return(list(dmu1 = dmu1, dmu2 = dmu2))
}

Q.phiga = function(zeta, Sigma, E.hat, cumsum.ni)
{
   pp = length(zeta)
   if(pp == 1){ phi = zeta; ga=1}
   else{ phi=zeta[1]; ga=zeta[2]}
   N = length(unique(Data$Subject))
   r = nrow(Sigma)
   sum1 = log(det(DEC(phi,ga, Data$Time[Data$Subject == 1])))
   sum2 = sum(solve(kronecker(Sigma, DEC(phi,ga, Data$Time[Data$Subject == 1])))*E.hat[1: cumsum.ni[1], 1: cumsum.ni[1]])
   for(i in 2: N){
     sum1 = sum1 + log(det(DEC(phi,ga, Data$Time[Data$Subject == i])))
     sum2 = sum2 + sum(solve(kronecker(Sigma, DEC(phi,ga, Data$Time[Data$Subject == i])))*E.hat[(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]])
   }
  opt.Qfn = r * sum1 + sum2
  return(opt.Qfn)
}

# ECM algorithm:
MNLMMcm.ECM = function(init.par, tol=1e-6, typeCi=c('UNC','AR1','DEC'), max.iter = 100000, per=1)
{
   begin = proc.time()[1]
 # initial values of parameters:
   N = length(unique(Data$Subject))
   beta = init.par$beta
   DD = init.par$DD
   Sigma = init.par$Sigma
   phi = init.par$phi
   ga = init.par$ga
   p = length(beta); q = ncol(DD); r = ncol(Sigma)
   vech.D = vech.posi(q)
   vech.Sig = vech.posi(r)
   b = rmvnorm(N, rep(0, q), DD)
   iter = 0
   theta.old = c(beta, DD[vech.D], Sigma[vech.Sig], phi, ga)

 # base settings:
   y = vecData$Resp
   n = length(y)
   na.ind = which(is.na(y))
   y[na.ind] = 999

   Data.o = vecData[which(vecData$Censor==0),]
   Data.c = vecData[which(vecData$Censor==1),]
   Data.m = vecData[na.ind,]
   Data.p = vecData[-na.ind, ]
   cen.ind = which(Data.p$Censor == 1)

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
   Nc = length(cen.subj)
   ci = numeric(Nc)
   for(i in 1: Nc) ci[i] = length(Data.c$Time[Data.c$Subject == cen.subj[i]])
   cumsum.nc = cumsum(ci)
   num.cen = length(cen.ind)
   num.na = length(na.ind)

   A = matrix(0, 6, p)
   B = diag(6)
   TrZ = matrix(0, ncol=N*q, nrow=n)
   MU = Xtilde = NULL
   for(i in 1: N){
     A[1, 1:2] = A[2, 3:4] = A[5, 9:10] = A[6, 11:12] = c(1, Data$arm[which(Data$Subject==i)][1])
     A[3, 5:6] = A[4, 7:8] = c(1, Data$lcd4[which(Data$Subject==i)][1])
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
   np = length(Yp.tilde)
   Yo.tilde = Yp.tilde[-cen.ind]
   no = length(Yo.tilde)
   Utilde = Yp.tilde[cen.ind]
   TP = diag(n)[-na.ind, ]
   TM = diag(n)[na.ind, ]
   TO = diag(np)[-cen.ind, ]
   TC = diag(np)[cen.ind, ]
   if(num.cen == 1) TC = t(TC)

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
#   if(num.cen == 1){
#    log.cdf = log(pnorm(Utilde, mean=mu.co, sd=sqrt(Sig.cc.o)))
#   } else log.cdf = log(pmvnorm(lower=-Inf, upper=Utilde, mean=c(mu.co), sigma=Sig.cc.o))[1]
   if(ci[1]==1){ log.cdf = log(pnorm(Utilde[1:cumsum.nc[1]], mean=mu.co[1:cumsum.nc[1]], sd=sqrt(Sig.cc.o[1:cumsum.nc[1], 1:cumsum.nc[1]])))
   } else log.cdf = log(pmvnorm(lower=rep(-Inf, ci[1]), upper=Utilde[1:cumsum.nc[1]], mean=c(mu.co[1:cumsum.nc[1]]), sigma=Sig.cc.o[1:cumsum.nc[1], 1:cumsum.nc[1]]))[1]
   if(Nc != 1){
   for(i in 2: Nc){
   if(ci[i]==1){ log.cdf = log.cdf + log(pnorm(Utilde[(cumsum.nc[i-1]+1): cumsum.nc[i]], mean=mu.co[(cumsum.nc[i-1]+1): cumsum.nc[i]], sd=sqrt(Sig.cc.o[(cumsum.nc[i-1]+1): cumsum.nc[i], (cumsum.nc[i-1]+1): cumsum.nc[i]])))
   } else log.cdf = log.cdf + log(pmvnorm(lower=rep(-Inf, ci[i]), upper=Utilde[(cumsum.nc[i-1]+1): cumsum.nc[i]], mean=c(mu.co[(cumsum.nc[i-1]+1): cumsum.nc[i]]), sigma=Sig.cc.o[(cumsum.nc[i-1]+1): cumsum.nc[i], (cumsum.nc[i-1]+1): cumsum.nc[i]]))[1]
   }}
   log.det.Sig.inv = log(det(TLam.oo.inv[1:cumsum.oi[1], 1:cumsum.oi[1]]))
   for(i in 2: N) log.det.Sig.inv = log.det.Sig.inv + log(det(as.matrix(TLam.oo.inv[(cumsum.oi[i-1]+1): cumsum.oi[i], (cumsum.oi[i-1]+1): cumsum.oi[i]])))
   old.loglik = log.cdf -0.5*log(2*pi)*no + 0.5*log.det.Sig.inv - 0.5*t(yo.cent.tilde) %*% TLam.oo.inv %*% yo.cent.tilde
   iter.lnL = old.loglik
   cat(rep('=', 25), 'MNLMM with ', typeCi, ' errors; ', 'censoring = ', num.cen/n*100, '%; missing = ', num.na/n*100, '%', rep('=', 25), sep = '', '\n')
   cat('iter = ', iter, ',\t obs.logli = ', old.loglik, sep = '', '\n')
   repeat
   {
     iter = iter+1
 # E-Step:
     EX = matrix(0, nrow=num.cen, ncol=1)
     EXX = matrix(0, nrow=num.cen, ncol=num.cen)
     Yci.hat = DTN.moment(mu=mu.co[1:cumsum.nc[1]], sigma=Sig.cc.o[1:cumsum.nc[1],1:cumsum.nc[1]], a.ast=rep(-Inf, ci[1]), b.ast=Utilde[1:cumsum.nc[1]])
     EX[1:cumsum.nc[1]] = Yci.hat$EX
     EXX[1:cumsum.nc[1], 1:cumsum.nc[1]] = Yci.hat$EXX
     if(Nc != 1){
     for(i in 2: Nc){
       Yci.hat = DTN.moment(mu=c(mu.co[(cumsum.nc[i-1]+1):cumsum.nc[i]]), sigma=Sig.cc.o[(cumsum.nc[i-1]+1):cumsum.nc[i],(cumsum.nc[i-1]+1):cumsum.nc[i]], a.ast=rep(-Inf, ci[i]), b.ast=Utilde[(cumsum.nc[i-1]+1):cumsum.nc[i]])
       EX[(cumsum.nc[i-1]+1):cumsum.nc[i]] = Yci.hat$EX
       EXX[(cumsum.nc[i-1]+1):cumsum.nc[i], (cumsum.nc[i-1]+1):cumsum.nc[i]] = Yci.hat$EXX
     }}
     TLam.pp.inv = solve(TLam.pp)
     Yhat.cent = t(TO)%*%Yo.tilde + t(TC)%*%EX - Xp.tilde%*%beta
     ym.hat = Xm.tilde %*% beta + TLam.mp %*% TLam.pp.inv %*% Yhat.cent
     y.hat = (t(TP) + t(TM)%*%TLam.mp%*%TLam.pp.inv) %*% (t(TO)%*%Yo.tilde + t(TC)%*%EX) + t(TM)%*%TM%*%(diag(n)-TLam%*%t(TP)%*%TLam.pp.inv%*%TP) %*% Xtilde %*% beta
     Sig.mm.p = TLam.mm - TLam.mp %*% TLam.pp.inv %*% t(TLam.mp)
     ym2.hat = ym.hat %*% t(ym.hat) + Sig.mm.p
     ycm.hat = EX %*% t(beta) %*% t(Xm.tilde) + (EX%*%t(Yo.tilde)%*%TO + EXX %*% TC - EX%*%t(beta)%*%t(Xp.tilde)) %*% TLam.pp.inv %*% t(TLam.mp)
     y2.hat = t(TP)%*%t(TO)%*% Yo.tilde %*% (t(Yo.tilde)%*%TO%*%TP + t(EX)%*%TC%*%TP + t(ym.hat)%*%TM) + t(TP)%*%t(TC) %*% (EX%*%t(Yo.tilde)%*%TO%*%TP + EXX%*%TC%*%TP + ycm.hat %*% TM) + t(TM)%*%(ym.hat%*%t(Yo.tilde)%*%TO%*%TP + t(ycm.hat)%*%TC%*%TP + ym2.hat%*%TM)

     TD = kronecker(diag(N), DD)   
     bb = TD %*% t(TrZ) %*% t(TP) %*% TLam.pp.inv %*% Yhat.cent
     TSig.b = TD - TD %*% t(TrZ) %*% t(TP) %*% TLam.pp.inv %*% TP %*% TrZ %*% TD
     b2 = TD %*% t(TrZ) %*% t(TP) %*% TLam.pp.inv %*% (Yhat.cent %*% t(Yhat.cent) - t(TC) %*% EX %*% t(EX) %*% TC + t(TC) %*% EXX %*% TC) %*% TLam.pp.inv %*% TP %*% TrZ %*% TD + TSig.b
     yb = (t(TP) %*% (t(TO)%*%Yo.tilde%*%t(Yhat.cent) + t(TC)%*%EX%*%t(t(TO)%*%Yo.tilde-Xp.tilde%*%beta) + t(TC)%*%EXX%*%TC) + t(TM) %*% (ym.hat%*%t(t(TO)%*%Yo.tilde-Xp.tilde%*%beta) + t(ycm.hat)%*%TC)) %*% TLam.pp.inv %*% TP %*% TrZ %*% TD    
     Xtil.beta = Xtilde %*% beta
     y_zbXbeta = (y.hat-TrZ%*%bb) %*% t(Xtil.beta)
     E.hat = y2.hat - yb%*%t(TrZ) - TrZ%*%t(yb) - t(y_zbXbeta) - y_zbXbeta + Xtil.beta%*%t(Xtil.beta) + TrZ%*%b2%*%t(TrZ)

 # CM-Steps:
     TOme.inv = solve(TOme)
     beta = solve(t(Xtilde) %*% TOme.inv %*% Xtilde) %*% (t(Xtilde) %*% TOme.inv %*% (y.hat - TrZ %*% bb))    # ECM
#     TLam.inv = solve(TLam)
#     beta = solve(t(Xtilde) %*% TLam.inv %*% Xtilde) %*% (t(Xtilde) %*% TLam.inv %*% y.hat)                  # AECM
     sum.b2 = 0
     for(i in 1: N) sum.b2 = sum.b2 + b2[((i-1)*q+1):(i*q), ((i-1)*q+1):(i*q)]
     DD = as.matrix(sum.b2 / N)
     if(r == 1){
     Ce=sum(solve(DEC(phi, ga, unique(Data$Time[Data$Subject == 1]))) * E.hat[1: cumsum.ni[1], 1: cumsum.ni[1]])
     for(i in 2: N) Ce = Ce + sum(solve(DEC(phi, ga, unique(Data$Time[Data$Subject == i]))) * E.hat[(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]])
     Sigma[1,1] = Ce / n
     } else{
     for(j in 1: r)for(l in 1:r)
     {
       Ce=0
       Ce=sum(solve(DEC(phi, ga, unique(Data$Time[Data$Subject == 1]))) * E.hat[1: cumsum.ni[1], 1: cumsum.ni[1]][((j-1)*si[1]+1): (j*si[1]), ((l-1)*si[1]+1): (l*si[1])])
       for(i in 2: N){ Ce = Ce + sum(solve(DEC(phi, ga, unique(Data$Time[Data$Subject == i]))) * E.hat[(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]][((j-1)*si[i]+1): (j*si[i]), ((l-1)*si[i]+1): (l*si[i])])
       }
       Sigma[j, l] = Ce /(n / r)
     }}
     if(typeCi == 'UNC'){ phi=1e-6; ga = 1}
     if(typeCi == 'AR1'){
        ga = 1
        phi = optim(par = phi, fn = Q.phiga, method = "L-BFGS-B", lower = 1e-6, upper = 1-1e-6, Sigma=Sigma, E.hat=E.hat, cumsum.ni=cumsum.ni)$par
      }
     if(typeCi == 'DEC'){
        par.DEC = optim(par = c(phi, ga), fn = Q.phiga, method = "L-BFGS-B", lower = c(1e-6, 1e-6), upper = c(1-1e-6, 2), Sigma=Sigma, E.hat=E.hat, cumsum.ni=cumsum.ni)$par
        phi = par.DEC[1]; ga = par.DEC[2]
     }
# evaluate new log-likelihood
     b = matrix(bb, ncol=q, byrow=T)
     MU = Xtilde = NULL
     for(i in 1: N){
       A[1, 1:2] = A[2, 3:4] = A[5, 9:10] = A[6, 11:12] = c(1, Data$arm[which(Data$Subject==i)][1])
       A[3, 5:6] = A[4, 7:8] = c(1, Data$lcd4[which(Data$Subject==i)][1])
       eta.i = A %*% beta + B %*% b[i, ]
       mu.i = mu.fn(eta.i, Data$Time[Data$Subject == i])
       MU = c(MU, mu.i$mu1, mu.i$mu2)
       dmu.i = dmu(eta.i, Data$Time[Data$Subject == i])
       Xtilde = rbind(Xtilde, dmu.i$dmu1 %*% A, dmu.i$dmu2 %*% A)
       if(i == 1) idx = 1: cumsum.ni[i]
       else idx = (cumsum.ni[i-1]+1): cumsum.ni[i]
       TrZ[idx, ((i-1)*q+1): (i*q)] = rbind(dmu.i$dmu1 %*% B, dmu.i$dmu2 %*% B)
     }
     Ytilde = y - MU + Xtilde %*% beta + TrZ %*% bb
     Xp.tilde = Xtilde[-na.ind, ]
     Xm.tilde = Xtilde[na.ind, ]
     Xc.tilde = Xp.tilde[cen.ind, ]
     Xo.tilde = Xp.tilde[-cen.ind, ]
     Yp.tilde = Ytilde[-na.ind]
     Yo.tilde = Yp.tilde[-cen.ind]
     Utilde = Yp.tilde[cen.ind]
     TLam = TOme = matrix(0, ncol=n, nrow=n)
     TOme[1: cumsum.ni[1], 1: cumsum.ni[1]] = kronecker(Sigma, DEC(phi,ga, Data$Time[Data$Subject == 1]))
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
#     if(num.cen == 1){
#       log.cdf = log(pnorm(Utilde, mean=mu.co, sd=sqrt(Sig.cc.o)))
#     } else log.cdf = log(pmvnorm(lower=-Inf, upper=Utilde, mean=c(mu.co), sigma=Sig.cc.o))[1]
     if(ci[1]==1){ log.cdf = log(pnorm(Utilde[1:cumsum.nc[1]], mean=mu.co[1:cumsum.nc[1]], sd=sqrt(Sig.cc.o[1:cumsum.nc[1], 1:cumsum.nc[1]])))
     } else log.cdf = log(pmvnorm(lower=rep(-Inf, ci[1]), upper=Utilde[1:cumsum.nc[1]], mean=c(mu.co[1:cumsum.nc[1]]), sigma=Sig.cc.o[1:cumsum.nc[1], 1:cumsum.nc[1]]))[1]
     if(Nc != 1){
     for(i in 2: Nc){
     if(ci[i]==1){ log.cdf = log.cdf + log(pnorm(Utilde[(cumsum.nc[i-1]+1): cumsum.nc[i]], mean=mu.co[(cumsum.nc[i-1]+1): cumsum.nc[i]], sd=sqrt(Sig.cc.o[(cumsum.nc[i-1]+1): cumsum.nc[i], (cumsum.nc[i-1]+1): cumsum.nc[i]])))
     } else log.cdf = log.cdf + log(pmvnorm(lower=rep(-Inf, ci[i]), upper=Utilde[(cumsum.nc[i-1]+1): cumsum.nc[i]], mean=c(mu.co[(cumsum.nc[i-1]+1): cumsum.nc[i]]), sigma=Sig.cc.o[(cumsum.nc[i-1]+1): cumsum.nc[i], (cumsum.nc[i-1]+1): cumsum.nc[i]]))[1]
     }}
     log.det.Sig.inv = log(det(TLam.oo.inv[1:cumsum.oi[1], 1:cumsum.oi[1]]))
     for(i in 2: N) log.det.Sig.inv = log.det.Sig.inv + log(det(as.matrix(TLam.oo.inv[(cumsum.oi[i-1]+1): cumsum.oi[i], (cumsum.oi[i-1]+1): cumsum.oi[i]])))
     new.loglik = log.cdf -0.5*log(2*pi)*no + 0.5*log.det.Sig.inv - 0.5*t(yo.cent.tilde) %*% TLam.oo.inv %*% yo.cent.tilde
     iter.lnL = c(iter.lnL, new.loglik)
     theta.new = c(beta, DD[vech.D], Sigma[vech.Sig], phi, ga)
     diff.lnL = (new.loglik - old.loglik)
     diff = sqrt(t(theta.new-theta.old) %*% (theta.new-theta.old))
     if(iter%%per == 0) cat('iter = ', iter, ',\t obs.loglik = ', new.loglik, ',\t diff.lnL=', diff.lnL, ',\t diff=', diff, ',\t beta=', beta, ',\t phi=', phi, ',\t ga=', ga, sep = ' ', '\n')
     if((diff < tol || diff.lnL < tol) || iter>=max.iter) break
     old.loglik = new.loglik; theta.old = theta.new
   }
   cat('iter = ', iter, ',\t obs.loglik = ', new.loglik, sep = '', '\n')
   end = proc.time()[1]
   run.sec = end - begin
   cat('It took', run.sec, 'seconds.\n')
   cat(paste(rep('-', 50), sep = '', collapse = ''), '\n')
   cat('beta =', beta, '\n')
   cat('DD =', DD[vech.D], '\n')
   cat('Sigma =', Sigma[vech.Sig], '\n')
   cat('phi =', phi, ',\t ga=', ga, '\n')
   cat(paste(rep('=', 50), sep = '', collapse = ''), '\n')
   yimp.hat = y.hat + MU - Xtilde %*% beta - TrZ %*% bb
   e = vecData$Resp - yimp.hat
   yfit = MU
   yfit[na.ind] = yimp.hat[na.ind]
   yfit[which(vecData$Censor==1)] = yimp.hat[which(vecData$Censor==1)]
   EST = c(beta, DD[vech.D], Sigma[vech.Sig], phi, ga)
   para = list(beta = beta, DD=DD, Sigma = Sigma, phi = phi, ga = ga, b=b)
   SD = IM.MNLMMcm(para, mc.size=500, typeCi=typeCi[1])
   if(typeCi == 'UNC') m = length(EST) - 2
   if(typeCi == 'AR1') m = length(EST) - 1
   if(typeCi == 'DEC') m = length(EST)
   aic = 2 * m - 2 * new.loglik
   bic = m * log(N) - 2 * new.loglik
   data.inf = list(cen.subj = cen.subj, na.subj=na.subj, cen.rate=num.cen/n*100, na.rate=num.na/n*100)
   model.inf = list(loglik = new.loglik, iter.lnL = iter.lnL, aic = aic, bic = bic, m=m)
   return(list(run.sec = run.sec, iter = iter, model.inf = model.inf, para = para, est = EST, SD = SD, yfit=yfit, yimp=yimp.hat, ym.hat=yimp.hat[na.ind], yc.hat=yimp.hat[-na.ind][cen.ind], error = e, data.inf=data.inf))
}

# Moment for truncated multivariate normal
DTN.moment = function(mu, sigma, a.ast, b.ast)
{
  require(mvtnorm)
  GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
  p = length(mu)
  a.ast = ifelse(a.ast==-Inf,rep(-1e12, p), a.ast)
  b.ast = ifelse(b.ast==Inf,rep(1e12, p), b.ast)
  if(p == 1){
    Sd = sqrt(sigma)
    lower.std=(a.ast - mu)/Sd
    upper.std=(b.ast - mu)/Sd
    EY = mu+Sd*(dnorm(lower.std)-dnorm(upper.std))/ (pnorm(upper.std)-pnorm(lower.std))
    variance=Sd^2*(1+(lower.std*(dnorm(lower.std))-upper.std*dnorm(upper.std))/(pnorm(upper.std)-pnorm(lower.std))-((dnorm(lower.std)-dnorm(upper.std))/(pnorm(upper.std)-pnorm(lower.std)))^2)
    EYY = variance + EY^2
  }
  else{
  Lambda = diag(sqrt(diag(sigma)))
  Lambda.inv=diag(1/sqrt(diag(sigma)))
  #Lambda.nh=diag(1/sqrt(diag(sigma)))
  #R=Lambda.nh %*%sigma%*% Lambda.nh
  R=Lambda.inv %*%sigma%*% Lambda.inv
  a= c(Lambda.inv %*% (a.ast-mu))
  b= c(Lambda.inv %*% (b.ast-mu))
  al0 = pmvnorm(lower = a, upper = b, sigma = R, algorithm = GB)[1]
### pdf & cdf
  f1a = dnorm(a)
  f1b = dnorm(b)
  f2 = matrix(NA, p, p)
  G1a = G1b = rep(NA, p)
  G2 = matrix(NA, p, p)
  for(r in 1:p)
  {
    temp = R[-r,r]
    S1 = R[-r,-r] - temp %*% t(R[r,-r])
    mua = temp * a[r]; low = a[-r]-mua; upp = b[-r]-mua
    G1a[r] = pmvnorm(lower = low, upper = upp, sigma = S1, algorithm = GB)[1]
    mub = temp * b[r]; low = a[-r]-mub; upp = b[-r]-mub
    G1b[r] = pmvnorm(lower = low, upper = upp, sigma = S1, algorithm = GB)[1]
  }
  qa = f1a*G1a; qb = f1b*G1b
  EX = c(R %*% (qa-qb)) / al0
  EY = c(mu + Lambda %*% EX)

  H = matrix(0,p,p)
  for(r in 1:(p-1))
  {
    for(s in (r+1):p)
    {
      rs = c(r,s)
      pdf.aa = dmvnorm(c(a[r],a[s]),sigma=round(R[rs,rs],5), log =F)
      pdf.ab = dmvnorm(c(a[r],b[s]),sigma=round(R[rs,rs],5), log =F)
      pdf.ba = dmvnorm(c(b[r],a[s]),sigma=round(R[rs,rs],5), log =F)
      pdf.bb = dmvnorm(c(b[r],b[s]),sigma=round(R[rs,rs],5), log =F)
      if(p==2){cdf.aa=cdf.ab=cdf.ba=cdf.bb=1}
      if(p>2)
      {
        tmp = R[-rs,rs]%*%solve(R[rs,rs])
        mu.aa = c(tmp%*%c(a[r],a[s]))
        mu.ab = c(tmp%*%c(a[r],b[s]))
        mu.ba = c(tmp%*%c(b[r],a[s]))
        mu.bb = c(tmp%*%c(b[r],b[s]))
        R21 = R[-rs,-rs] - R[-rs,rs]%*%solve(R[rs,rs]) %*% R[rs,-rs]
        cdf.aa = pmvnorm(lower = a[-rs], upper = b[-rs], mean=mu.aa, sigma = R21, algorithm = GB)[1]
        cdf.ab = pmvnorm(lower = a[-rs], upper = b[-rs], mean=mu.ab, sigma = R21, algorithm = GB)[1]
        cdf.ba = pmvnorm(lower = a[-rs], upper = b[-rs], mean=mu.ba, sigma = R21, algorithm = GB)[1]
        cdf.bb = pmvnorm(lower = a[-rs], upper = b[-rs], mean=mu.bb, sigma = R21, algorithm = GB)[1]
      }
      H[r,s] = H[s,r] = pdf.aa*cdf.aa - pdf.ab*cdf.ab - pdf.ba*cdf.ba + pdf.bb*cdf.bb
    }
  }
  D = matrix(0,p,p)
  diag(D) = a * qa - b * qb - diag(R%*%H)
  EXX = R + R %*% (H + D) %*% R / al0
  EYY = mu%*%t(mu) + Lambda%*%EX%*%t(mu) + mu%*%t(EX)%*%Lambda + Lambda%*%EXX%*%Lambda
  }
  CovY=EYY-(EY)%*%t(EY)
  return(list(EX=EY,EXX=EYY, CovX=CovY))
}

IM.MNLMMcm = function(para.est, mc.size=500, typeCi=c('UNC','AR1','DEC'))
{
# setting
   beta = para.est$beta
   DD = para.est$DD
   Sigma = para.est$Sigma 
   phi = para.est$phi
   ga = para.est$ga
   b = para.est$b
   N = length(unique(Data$Subject))
   y = vecData$Resp
   na.ind = which(is.na(y))
   n = length(y)
   p = length(beta); q = ncol(DD); r = ncol(Sigma)

   Data.o = vecData[which(vecData$Censor==0),]
   Data.c = vecData[which(vecData$Censor==1),]
   Data.m = vecData[na.ind,]
   Data.p = vecData[-na.ind, ]
   cen.ind = which(Data.p$Censor == 1)

   si = oi = pii = numeric(N)
   for(i in 1: N) si[i] = length(unique(Data$Time[Data$Subject == i]))
   ni = si * r
   cumsum.ni = cumsum(ni)
   for(i in 1: N) oi[i] = length(Data.o$Time[Data.o$Subject == i])
   cumsum.oi = cumsum(oi)
   for(i in 1: N) pii[i] = length(Data.p$Time[Data.p$Subject == i])
   cumsum.pi = cumsum(pii)
   IDX = as.list(N)
   IDX[[1]] = 1: cumsum.ni[1] 
   for(i in 2: N) IDX[[i]] = (cumsum.ni[i-1]+1): cumsum.ni[i]

   cen.subj = unique(Data.c$Subject)
   cen.subj = as.numeric(levels(cen.subj))[cen.subj]
   na.subj = unique(Data.m$Subject)
   na.subj = as.numeric(levels(na.subj))[na.subj]
   cm.subj = sort(unique(c(cen.subj, na.subj)))
   Nc = length(cen.subj)
   ci = numeric(Nc)
   for(i in 1: Nc) ci[i] = length(Data.c$Time[Data.c$Subject == cen.subj[i]])
   cumsum.nc = cumsum(ci)
   Nm = length(na.subj)
   mi = numeric(Nm)
   for(i in 1: Nm) mi[i] = length(Data.m$Time[Data.m$Subject == na.subj[i]])
   cumsum.na = cumsum(mi)
   num.cen = length(cen.ind)
   num.na = length(na.ind)

   cen.idx = na.idx = as.list(N)
   for(i in 1: N) cen.idx[[i]] = na.idx[[i]] = NA
   cen.idx[[cen.subj[1]]] = 1:cumsum.nc[1]
   for(i in 2:Nc) cen.idx[[cen.subj[i]]] = (cumsum.nc[i-1]+1): cumsum.nc[i]
   na.idx[[na.subj[[1]]]] = 1: cumsum.na[1]
   for(i in 2:Nm) na.idx[[na.subj[i]]] = (cumsum.na[i-1]+1): cumsum.na[i] 
   CM = as.list(N)
   for(i in 1: N) CM[[i]] = vecData$Censor[vecData$Subject==i]

   A = matrix(0, 6, p)
   B = diag(6)
   TrZ = matrix(0, ncol=N*q, nrow=n)
   MU = Xtilde = NULL
   for(i in 1: N){
     A[1, 1:2] = A[2, 3:4] = A[5, 9:10] = A[6, 11:12] = c(1, Data$arm[which(Data$Subject==i)][1])
     A[3, 5:6] = A[4, 7:8] = c(1, Data$lcd4[which(Data$Subject==i)][1])
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
   np = length(Yp.tilde)
   Yo.tilde = Yp.tilde[-cen.ind]
   no = length(Yo.tilde)
   Utilde = Yp.tilde[cen.ind]
   TP = diag(n)[-na.ind, ]
   TM = diag(n)[na.ind, ]
   TO = diag(np)[-cen.ind, ]
   TC = diag(np)[cen.ind, ]
   if(num.cen == 1) TC = t(TC)

   TLam = TOme = matrix(0, ncol=n, nrow=n)
   TOme[1: cumsum.ni[1], 1: cumsum.ni[1]] = kronecker(Sigma, DEC(phi,ga, unique(Data$Time[Data$Subject == 1])))
   Zi = matrix(TrZ[1: cumsum.ni[1], 1:q], ncol=q)
   TLam[1: cumsum.ni[1], 1: cumsum.ni[1]] = as.matrix(TrZ[1: cumsum.ni[1], 1:q]) %*% DD %*% t(TrZ[1: cumsum.ni[1], 1:q]) + TOme[1: cumsum.ni[1], 1: cumsum.ni[1]]
   for(i in 2: N){
    Zi = matrix(TrZ[(cumsum.ni[i-1]+1): cumsum.ni[i], ((i-1)*q+1): (i*q)], ncol=q)
    TOme[(cumsum.ni[i-1]+1): cumsum.ni[i], (cumsum.ni[i-1]+1): cumsum.ni[i]] = kronecker(Sigma, DEC(phi,ga, unique(Data$Time[Data$Subject == i])))
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
   TLam.inv = solve(TLam)
   TLam.pp.inv = solve(TLam.pp)
   TLam.oo.inv = solve(TLam.oo)

 #  Information matrix
   p = length(beta)
   g1 = q*(q+1)/2
   g2 = r*(r+1)/2
   g = g1+g2+2
   dot.L = as.list(numeric(g))
   dot.L = as.list(matrix(0, N, g))
   vechD = vech.posi(q) 
   vechS = vech.posi(r) 
 
   for(l in 1: g1)
   {
     dot.DD = matrix(0, q, q)
     dot.DD[matrix(vechD[l, ], 1)] = dot.DD[matrix(rev(vechD[l, ]), 1)] = 1
     dot.L[[(l-1)*N+1]] = TrZ[IDX[[1]], 1:q] %*% dot.DD %*% t(TrZ[IDX[[1]], 1:q])
     for(i in 2: N) dot.L[[(l-1)*N+i]] = TrZ[IDX[[i]], ((i-1)*q+1): (i*q)] %*% dot.DD %*% t(TrZ[IDX[[i]], ((i-1)*q+1): (i*q)])
   }
   for(l in 1: g2)
   {
     dot.Sig = matrix(0, r, r)
     dot.Sig[matrix(vechS[l, ], 1)] = dot.Sig[matrix(rev(vechS[l, ]), 1)] = 1
     for(i in 1: N) dot.L[[g1*N+(l-1)*N+i]] = kronecker(dot.Sig, DEC(phi, ga, Data$Time[Data$Subject == i]))
   }
   for(i in 1: N) dot.L[[(g1+g2)*N+i]] = kronecker(Sigma, DEC.dot.phi(phi, ga, Data$Time[Data$Subject == i]))
   for(i in 1: N) dot.L[[(g-1)*N+i]] = kronecker(Sigma, DEC.dot.ga(phi, ga, Data$Time[Data$Subject == i]))

# Generate Monte Carlo sample for yc and ym
   yo.cent.tilde = Yo.tilde - Xo.tilde %*% beta
   mu.co = Xc.tilde %*% beta + TLam.co %*% TLam.oo.inv %*% yo.cent.tilde
   Sig.cc.o = TLam.cc - TLam.co %*% TLam.oo.inv %*% t(TLam.co)

   yp.cent.tilde = t(TO)%*%Yo.tilde + t(TC)%*%mu.co - Xp.tilde%*%beta
   mu.mp = Xm.tilde %*% beta + TLam.mp %*% TLam.pp.inv %*% yp.cent.tilde
   Sig.mm.p = TLam.mm - TLam.mp %*% TLam.pp.inv %*% t(TLam.mp)

   y.samp = matrix(rep(y, mc.size), nrow=mc.size, ncol=n, byrow=T)
   for(m in 1: mc.size){
   for(i in cm.subj){
    # Yc #
   if(!is.na(cen.idx[i])){ 
      yc.hat = c(rtmvnorm(1, mean=c(mu.co[cen.idx[[i]]]), sigma=Sig.cc.o[cen.idx[[i]], cen.idx[[i]]], lower=rep(-Inf, length(cen.idx[[i]])), upper=Utilde[cen.idx[[i]]], algorithm="gibbs", burn.in.samples=10))
      y.samp[m, IDX[[i]][which(CM[[i]]==1)]] = yc.hat
   }
   # Ym #
   if(!is.na(na.idx[[i]])[1]){ 
      lna = length(na.idx[[i]])
      if(lna == 1) ym.hat = c(rnorm(1, mean=mu.mp[na.idx[[i]]], sd=sqrt(Sig.mm.p[na.idx[[i]], na.idx[[i]]])))
      if(lna > 1) ym.hat = c(rmvnorm(1, mu.mp[na.idx[[i]]], round(Sig.mm.p[na.idx[[i]], na.idx[[i]]],6)))
      y.samp[m, IDX[[i]][which(is.na(CM[[i]]))]] = ym.hat
   }}}

# Compute the score vector and Hessian matrix 
   H = SS = matrix(0, ncol=(p+g), nrow=(p+g))
   S = numeric((p+g))
   for(m in 1: mc.size){
# H.beta
   H[1:p, 1:p] = H[1:p, 1:p] - t(Xtilde) %*% TLam.inv %*% Xtilde
   ys.cent = y.samp[m, ] - Xtilde%*%beta
# S.beta
   sb = t(Xtilde)%*%TLam.inv%*%ys.cent  
   S[1:p] = S[1:p] + sb    
   SS[1:p, 1:p] = SS[1:p, 1:p] + sb %*% t(sb)

# H.xi
   Linv.dotL = as.list(numeric(N*g))
   for(s in 1: g)for(l in 1: s){
     for(i in 1: N){
       Linv.dotL[[(s-1)*N+i]] = TLam.inv[IDX[[i]], IDX[[i]]] %*% dot.L[[(s-1)*N+i]]
       H[1:p, (p+s)] = H[1:p, (p+s)] + t(Xtilde[IDX[[i]], ])%*%Linv.dotL[[(s-1)*N+i]]%*%TLam.inv[IDX[[i]],IDX[[i]]]%*%ys.cent[IDX[[i]]]
       H[p+s,p+l] = H[p+s,p+l] - 0.5*sum(diag(Linv.dotL[[(s-1)*N+i]] %*% Linv.dotL[[(l-1)*N+i]]))-0.5*sum(diag((ys.cent[IDX[[i]]])%*%t(ys.cent[IDX[[i]]])%*%(Linv.dotL[[(s-1)*N+i]] %*% Linv.dotL[[(l-1)*N+i]]%*%TLam.inv[IDX[[i]],IDX[[i]]]+Linv.dotL[[(l-1)*N+i]]%*%Linv.dotL[[(s-1)*N+i]]%*%TLam.inv[IDX[[i]],IDX[[i]]])))
     }}
   for(l in (p+1):(p+g-1))for(s in (l+1):(p+g)) H[l, s] = H[s, l]
   H[(p+1):(p+g), 1:p] = t(H[1:p, (p+1):(p+g)])
   for(s in 1: g){
     for(i in 1: N){
       Linv.dotL[[(s-1)*N+i]] = TLam.inv[IDX[[i]], IDX[[i]]] %*% dot.L[[(s-1)*N+i]]
       S[p+s] = S[p+s] -0.5*sum(diag(Linv.dotL[[(s-1)*N+i]])) + 0.5 * sum(diag(ys.cent[IDX[[i]]]%*%t(ys.cent[IDX[[i]]])%*%Linv.dotL[[(s-1)*N+i]]%*%TLam.inv[IDX[[i]], IDX[[i]]]))
   }}
   SS[(p+1):(p+g), (p+1):(p+g)] = SS[(p+1):(p+g), (p+1):(p+g)] + S[-c(1:p)] %*% t(S[-c(1:p)])
  }
  H.mean = H/mc.size
  SS.mean = SS/mc.size
  S.mean = S/mc.size

  I.theta = -H.mean - SS.mean  + S.mean%*%t(S.mean)
  if(typeCi == 'UNC'){
   V.theta = solve(I.theta[1:(p+g-2),1:(p+g-2)])
   sd.theta = c(sqrt(diag(V.theta)), 0, 0)
  }
  if(typeCi == 'AR1'){
   V.theta = solve(I.theta[1:(p+g-1),1:(p+g-1)])
   sd.theta = c(sqrt(diag(V.theta)), 0)
  }
  if(typeCi == 'DEC'){
   V.theta = solve(I.theta)
   sd.theta = sqrt(diag(V.theta))
  }
  EST = c(beta, DD[vechD], Sigma[vechS], phi, ga)
  out = rbind(EST, sd.theta)
  SD = list(out=out, I.theta=I.theta, V.theta=V.theta)
  FI=list(Itheta = round(I.theta, 4), Ibeta = I.theta[1:p, 1:p], Iomega = I.theta[-(1:p), -(1:p)])
  return(list(out=out, FI=FI, SD = SD))
}
