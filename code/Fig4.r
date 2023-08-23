# Reproduct Figure 4
load(paste(PATH, 'data/A5055FitResult.RData',sep=""))

cenna = intersect(IBCM0$data.inf$cen.subj, IBCM0$data.inf$na.subj)
yfit = IBCM1$pos.inf$yfit.out
ym.hat = IBCM1$pos.inf$Ym.out
yfitM = IBM2$pos.summ$yfit.out
ym.hatM = IBM2$pos.summ$Ym.out

Nna.i = numeric(N)
for(i in 1: N) Nna.i[i] = sum(is.na(Data$cd4.cd8[Data$Subj==i]))
cumsum.na = c(0, cumsum(Nna.i))

sel = c(2, 27, 14, 34, 17, 44, 22)
PATH = paste(getwd(),"/Data_and_Code_METR23/",sep="")
postscript(paste(PATH, 'results/fig4.eps', sep=''), width=15, height=12)
par(mfcol=c(4,4))
for(i in 1: 7){
# log10(RNA)
par(mar=c(0,4,4,0.5))
  a1 = c(Data$lgcopy[Data$Subj==sel[i]], yfit[c(1,3,4), which(vecData$Subj==sel[i] & vecData$Var==1)])
  plot(Data$Day[Data$Subj==sel[i]], Data$lgcopy[Data$Subj==sel[i]], type='n', xlab='Day', ylab=expression(log[10](RNA)), main=paste('Patient ', sel[i], sep=''), ylim=c(min(a1), max(a1)), las=1, xaxt='n', cex.main=2, font.main=2)
  axis(1, Data$Day[Data$Subj==sel[i]], labels=F)
  polygon(c(Data$Day[Data$Subj==sel[i]], rev(Data$Day[Data$Subj==sel[i]])), c(yfit[3, which(vecData$Subj==sel[i] & vecData$Var==1)], rev(yfit[4, which(vecData$Subj==sel[i] & vecData$Var==1)])), col="lightcyan", border="cyan")
  polygon(c(Data$Day[Data$Subj==sel[i]], rev(Data$Day[Data$Subj==sel[i]])), c(yfitM[3, which(vecData$Subj==sel[i] & vecData$Var==1)], rev(yfitM[4, which(vecData$Subj==sel[i] & vecData$Var==1)])), col="pink", border="pink3")

  points(Data$Day[Data$Subj==sel[i]], yfitM[1, which(vecData$Subj==sel[i] & vecData$Var==1)], pch=9, col='magenta', cex=1.5)
  lines(Data$Day[Data$Subj==sel[i]], yfitM[1, which(vecData$Subj==sel[i] & vecData$Var==1)], lty=1, col='magenta', lwd=1)
  points(Data$Day[Data$Subj==sel[i]], yfit[1, which(vecData$Subj==sel[i] & vecData$Var==1)], pch=7, col='blue', cex=1.5)
  lines(Data$Day[Data$Subj==sel[i]], yfit[1, which(vecData$Subj==sel[i] & vecData$Var==1)], lty=1, col='blue', lwd=1)

  points(Data$Day[Data$Subj==sel[i]], Data$lgcopy[Data$Subj==sel[i]], pch=20, col=1, cex=1.5)
  abline(h=log10(50), lty=3, col='gray25', lwd=1.5)

# CD4/CD8
par(mar=c(4,4,0,0.5))
  a2 = na.omit(c(Data$cd4.cd8[Data$Subj==sel[i]], yfit[c(1,3,4), which(vecData$Subj==sel[i] & vecData$Var==2)], ym.hat[3:4, (cumsum.na[sel[i]]+1):cumsum.na[sel[i]+1]], ym.hatM[3:4, (cumsum.na[sel[i]]+1):cumsum.na[sel[i]+1]]))
  plot(Data$Day[Data$Subj==sel[i]], Data$cd4.cd8[Data$Subj==sel[i]], type='n', xlab='Day', ylab='CD4/CD8', ylim=c(min(a2), max(a2)), las=1, xaxt='n')
  axis(1, Data$Day[Data$Subj==sel[i]])
  polygon(c(Data$Day[Data$Subj==sel[i]], rev(Data$Day[Data$Subj==sel[i]])), c(yfit[3, which(vecData$Subj==sel[i] & vecData$Var==2)], rev(yfit[4, which(vecData$Subj==sel[i] & vecData$Var==2)])), col="lightcyan", border=5)
  polygon(c(Data$Day[Data$Subj==sel[i]], rev(Data$Day[Data$Subj==sel[i]])), c(yfitM[3, which(vecData$Subj==sel[i] & vecData$Var==2)], rev(yfitM[4, which(vecData$Subj==sel[i] & vecData$Var==2)])), col="pink", border="pink3")

  points(Data$Day[Data$Subj==sel[i]], yfitM[1, which(vecData$Subj==sel[i] & vecData$Var==2)], pch=9, col='magenta', cex=1.5)
  lines(Data$Day[Data$Subj==sel[i]], yfitM[1, which(vecData$Subj==sel[i] & vecData$Var==2)], lty=1, col='magenta', lwd=1)
  points(Data$Day[Data$Subj==sel[i]], yfit[1, which(vecData$Subj==sel[i] & vecData$Var==2)], pch=7, col='blue', cex=1.5)
  lines(Data$Day[Data$Subj==sel[i]], yfit[1, which(vecData$Subj==sel[i] & vecData$Var==2)], lty=1, col='blue', lwd=1)

  loc.na = which(is.na(Data$cd4.cd8[Data$Subj==sel[i]]))
  points(Data$Day[Data$Subj==sel[i]][loc.na], ym.hatM[1, (cumsum.na[sel[i]]+1):cumsum.na[sel[i]+1]], pch=18, col=2, cex=2)
  points(rep(Data$Day[Data$Subj==sel[i]][loc.na], each=2), ym.hatM[3:4, (cumsum.na[sel[i]]+1):cumsum.na[sel[i]+1]], pch='--', col=2, cex=2)
  segments(Data$Day[Data$Subj==sel[i]][loc.na], ym.hatM[3, (cumsum.na[sel[i]]+1):cumsum.na[sel[i]+1]], Data$Day[Data$Subj==sel[i]][loc.na], ym.hatM[4, (cumsum.na[sel[i]]+1):cumsum.na[sel[i]+1]], lty=2, col=2, lwd=1.5)
  points(Data$Day[Data$Subj==sel[i]][loc.na], ym.hat[1, (cumsum.na[sel[i]]+1):cumsum.na[sel[i]+1]], pch=15, col="blue3", cex=1.5)
  points(rep(Data$Day[Data$Subj==sel[i]][loc.na], each=2), ym.hat[3:4, (cumsum.na[sel[i]]+1):cumsum.na[sel[i]+1]], pch='--', col="blue3", cex=2)
  segments(Data$Day[Data$Subj==sel[i]][loc.na], ym.hat[3, (cumsum.na[sel[i]]+1):cumsum.na[sel[i]+1]], Data$Day[Data$Subj==sel[i]][loc.na], ym.hat[4, (cumsum.na[sel[i]]+1):cumsum.na[sel[i]+1]], lty=2, col="blue3", lwd=2)

  points(Data$Day[Data$Subj==sel[i]], Data$cd4.cd8[Data$Subj==sel[i]], pch=20, col=1, cex=1.5)
}
plot(0:1, 0:1, type='n', axes=F, bty='n', xlab='', ylab='')
legend(0, 0.8, c('Observations', 'Fitted values under MNLMM-CM', 'Fitted values under MNLMM',
       '95% PI for missing under MNLMM-CM','95% PI for missing under MNLMM'),
        col=c('black', 'blue', 'magenta','blue3','red'), lty=c(NA, 1, 1, 2, 2), pch=c(16, 7, 9, 15, 18), bty='n', cex=0.93)
legend(0, 0.3, c(' 95% PB under MNLMM-CM', ' 95% PB under MNLMM'), fill=c('lightcyan', 'pink'), border=c('cyan', 'pink3'), bty='n', cex=0.93)
legend(0, 0.1, 'Detection limit = log10(50)', lty=3, pch=NA, col='gray25', bty='n', cex=0.95)
dev.off()
