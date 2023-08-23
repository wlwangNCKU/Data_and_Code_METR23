library("vioplot")
PATH1 = paste(PATH, 'results/N=25/', sep='')
PATH2 = paste(PATH, 'results/N=50/', sep='')
PATH3 = paste(PATH, 'results/N=75/', sep='')
PATH4 = paste(PATH, 'results/N=100/', sep='')

#N=25
mseyfit25<-read.table(paste(PATH1,'mseyfit.txt',sep=""),na.strings="NA",sep="")[,-1]
n25<-rep(25, 600)
mseyfit25<-cbind(mseyfit25, n25)
colnames(mseyfit25)<-c("case","mleCM","bayCM","mleM","bayM","n")

#N=50
mseyfit50<-read.table(paste(PATH2,'mseyfit.txt',sep=""),na.strings="NA",sep="")[,-1]
n50<-rep(50, 600)
mseyfit50<-cbind(mseyfit50, n50)
colnames(mseyfit50)<-c("case","mleCM","bayCM","mleM","bayM","n")

#N=75
mseyfit75<-read.table(paste(PATH3,'mseyfit.txt',sep=""),na.strings="NA",sep="")[,-1]
n75<-rep(75, 600)
mseyfit75<-cbind(mseyfit75, n75)
colnames(mseyfit75)<-c("case","mleCM","bayCM","mleM","bayM","n")

#N=100
mseyfit100<-read.table(paste(PATH4,'mseyfit.txt',sep=""),na.strings="NA",sep="")[,-1]
n100<-rep(100, 600)
mseyfit100<-cbind(mseyfit100, n100)
colnames(mseyfit100)<-c("case","mleCM","bayCM","mleM","bayM","n")

mse.yfit<-rbind(mseyfit25, mseyfit50, mseyfit75, mseyfit100)
MCAR10.C10<-mse.yfit[mse.yfit$case=='MCAR10.C10',]
MAR.C10<-mse.yfit[mse.yfit$case=='MAR.C10',]
MCAR10.C30<-mse.yfit[mse.yfit$case=='MCAR10.C30',]
MAR.C30<-mse.yfit[mse.yfit$case=='MAR.C30',]
MCAR10.C50<-mse.yfit[mse.yfit$case=='MCAR10.C50',]
MAR.C50<-mse.yfit[mse.yfit$case=='MAR.C50',]

################################################################################
# --------------------------- Figure 2 for Simulation --------------------------
################################################################################
postscript(paste(PATH, 'results/SIMfig2.eps', sep=''), height=4, width=12)
layout(matrix(c(1,1,1,8,8,8,2,3,4,9,10,11,5,6,7,12,13,14), 3, 6, byrow=T),
       widths=c(6,4.5,4.5,6,4.5,4.5), heights=c(0.75,0.75,4))
N=c(25,50,75,100)
par(mar=c(0, 5.2, 0.5, 0.5), cex.lab=2)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', ylab='', xlab='')
mtext('MAR', side=1, line=-2, cex=1.5, font=2)
par(mar=c(0, 5.2, 0, 0), cex.lab=2)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', ylab='', xlab='')
mtext('10% Censoring', 1, line=-2.3, cex=1.2, font=2)
par(mar=c(0, 0, 0, 0), cex.lab=2)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', ylab='', xlab='')
mtext('30% Censoring', 1, line=-2.3, cex=1.2, font=2)
par(mar=c(0, 0, 0, 0.5), cex.lab=2)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', ylab='', xlab='')
mtext('50% Censoring', 1, line=-2.3, cex=1.2, font=2)

par(mar=c(4, 5.2, 0, 0), cex.lab=2)
vioplot(bayM~n, data=MAR.C10, col = 0, plotCentre = "line", side = "left",ylim=c(0.6,3),yaxt="n",xlab="", ylab="",las=1,border="blue",cex.axis=1.2)
vioplot(bayCM~n, data=MAR.C10, col = "palevioletred", plotCentre = "line", side = "right", add = T,border="red")
axis(2, at=seq(1,3,0.5))
title(xlab = "Sample Size", ylab = expression(MSE~(hat(y))),cex.lab=1.2)
legend("topright", fill = c(0,"palevioletred"), legend = c("MNLMM","MNLMM-CM"), cex=1, bty='n',border=c("blue","red"))

par(mar=c(4, 0, 0, 0), cex.lab=2)
vioplot(bayM~n, data=MAR.C30, col = 0, plotCentre = "line", side = "left",ylim=c(0.6,3),yaxt="n",xlab="", ylab="",las=1,border="blue")
axis(2,labels=F, at=seq(1,3,0.5))
vioplot(bayCM~n, data=MAR.C30, col = "palevioletred", plotCentre = "line", side = "right", add = T,border="red")
title(xlab = "Sample Size", ylab = expression(MSE~(hat(y))),cex.lab=1.2)
legend("topright", fill = c(0,"palevioletred"), legend = c("MNLMM","MNLMM-CM"), cex=1, bty='n',border=c("blue","red"))

par(mar=c(4, 0, 0, 0.5), cex.lab=2)
vioplot(bayM~n, data=MAR.C50, col = 0, plotCentre = "line", side = "left",ylim=c(0.6,3),yaxt="n",xlab="", ylab="",las=1,border="blue")
axis(2,labels=F, at=seq(1,3,0.5))
vioplot(bayCM~n, data=MAR.C50, col = "palevioletred", plotCentre = "line", side = "right", add = T,border="red")
title(xlab = "Sample Size", ylab = expression(MSE~(hat(y))),cex.lab=1.2)
legend("topright", fill = c(0,"palevioletred"), legend = c("MNLMM","MNLMM-CM"), cex=1, bty='n',border=c("blue","red"))

par(mar=c(0, 5.2, 0.5, 0.5), cex.lab=2)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', ylab='', xlab='')
mtext('MCAR', side=1, line=-2, cex=1.5, font=2)
par(mar=c(0, 5.2, 0, 0), cex.lab=2)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', ylab='', xlab='')
mtext('10% Censoring', 1, line=-2.3, cex=1.2, font=2)
par(mar=c(0, 0, 0, 0), cex.lab=2)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', ylab='', xlab='')
mtext('30% Censoring', 1, line=-2.3, cex=1.2, font=2)
par(mar=c(0, 0, 0, 0.5), cex.lab=2)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', ylab='', xlab='')
mtext('50% Censoring', 1, line=-2.3, cex=1.2, font=2)

par(mar=c(4, 5.2, 0, 0), cex.lab=2)
vioplot(bayM~n, data=MCAR10.C10, col = 0, plotCentre = "line", side = "left",ylim=c(0.6,3),xlab="", ylab="",las=1,border="blue",cex.axis=1.2)
vioplot(bayCM~n, data=MCAR10.C10, col = "palevioletred", plotCentre = "line", side = "right", add = T,border="red")
title(xlab = "Sample Size", ylab = expression(MSE~~(hat(y))),cex.lab=1.2)
legend("topright", fill = c(0,"palevioletred"), legend = c("MNLMM","MNLMM-CM"), cex=1, bty='n',border=c("blue","red"))
 
par(mar=c(4, 0, 0, 0), cex.lab=2)
vioplot(bayM~n, data=MCAR10.C30, col = 0, plotCentre = "line", side = "left",ylim=c(0.6,3),yaxt="n",xlab="", ylab="",las=1,border="blue")
axis(2,labels=F, at=seq(1,3,0.5))
vioplot(bayCM~n, data=MCAR10.C30, col = "palevioletred", plotCentre = "line", side = "right", add = T,border="red")
title(xlab = "Sample Size", ylab = expression(MSE~~(hat(y))),cex.lab=1.2)
legend("topright", fill = c(0,"palevioletred"), legend = c("MNLMM","MNLMM-CM"), cex=1, bty='n',border=c("blue","red"))

par(mar=c(4, 0, 0, 0.5), cex.lab=2)
vioplot(bayM~n, data=MCAR10.C50, col = 0, plotCentre = "line", side = "left",ylim=c(0.6,3),yaxt="n",xlab="", ylab="",las=1,border="blue")
axis(2,labels=F, at=seq(1,3,0.5))
vioplot(bayCM~n, data=MCAR10.C50, col = "palevioletred", plotCentre = "line", side = "right", add = T,border="red")
title(xlab = "Sample Size", ylab = expression(MSE~~(hat(y))),cex.lab=1.2)
legend("topright", fill = c(0,"palevioletred"), legend = c("MNLMM","MNLMM-CM"), cex=1, bty='n',border=c("blue","red"))
dev.off()
