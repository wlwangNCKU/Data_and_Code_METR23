PATH1 = paste(PATH, 'results/N=25/', sep='')
PATH2 = paste(PATH, 'results/N=50/', sep='')
PATH3 = paste(PATH, 'results/N=75/', sep='')
PATH4 = paste(PATH, 'results/N=100/', sep='')

MSE = function(model) (model$Est-model$real.data)^2
n<-rep(c(25,50,75,100), each=7800)

#################################### bayCM #####################################
bayCM25<-read.table(paste(PATH1, 'table(bayCM).txt', sep=""), header=F, sep="")
bayCM50<-read.table(paste(PATH2, 'table(bayCM).txt', sep=""), header=F, sep="")
bayCM75<-read.table(paste(PATH3, 'table(bayCM).txt', sep=""), header=F, sep="")
bayCM100<-read.table(paste(PATH4, 'table(bayCM).txt', sep=""), header=F, sep="")

bayCM<-rbind(bayCM25, bayCM50, bayCM75, bayCM100)
bayCM<-cbind(bayCM, n)
colnames(bayCM)<-c("loop","para","case","real.data","Est","SE","LCI","UCI","Len","CP","n")

bayCM.b1<-bayCM[which(bayCM$para=="b1"),]; bayCM.b2<-bayCM[which(bayCM$para=="b2"),]
bayCM.b3<-bayCM[which(bayCM$para=="b3"),]; bayCM.b4<-bayCM[which(bayCM$para=="b4"),]
bayCM.b5<-bayCM[which(bayCM$para=="b5"),]; bayCM.d11<-bayCM[which(bayCM$para=="d11"),]
bayCM.d21<-bayCM[which(bayCM$para=="d21"),]; bayCM.d22<-bayCM[which(bayCM$para=="d22"),]
bayCM.sigma11<-bayCM[which(bayCM$para=="sigma11"),]; bayCM.sigma21<-bayCM[which(bayCM$para=="sigma21"),]
bayCM.sigma22<-bayCM[which(bayCM$para=="sigma22"),]

### para MSE
bayCM.b1mse<-MSE(bayCM.b1); bayCM.b2mse<-MSE(bayCM.b2); bayCM.b3mse<-MSE(bayCM.b3)
bayCM.b4mse<-MSE(bayCM.b4); bayCM.b5mse<-MSE(bayCM.b5); bayCM.d11mse<-MSE(bayCM.d11)
bayCM.d21mse<-MSE(bayCM.d21); bayCM.d22mse<-MSE(bayCM.d22); bayCM.sigma11mse<-MSE(bayCM.sigma11)
bayCM.sigma21mse<-MSE(bayCM.sigma21); bayCM.sigma22mse<-MSE(bayCM.sigma22)

bayCM.mse<-cbind(bayCM.b1$case, bayCM.b1mse, bayCM.b2mse, bayCM.b3mse, bayCM.b4mse, bayCM.b5mse,
                 bayCM.d11mse, bayCM.d21mse, bayCM.d22mse, 
                 bayCM.sigma11mse, bayCM.sigma21mse, bayCM.sigma22mse, bayCM.b1$n)
bayCM.mse<-data.frame(bayCM.mse)
colnames(bayCM.mse)<-c("case","b1","b2","b3","b4","b5","d11","d21","d22","sigma11","sigma21","sigma22","n")

#################################### bayM #####################################
bayM25<-read.table(paste(PATH1, 'table(bayM).txt', sep=""), header=F, sep="")
bayM50<-read.table(paste(PATH2, 'table(bayM).txt', sep=""), header=F, sep="")
bayM75<-read.table(paste(PATH3, 'table(bayM).txt', sep=""), header=F, sep="")
bayM100<-read.table(paste(PATH4, 'table(bayM).txt', sep=""), header=F, sep="")

bayM<-rbind(bayM25, bayM50, bayM75, bayM100)
bayM<-cbind(bayM, n)
colnames(bayM)<-c("loop","para","case","real.data","Est","SE","LCI","UCI","Len","CP","n")

bayM.b1<-bayM[which(bayM$para=="b1"),]; bayM.b2<-bayM[which(bayM$para=="b2"),]
bayM.b3<-bayM[which(bayM$para=="b3"),]; bayM.b4<-bayM[which(bayM$para=="b4"),]
bayM.b5<-bayM[which(bayM$para=="b5"),]; bayM.d11<-bayM[which(bayM$para=="d11"),]
bayM.d21<-bayM[which(bayM$para=="d21"),]; bayM.d22<-bayM[which(bayM$para=="d22"),]
bayM.sigma11<-bayM[which(bayM$para=="sigma11"),]; bayM.sigma21<-bayM[which(bayM$para=="sigma21"),]
bayM.sigma22<-bayM[which(bayM$para=="sigma22"),]

### para MSE
bayM.b1mse<-MSE(bayM.b1); bayM.b2mse<-MSE(bayM.b2); bayM.b3mse<-MSE(bayM.b3)
bayM.b4mse<-MSE(bayM.b4); bayM.b5mse<-MSE(bayM.b5); bayM.d11mse<-MSE(bayM.d11)
bayM.d21mse<-MSE(bayM.d21); bayM.d22mse<-MSE(bayM.d22); bayM.sigma11mse<-MSE(bayM.sigma11)
bayM.sigma21mse<-MSE(bayM.sigma21); bayM.sigma22mse<-MSE(bayM.sigma22)

bayM.mse<-cbind(bayM.b1$case, bayM.b1mse, bayM.b2mse, bayM.b3mse, bayM.b4mse, bayM.b5mse,
                bayM.d11mse, bayM.d21mse, bayM.d22mse, 
                bayM.sigma11mse, bayM.sigma21mse, bayM.sigma22mse, bayM.b1$n)
bayM.mse<-data.frame(bayM.mse)
colnames(bayM.mse)<-c("case","b1","b2","b3","b4","b5","d11","d21","d22","sigma11","sigma21","sigma22","n")

#######################################################################################################
n<-c(25,50,75,100)
model.case = function(model){data.frame(rbind(t(colMeans(model[which(model$n=="25"),c(2:12)],na.rm=T)),
                                              t(colMeans(model[which(model$n=="50"),c(2:12)],na.rm=T)),
                                              t(colMeans(model[which(model$n=="75"),c(2:12)],na.rm=T)),
                                              t(colMeans(model[which(model$n=="100"),c(2:12)],na.rm=T))))}

### case1=MAR.C10
bayCM.C1<-bayCM.mse[which(bayCM.mse$case=="1"),]; bayM.C1<-bayM.mse[which(bayM.mse$case=="1"),]
bayCM.case1<-cbind(model.case(bayCM.C1), n); bayM.case1<-cbind(model.case(bayM.C1), n)

### case2=MACR.C10
bayCM.C2<-bayCM.mse[which(bayCM.mse$case=="2"),]; bayM.C2<-bayM.mse[which(bayM.mse$case=="2"),]
bayCM.case2<-cbind(model.case(bayCM.C2), n); bayM.case2<-cbind(model.case(bayM.C2), n)

### case3=MAR.C30
bayCM.C3<-bayCM.mse[which(bayCM.mse$case=="3"),]; bayM.C3<-bayM.mse[which(bayM.mse$case=="3"),]
bayCM.case3<-cbind(model.case(bayCM.C3), n); bayM.case3<-cbind(model.case(bayM.C3), n)

### case4=MACR.C30
bayCM.C4<-bayCM.mse[which(bayCM.mse$case=="4"),]; bayM.C4<-bayM.mse[which(bayM.mse$case=="4"),]
bayCM.case4<-cbind(model.case(bayCM.C4), n); bayM.case4<-cbind(model.case(bayM.C4), n)

### case5=MAR.C50
bayCM.C5<-bayCM.mse[which(bayCM.mse$case=="5"),]; bayM.C5<-bayM.mse[which(bayM.mse$case=="5"),]
bayCM.case5<-cbind(model.case(bayCM.C5), n); bayM.case5<-cbind(model.case(bayM.C5), n)

### case6=MACR.C50
bayCM.C6<-bayCM.mse[which(bayCM.mse$case=="6"),]; bayM.C6<-bayM.mse[which(bayM.mse$case=="6"),]
bayCM.case6<-cbind(model.case(bayCM.C6), n); bayM.case6<-cbind(model.case(bayM.C6), n)

################################################################################
# ------------------------- Figure S.1 for Simulation --------------------------
################################################################################
postscript(paste(PATH, 'results/SIMfigS1.eps', sep=''), width=6.5, height=40)
N = c(25,50,75,100)
layout(matrix(c(1,2,2,2,3,3,3,1,4:86), 13, 7, byrow=T),
       widths=c(1.8,rep(2,6)), heights=c(2,2,rep(4.5,10),8.5))

par(mar=c(0, 1, 0, 0), cex.lab=1)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
legend("center",lty=c(1,2),pch=c(19,15), col=c("red", "blue"),
       legend = c("MNLMM-CM", "MNLMM"), bty='n', cex=0.7, lwd=1, pt.cex=1)

par(mar=c(0, 0, 0, 0), cex.lab=2)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', ylab='', xlab='')
mtext('MAR', side=1, line=-1.2, cex=1.1, font=2)
par(mar=c(0, 0, 0, 0), cex.lab=2)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', ylab='', xlab='')
mtext('MCAR', side=1, line=-1.2, cex=1.1, font=2)

par(mar=c(0, 0, 0, 0), cex.lab=0.9)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', ylab='', xlab='')
mtext('10% Censoring', side=1, line=-1.45, cex=0.6, font=2)
par(mar=c(0, 0, 0, 0), cex.lab=0.9)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', ylab='', xlab='')
mtext('30% Censoring', side=1, line=-1.45, cex=0.6, font=2)
par(mar=c(0, 0, 0, 0), cex.lab=0.9)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', ylab='', xlab='')
mtext('50% Censoring', side=1, line=-1.45, cex=0.6, font=2)
par(mar=c(0, 0, 0, 0), cex.lab=0.9)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', ylab='', xlab='')
mtext('10% Censoring', side=1, line=-1.45, cex=0.6, font=2)
par(mar=c(0, 0, 0, 0), cex.lab=0.9)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', ylab='', xlab='')
mtext('30% Censoring', side=1, line=-1.45, cex=0.6, font=2)
par(mar=c(0, 0, 0, 0), cex.lab=0.9)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', ylab='', xlab='')
mtext('50% Censoring', side=1, line=-1.45, cex=0.6, font=2)

### beta 1
par(mar=c(0, 1, 0, 0), cex.lab=1)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', ylab='', xlab='')
mtext(expression(bold(beta)[1]), side=1, line=-2.4,at=0.25, cex=1, font=2)

par(mar=c(0, 0, 0, 0), cex.lab=1)
plot(bayCM.case1$b1,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.055),xlim=c(1,5))
axis(2,labels=c(0,0.025,0.05),at=c(0,0.025,0.05),cex.axis=0.7,las=1)
lines(bayM.case1$b1,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case3$b1,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.129),xlim=c(1,5))
axis(2,labels=c(0,0.06,0.12),at=c(0,0.06,0.12),cex.axis=0.7,las=1)
lines(bayM.case3$b1,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case5$b1,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.85),xlim=c(1,5))
axis(2,labels=c(0,0.4,0.8),at=c(0,0.4,0.8),cex.axis=0.7,las=1)
lines(bayM.case5$b1,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case2$b1,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.055),xlim=c(1,5))
axis(2,labels=c(0,0.025,0.05),at=c(0,0.025,0.05),cex.axis=0.7,las=1)
lines(bayM.case2$b1,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case4$b1,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.275),xlim=c(1,5))
axis(2,labels=c(0,0.125,0.25),at=c(0,0.125,0.25),cex.axis=0.7,las=1)
lines(bayM.case4$b1,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case6$b1,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.85),xlim=c(1,5))
axis(2,labels=c(0,0.4,0.8),at=c(0,0.4,0.8),cex.axis=0.7,las=1)
lines(bayM.case6$b1,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

### beta 2
par(mar=c(0, 1, 0, 0), cex.lab=1)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', ylab='', xlab='')
mtext(expression(bold(beta)[2]), side=1, line=-2.4,at=0.25, cex=1, font=2)

par(mar=c(0, 0, 0, 0), cex.lab=1)
plot(bayCM.case1$b2,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,9),xlim=c(1,5))
axis(2,labels=c(0,4,8),at=c(0,4,8),cex.axis=0.7,las=1)
lines(bayM.case1$b2,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case3$b2,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,140),xlim=c(1,5))
axis(2,labels=c(0,66,132),at=c(0,66,132),cex.axis=0.7,las=1)
lines(bayM.case3$b2,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case5$b2,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,11800),xlim=c(1,5))
axis(2,labels=c(0,500,1000),at=c(0,5000,10000),cex.axis=0.7,las=1)
lines(bayM.case5$b2,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case2$b2,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,9),xlim=c(1,5))
axis(2,labels=c(0,4.1,8.2),at=c(0,4.1,8.2),cex.axis=0.7,las=1)
lines(bayM.case2$b2,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case4$b2,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,470),xlim=c(1,5))
axis(2,labels=c(0,210,420),at=c(0,210,420),cex.axis=0.7,las=1)
lines(bayM.case4$b2,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case6$b2,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,11800),xlim=c(1,5))
axis(2,labels=c(0,500,1000),at=c(0,5000,10000),cex.axis=0.7,las=1)
lines(bayM.case6$b2,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

### beta 3
par(mar=c(0, 1, 0, 0), cex.lab=1)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', ylab='', xlab='')
mtext(expression(bold(beta)[3]), side=1, line=-2.4,at=0.25, cex=1, font=2)

par(mar=c(0, 0, 0, 0), cex.lab=1)
plot(bayCM.case1$b3,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,7.3),xlim=c(1,5))
axis(2,labels=c(0,3.5,7),at=c(0,3.5,7),cex.axis=0.7,las=1)
lines(bayM.case1$b3,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case3$b3,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,60),xlim=c(1,5))
axis(2,labels=c(0,27,54),at=c(0,27,54),cex.axis=0.7,las=1)
lines(bayM.case3$b3,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case5$b3,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,4200),xlim=c(1,5))
axis(2,labels=c(0,2000,4000),at=c(0,2000,4000),cex.axis=0.7,las=1)
lines(bayM.case5$b3,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case2$b3,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,7.3),xlim=c(1,5))
axis(2,labels=c(0,3.5,7),at=c(0,3.5,7),cex.axis=0.7,las=1)
lines(bayM.case2$b3,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case4$b3,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,360),xlim=c(1,5))
axis(2,labels=c(0,158,316),at=c(0,158,316),cex.axis=0.7,las=1)
lines(bayM.case4$b3,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case6$b3,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,4200),xlim=c(1,5))
axis(2,labels=c(0,2000,4000),at=c(0,2000,4000),cex.axis=0.7,las=1)
lines(bayM.case6$b3,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

### beta 4
par(mar=c(0, 1, 0, 0), cex.lab=1)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', ylab='', xlab='')
mtext(expression(bold(beta)[4]), side=1, line=-2.4,at=0.25, cex=1, font=2)

par(mar=c(0, 0, 0, 0), cex.lab=1)
plot(bayCM.case1$b4,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.085),xlim=c(1,5))
axis(2,labels=c(0,0.027,0.054),at=c(0,0.038,0.076),cex.axis=0.7,las=1)
lines(bayM.case1$b4,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case3$b4,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.078),xlim=c(1,5))
axis(2,labels=c(0,0.035,0.07),at=c(0,0.035,0.07),cex.axis=0.7,las=1)
lines(bayM.case3$b4,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case5$b4,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.121),xlim=c(1,5))
axis(2,labels=c(0,0.05,0.1),at=c(0,0.05,0.1),cex.axis=0.7,las=1)
lines(bayM.case5$b4,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case2$b4,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.075),xlim=c(1,5))
axis(2,labels=c(0,0.034,0.068),at=c(0,0.034,0.068),cex.axis=0.7,las=1)
lines(bayM.case2$b4,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case4$b4,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.082),xlim=c(1,5))
axis(2,labels=c(0,0.037,0.074),at=c(0,0.037,0.074),cex.axis=0.7,las=1)
lines(bayM.case4$b4,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case6$b4,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.121),xlim=c(1,5))
axis(2,labels=c(0,0.05,0.1),at=c(0,0.05,0.1),cex.axis=0.7,las=1)
lines(bayM.case6$b4,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

### beta 5
par(mar=c(0, 1, 0, 0), cex.lab=1)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', ylab='', xlab='')
mtext(expression(bold(beta)[5]), side=1, line=-2.4,at=0.25, cex=1, font=2)

par(mar=c(0, 0, 0, 0), cex.lab=1)
plot(bayCM.case1$b5,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,7.3e-06),xlim=c(1,5))
axis(2,labels=c(0,expression(2%*%10^-6),expression(5%*%10^-6)),at=c(0,2e-6,5e-6),cex.axis=0.7,las=1)
lines(bayM.case1$b5,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case3$b5,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,6e-06),xlim=c(1,5))
axis(2,labels=c(0,expression(2%*%10^-6),expression(5%*%10^-6)),at=c(0,2e-6,5e-6),cex.axis=0.7,las=1)
lines(bayM.case3$b5,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case5$b5,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,2e-05),xlim=c(1,5))
axis(2,labels=c(0,expression(1%*%10^-5),expression(2%*%10^-5)),at=c(0,1e-5,2e-5),cex.axis=0.7,las=1)
lines(bayM.case5$b5,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case2$b5,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,5.8e-6),xlim=c(1,5))
axis(2,labels=c(0,expression(2.6%*%10^-6),expression(5.2%*%10^-6)),at=c(0,2.6e-6,5.2e-6),cex.axis=0.7,las=1)
lines(bayM.case2$b5,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case4$b5,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,1.1e-05),xlim=c(1,5))
axis(2,labels=c(0,expression(5%*%10^-6),expression(10^-5)),at=c(0,5e-6,1e-5),cex.axis=0.7,las=1)
lines(bayM.case4$b5,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case6$b5,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,2e-05),xlim=c(1,5))
axis(2,labels=c(0,expression(1%*%10^-5),expression(2%*%10^-5)),at=c(0,1e-5,2e-5),cex.axis=0.7,las=1)
lines(bayM.case6$b5,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

### d11
par(mar=c(0, 1, 0, 0), cex.lab=1)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', ylab='', xlab='')
mtext(expression(d[11]), side=1, line=-2.4,at=0.25, cex=1, font=2)

par(mar=c(0, 0, 0, 0), cex.lab=1)
plot(bayCM.case1$d11,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.125),xlim=c(1,5))
axis(2,labels=c(0,0.06,0.12),at=c(0,0.06,0.12),cex.axis=0.7,las=1)
lines(bayM.case1$d11,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case3$d11,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.27),xlim=c(1,5))
axis(2,labels=c(0,0.13,0.26),at=c(0,0.13,0.26),cex.axis=0.7,las=1)
lines(bayM.case3$d11,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case5$d11,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.55),xlim=c(1,5))
axis(2,labels=c(0,0.25,0.5),at=c(0,0.25,0.5),cex.axis=0.7,las=1)
lines(bayM.case5$d11,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case2$d11,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.125),xlim=c(1,5))
axis(2,labels=c(0,0.06,0.12),at=c(0,0.06,0.12),cex.axis=0.7,las=1)
lines(bayM.case2$d11,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case4$d11,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.27),xlim=c(1,5))
axis(2,labels=c(0,0.13,0.26),at=c(0,0.13,0.26),cex.axis=0.7,las=1)
lines(bayM.case4$d11,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case6$d11,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.55),xlim=c(1,5))
axis(2,labels=c(0,0.25,0.5),at=c(0,0.25,0.5),cex.axis=0.7,las=1)
lines(bayM.case6$d11,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

### d21
par(mar=c(0, 1, 0, 0), cex.lab=1)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', ylab='', xlab='')
mtext(expression(d[21]), side=1, line=-2.4,at=0.25, cex=1, font=2)

par(mar=c(0, 0, 0, 0), cex.lab=1)
plot(bayCM.case1$d21,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.065),xlim=c(1,5))
axis(2,labels=c(0,0.029,0.058),at=c(0,0.029,0.058),cex.axis=0.7,las=1)
lines(bayM.case1$d21,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case3$d21,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.065),xlim=c(1,5))
axis(2,labels=c(0,0.029,0.058),at=c(0,0.029,0.058),cex.axis=0.7,las=1)
lines(bayM.case3$d21,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case5$d21,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.065),xlim=c(1,5))
axis(2,labels=c(0,0.029,0.058),at=c(0,0.029,0.058),cex.axis=0.7,las=1)
lines(bayM.case5$d21,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case2$d21,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.065),xlim=c(1,5))
axis(2,labels=c(0,0.029,0.058),at=c(0,0.029,0.058),cex.axis=0.7,las=1)
lines(bayM.case2$d21,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case4$d21,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.065),xlim=c(1,5))
axis(2,labels=c(0,0.029,0.058),at=c(0,0.029,0.058),cex.axis=0.7,las=1)
lines(bayM.case4$d21,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case6$d21,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.065),xlim=c(1,5))
axis(2,labels=c(0,0.029,0.058),at=c(0,0.029,0.058),cex.axis=0.7,las=1)
lines(bayM.case6$d21,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

### d22
par(mar=c(0, 1, 0, 0), cex.lab=1)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', ylab='', xlab='')
mtext(expression(d[22]), side=1, line=-2.4,at=0.25, cex=1, font=2)

par(mar=c(0, 0, 0, 0), cex.lab=1)
plot(bayCM.case1$d22,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.14),xlim=c(1,5))
axis(2,labels=c(0,0.063,0.126),at=c(0,0.063,0.126),cex.axis=0.7,las=1)
lines(bayM.case1$d22,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case3$d22,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.14),xlim=c(1,5))
axis(2,labels=c(0,0.063,0.126),at=c(0,0.063,0.126),cex.axis=0.7,las=1)
lines(bayM.case3$d22,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case5$d22,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.14),xlim=c(1,5))
axis(2,labels=c(0,0.063,0.126),at=c(0,0.063,0.126),cex.axis=0.7,las=1)
lines(bayM.case5$d22,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case2$d22,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.14),xlim=c(1,5))
axis(2,labels=c(0,0.063,0.126),at=c(0,0.063,0.126),cex.axis=0.7,las=1)
lines(bayM.case2$d22,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case4$d22,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.14),xlim=c(1,5))
axis(2,labels=c(0,0.063,0.126),at=c(0,0.063,0.126),cex.axis=0.7,las=1)
lines(bayM.case4$d22,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case6$d22,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.14),xlim=c(1,5))
axis(2,labels=c(0,0.063,0.126),at=c(0,0.063,0.126),cex.axis=0.7,las=1)
lines(bayM.case6$d22,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

### sigma11
par(mar=c(0, 1, 0, 0), cex.lab=1)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', ylab='', xlab='')
mtext(expression(bold(sigma)[11]), side=1, line=-2.4,at=0.25, cex=1, font=2)

par(mar=c(0, 0, 0, 0), cex.lab=1)
plot(bayCM.case1$sigma11,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.38),xlim=c(1,5))
axis(2,labels=c(0,0.17,0.34),at=c(0,0.17,0.34),cex.axis=0.7,las=1)
lines(bayM.case1$sigma11,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case3$sigma11,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.6),xlim=c(1,5))
axis(2,labels=c(0,0.26,0.52),at=c(0,0.26,0.52),cex.axis=0.7,las=1)
lines(bayM.case3$sigma11,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case5$sigma11,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,1.1),xlim=c(1,5))
axis(2,labels=c(0,0.5,1),at=c(0,0.5,1),cex.axis=0.7,las=1)
lines(bayM.case5$sigma11,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case2$sigma11,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.38),xlim=c(1,5))
axis(2,labels=c(0,0.17,0.34),at=c(0,0.17,0.34),cex.axis=0.7,las=1)
lines(bayM.case2$sigma11,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case4$sigma11,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,2.1),xlim=c(1,5))
axis(2,labels=c(0,0.9,1.8),at=c(0,0.9,1.8),cex.axis=0.7,las=1)
lines(bayM.case4$sigma11,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case6$sigma11,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,1.1),xlim=c(1,5))
axis(2,labels=c(0,0.5,1),at=c(0,0.5,1),cex.axis=0.7,las=1)
lines(bayM.case6$sigma11,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

### sigma21
par(mar=c(0, 1, 0, 0), cex.lab=1)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', ylab='', xlab='')
mtext(expression(bold(sigma)[21]), side=1, line=-2.4,at=0.25, cex=1, font=2)

par(mar=c(0, 0, 0, 0), cex.lab=1)
plot(bayCM.case1$sigma21,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.08),xlim=c(1,5))
axis(2,labels=c(0,0.033,0.066),at=c(0,0.033,0.066),cex.axis=0.7,las=1)
lines(bayM.case1$sigma21,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case3$sigma21,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.028),xlim=c(1,5))
axis(2,labels=c(0,0.012,0.024),at=c(0,0.012,0.024),cex.axis=0.7,las=1)
lines(bayM.case3$sigma21,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case5$sigma21,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.1),xlim=c(1,5))
axis(2,labels=c(0,0.04,0.08),at=c(0,0.04,0.08),cex.axis=0.7,las=1)
lines(bayM.case5$sigma21,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case2$sigma21,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.085),xlim=c(1,5))
axis(2,labels=c(0,0.036,0.072),at=c(0,0.036,0.072),cex.axis=0.7,las=1)
lines(bayM.case2$sigma21,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case4$sigma21,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.047),xlim=c(1,5))
axis(2,labels=c(0,0.02,0.04),at=c(0,0.02,0.04),cex.axis=0.7,las=1)
lines(bayM.case4$sigma21,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case6$sigma21,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,0.1),xlim=c(1,5))
axis(2,labels=c(0,0.04,0.08),at=c(0,0.04,0.08),cex.axis=0.7,las=1)
lines(bayM.case6$sigma21,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

### sigma22
par(mar=c(4, 1, 0, 0), cex.lab=1)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', ylab='', xlab='')
mtext(expression(bold(sigma)[22]), side=1, line=-2.4,at=0.25, cex=1, font=2)

par(mar=c(4, 0, 0, 0), cex.lab=1)
plot(bayCM.case1$sigma22,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,1.22),xlim=c(1,5))
axis(1,labels=N,at=1:4,cex.axis=0.75)
axis(2,labels=c(0,0.56,1.12),at=c(0,0.56,1.12),cex.axis=0.7,las=1)
lines(bayM.case1$sigma22,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case3$sigma22,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,1.42),xlim=c(1,5))
axis(1,labels=N,at=1:4,cex.axis=0.75)
axis(2,labels=c(0,0.66,1.32),at=c(0,0.66,1.32),cex.axis=0.7,las=1)
lines(bayM.case3$sigma22,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case5$sigma22,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,1.42),xlim=c(1,5))
axis(1,labels=N,at=1:4,cex.axis=0.75)
axis(2,labels=c(0,0.66,1.32),at=c(0,0.66,1.32),cex.axis=0.7,las=1)
lines(bayM.case5$sigma22,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case2$sigma22,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,1.22),xlim=c(1,5))
axis(1,labels=N,at=1:4,cex.axis=0.75)
axis(2,labels=c(0,0.56,1.12),at=c(0,0.56,1.12),cex.axis=0.7,las=1)
lines(bayM.case2$sigma22,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case4$sigma22,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,1.24),xlim=c(1,5))
axis(1,labels=N,at=1:4,cex.axis=0.75)
axis(2,labels=c(0,0.57,1.14),at=c(0,0.57,1.14),cex.axis=0.7,las=1)
lines(bayM.case4$sigma22,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)

plot(bayCM.case6$sigma22,type="o",col="red",xaxt="n",yaxt="n",xlab="Sample Size",ylab="MSE",cex.axis=0.7,lty=1,pch=19,cex=1.2,lwd=1.2,las=2,ylim=c(0,1.24),xlim=c(1,5))
axis(1,labels=N,at=1:4,cex.axis=0.75)
axis(2,labels=c(0,0.57,1.14),at=c(0,0.57,1.14),cex.axis=0.7,las=1)
lines(bayM.case6$sigma22,type="o",col="blue",lty=2,pch=15,cex=1.2,lwd=1.2,las=2)
dev.off()
