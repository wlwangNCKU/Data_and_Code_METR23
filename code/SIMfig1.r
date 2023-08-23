PATH1 = paste(PATH, 'results/N=25/', sep='')
PATH2 = paste(PATH, 'results/N=50/', sep='')
PATH3 = paste(PATH, 'results/N=75/', sep='')
PATH4 = paste(PATH, 'results/N=100/', sep='')

vech.posi=function(dim) cbind(rep(1:dim, 1:dim), unlist(mapply(':', 1, 1:dim)))
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
phi = 1e-6; ga = 1  
par = list(beta=beta, DD=DD, Sigma=Sigma, phi=phi, ga=ga)
par.true = c(beta, DD[vech.posi(q)], Sigma[vech.posi(r)], phi, ga)

##################################### N=25 #####################################
bay.m25<-read.table(paste(PATH1, "table(bayM).txt",sep=""))
colnames(bay.m25)<-c("loop","para","case","real data","Est","SE","LCI","UCI","Len","CP")

bay.cm25<-read.table(paste(PATH1, "table(bayCM).txt",sep=""))
colnames(bay.cm25)<-c("loop","para","case","real data","Est","SE","LCI","UCI","Len","CP")

baycm.c1.25<-bay.cm25[bay.cm25$case=='1',]; baycm.c2.25<-bay.cm25[bay.cm25$case=='2',]
baycm.c3.25<-bay.cm25[bay.cm25$case=='3',]; baycm.c4.25<-bay.cm25[bay.cm25$case=='4',]
baycm.c5.25<-bay.cm25[bay.cm25$case=='5',]; baycm.c6.25<-bay.cm25[bay.cm25$case=='6',]

baym.c1.25<-bay.m25[bay.m25$case=='1',]; baym.c2.25<-bay.m25[bay.m25$case=='2',]
baym.c3.25<-bay.m25[bay.m25$case=='3',]; baym.c4.25<-bay.m25[bay.m25$case=='4',]
baym.c5.25<-bay.m25[bay.m25$case=='5',]; baym.c6.25<-bay.m25[bay.m25$case=='6',]

##################################### N=50 #####################################
bay.m50<-read.table(paste(PATH2, "table(bayM).txt",sep=""))
colnames(bay.m50)<-c("loop","para","case","real data","Est","SE","LCI","UCI","Len","CP")

bay.cm50<-read.table(paste(PATH2, "table(bayCM).txt",sep=""))
colnames(bay.cm50)<-c("loop","para","case","real data","Est","SE","LCI","UCI","Len","CP")

baycm.c1.50<-bay.cm50[bay.cm50$case=='1',]; baycm.c2.50<-bay.cm50[bay.cm50$case=='2',]
baycm.c3.50<-bay.cm50[bay.cm50$case=='3',]; baycm.c4.50<-bay.cm50[bay.cm50$case=='4',]
baycm.c5.50<-bay.cm50[bay.cm50$case=='5',]; baycm.c6.50<-bay.cm50[bay.cm50$case=='6',]

baym.c1.50<-bay.m50[bay.m50$case=='1',]; baym.c2.50<-bay.m50[bay.m50$case=='2',]
baym.c3.50<-bay.m50[bay.m50$case=='3',]; baym.c4.50<-bay.m50[bay.m50$case=='4',]
baym.c5.50<-bay.m50[bay.m50$case=='5',]; baym.c6.50<-bay.m50[bay.m50$case=='6',]

##################################### N=75 #####################################
bay.m75<-read.table(paste(PATH3, "table(bayM).txt",sep=""))
colnames(bay.m75)<-c("loop","para","case","real data","Est","SE","LCI","UCI","Len","CP")

bay.cm75<-read.table(paste(PATH3, "table(bayCM).txt",sep=""))
colnames(bay.cm75)<-c("loop","para","case","real data","Est","SE","LCI","UCI","Len","CP")

baycm.c1.75<-bay.cm75[bay.cm75$case=='1',]; baycm.c2.75<-bay.cm75[bay.cm75$case=='2',]
baycm.c3.75<-bay.cm75[bay.cm75$case=='3',]; baycm.c4.75<-bay.cm75[bay.cm75$case=='4',]
baycm.c5.75<-bay.cm75[bay.cm75$case=='5',]; baycm.c6.75<-bay.cm75[bay.cm75$case=='6',]

baym.c1.75<-bay.m75[bay.m75$case=='1',]; baym.c2.75<-bay.m75[bay.m75$case=='2',]
baym.c3.75<-bay.m75[bay.m75$case=='3',]; baym.c4.75<-bay.m75[bay.m75$case=='4',]
baym.c5.75<-bay.m75[bay.m75$case=='5',]; baym.c6.75<-bay.m75[bay.m75$case=='6',]

##################################### N=100 #####################################
bay.m100<-read.table(paste(PATH4, "table(bayM).txt",sep=""))
colnames(bay.m100)<-c("loop","para","case","real data","Est","SE","LCI","UCI","Len","CP")

bay.cm100<-read.table(paste(PATH4, "table(bayCM).txt",sep=""))
colnames(bay.cm100)<-c("loop","para","case","real data","Est","SE","LCI","UCI","Len","CP")

baycm.c1.100<-bay.cm100[bay.cm100$case=='1',]; baycm.c2.100<-bay.cm100[bay.cm100$case=='2',]
baycm.c3.100<-bay.cm100[bay.cm100$case=='3',]; baycm.c4.100<-bay.cm100[bay.cm100$case=='4',]
baycm.c5.100<-bay.cm100[bay.cm100$case=='5',]; baycm.c6.100<-bay.cm100[bay.cm100$case=='6',]

baym.c1.100<-bay.m100[bay.m100$case=='1',]; baym.c2.100<-bay.m100[bay.m100$case=='2',]
baym.c3.100<-bay.m100[bay.m100$case=='3',]; baym.c4.100<-bay.m100[bay.m100$case=='4',]
baym.c5.100<-bay.m100[bay.m100$case=='5',]; baym.c6.100<-bay.m100[bay.m100$case=='6',]

MSE.table = function(baycm,baym){
  bayCM.MSE = cbind((baycm[which(baycm$para=='b1'),5]-par.true[1])^2, 
                    (baycm[which(baycm$para=='b2'),5]-par.true[2])^2,
                    (baycm[which(baycm$para=='b3'),5]-par.true[3])^2,
                    (baycm[which(baycm$para=='b4'),5]-par.true[4])^2,
                    (baycm[which(baycm$para=='b5'),5]-par.true[5])^2,
                    (baycm[which(baycm$para=='d11'),5]-par.true[6])^2,
                    (baycm[which(baycm$para=='d21'),5]-par.true[7])^2,
                    (baycm[which(baycm$para=='d22'),5]-par.true[8])^2,
                    (baycm[which(baycm$para=='sigma11'),5]-par.true[9])^2,
                    (baycm[which(baycm$para=='sigma21'),5]-par.true[10])^2,
                    (baycm[which(baycm$para=='sigma22'),5]-par.true[11])^2)
  
  bayM.MSE = cbind((baym[which(baym$para=='b1'),5]-par.true[1])^2, 
                   (baym[which(baym$para=='b2'),5]-par.true[2])^2,
                   (baym[which(baym$para=='b3'),5]-par.true[3])^2,
                   (baym[which(baym$para=='b4'),5]-par.true[4])^2,
                   (baym[which(baym$para=='b5'),5]-par.true[5])^2,
                   (baym[which(baym$para=='d11'),5]-par.true[6])^2,
                   (baym[which(baym$para=='d21'),5]-par.true[7])^2,
                   (baym[which(baym$para=='d22'),5]-par.true[8])^2,
                   (baym[which(baym$para=='sigma11'),5]-par.true[9])^2,
                   (baym[which(baym$para=='sigma21'),5]-par.true[10])^2,
                   (baym[which(baym$para=='sigma22'),5]-par.true[11])^2)
  
  Bay.CM<-round(as.vector(bayCM.MSE), 4)
  Bay.M<-round(as.vector(bayM.MSE), 4)
  return(list(Bay.M=Bay.M, Bay.CM=Bay.CM))
}

c1.25<-MSE.table(baycm.c1.25, baym.c1.25); c2.25<-MSE.table(baycm.c2.25, baym.c2.25)
c3.25<-MSE.table(baycm.c3.25, baym.c3.25); c4.25<-MSE.table(baycm.c4.25, baym.c4.25)
c5.25<-MSE.table(baycm.c5.25, baym.c5.25); c6.25<-MSE.table(baycm.c6.25, baym.c6.25)

c1.50<-MSE.table(baycm.c1.50, baym.c1.50); c2.50<-MSE.table(baycm.c2.50, baym.c2.50)
c3.50<-MSE.table(baycm.c3.50, baym.c3.50); c4.50<-MSE.table(baycm.c4.50, baym.c4.50)
c5.50<-MSE.table(baycm.c5.50, baym.c5.50); c6.50<-MSE.table(baycm.c6.50, baym.c6.50)

c1.75<-MSE.table(baycm.c1.75, baym.c1.75); c2.75<-MSE.table(baycm.c2.75, baym.c2.75)
c3.75<-MSE.table(baycm.c3.75, baym.c3.75); c4.75<-MSE.table(baycm.c4.75, baym.c4.75)
c5.75<-MSE.table(baycm.c5.75, baym.c5.75); c6.75<-MSE.table(baycm.c6.75, baym.c6.75)

c1.100<-MSE.table(baycm.c1.100, baym.c1.100); c2.100<-MSE.table(baycm.c2.100, baym.c2.100)
c3.100<-MSE.table(baycm.c3.100, baym.c3.100); c4.100<-MSE.table(baycm.c4.100, baym.c4.100)
c5.100<-MSE.table(baycm.c5.100, baym.c5.100); c6.100<-MSE.table(baycm.c6.100, baym.c6.100)

### Case 1: MAR + Censor 10% ###
Bay.MAR10<-data.frame(MSE=c(c1.25$Bay.M,c1.50$Bay.M,c1.75$Bay.M,c1.100$Bay.M,
                            c1.25$Bay.CM,c1.50$Bay.CM,c1.75$Bay.CM,c1.100$Bay.CM),
                            model=rep(c("AM","CM"),each=4400),
                            N=rep(rep(c(25,50,75,100),each=1100),2),stringsAsFactors = FALSE)
colnames(Bay.MAR10)<-c("MSE", "model", "N")

### Case 3: MAR + Censor 30% ###
Bay.MAR30<-data.frame(MSE=c(c3.25$Bay.M,c3.50$Bay.M,c3.75$Bay.M,c3.100$Bay.M,
                            c3.25$Bay.CM,c3.50$Bay.CM,c3.75$Bay.CM,c3.100$Bay.CM),
                      model=rep(c("AM","CM"),each=4400),
                      N=rep(rep(c(25,50,75,100),each=1100),2),stringsAsFactors = FALSE)
colnames(Bay.MAR30)<-c("MSE", "model", "N")

### Case 5: MAR + Censor 50% ###
Bay.MAR50<-data.frame(MSE=c(c5.25$Bay.M,c5.50$Bay.M,c5.75$Bay.M,c5.100$Bay.M,
                            c5.25$Bay.CM,c5.50$Bay.CM,c5.75$Bay.CM,c5.100$Bay.CM),
                      model=rep(c("AM","CM"),each=4400),
                      N=rep(rep(c(25,50,75,100),each=1100),2),stringsAsFactors = FALSE)
colnames(Bay.MAR50)<-c("MSE", "model", "N")

### Case 2: MCAR + Censor 10% ###
Bay.MCAR10<-data.frame(MSE=c(c2.25$Bay.M,c2.50$Bay.M,c2.75$Bay.M,c2.100$Bay.M,
                            c2.25$Bay.CM,c2.50$Bay.CM,c2.75$Bay.CM,c2.100$Bay.CM),
                      model=rep(c("AM","CM"),each=4400),
                      N=rep(rep(c(25,50,75,100),each=1100),2),stringsAsFactors = FALSE)
colnames(Bay.MCAR10)<-c("MSE", "model", "N")

### Case 4: MCAR + Censor 30% ###
Bay.MCAR30<-data.frame(MSE=c(c4.25$Bay.M,c4.50$Bay.M,c4.75$Bay.M,c4.100$Bay.M,
                            c4.25$Bay.CM,c4.50$Bay.CM,c4.75$Bay.CM,c4.100$Bay.CM),
                      model=rep(c("AM","CM"),each=4400),
                      N=rep(rep(c(25,50,75,100),each=1100),2),stringsAsFactors = FALSE)
colnames(Bay.MCAR30)<-c("MSE", "model", "N")

### Case 6: MCAR + Censor 50% ###
Bay.MCAR50<-data.frame(MSE=c(c6.25$Bay.M,c6.50$Bay.M,c6.75$Bay.M,c6.100$Bay.M,
                            c6.25$Bay.CM,c6.50$Bay.CM,c6.75$Bay.CM,c6.100$Bay.CM),
                      model=rep(c("AM","CM"),each=4400),
                      N=rep(rep(c(25,50,75,100),each=1100),2),stringsAsFactors = FALSE)
colnames(Bay.MCAR50)<-c("MSE", "model", "N")

################################################################################
# --------------------------- Figure 1 for Simulation --------------------------
################################################################################
postscript(paste(PATH, 'results/SIMfig1.eps', sep=''), height=4, width=12)
layout(matrix(c(1,1,1,8,8,8,2,3,4,9,10,11,5,6,7,12,13,14), 3, 6, byrow=T),
       widths=c(6,4.5,4.5,6,4.5,4.5), heights=c(0.75,0.75,4))
N = c(25,50,75,100)
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
boxplot(MSE ~ model+N, data = Bay.MAR10,ylim=c(0,2.7),las=1,cex.lab=1.6,cex.axis=1.2,
        at = c(1,2,4,5,7,8,10,11), col = c(0,"palevioletred"), border=c("blue","red") ,outline=FALSE,
        xaxt="n",xlab="Sample Size", ylab = expression(MSE~(hat(theta))))
axis(1,labels=c(25,50,75,100),at=c(1.5,4.5,7.5,10.5),cex.axis=1.2)
legend("topright", fill = c(0,"palevioletred"), border=c("blue","red"), legend = c("MNLMM","MNLMM-CM"), cex=1.2, bty='n')

par(mar=c(4, 0, 0, 0), cex.lab=2)
boxplot(MSE ~ model+N, data = Bay.MAR30,ylim=c(0,2.7), yaxt='n',cex.lab=1.6,
        at = c(1,2,4,5,7,8,10,11), col = c(0,"palevioletred"), border=c("blue","red"),outline=FALSE,
        xaxt="n",xlab="Sample Size")
axis(1,labels=c(25,50,75,100),at=c(1.5,4.5,7.5,10.5),cex.axis=1.2)
axis(2,labels=F,at=seq(0,2.5,0.5))
legend("topright", fill = c(0,"palevioletred"), border=c("blue","red"), legend = c("MNLMM","MNLMM-CM"), cex=1.2, bty='n')

par(mar=c(4, 0, 0, 0.5), cex.lab=2)
boxplot(MSE ~ model+N, data = Bay.MAR50,ylim=c(0,2.7), yaxt='n',cex.lab=1.6,
        at = c(1,2,4,5,7,8,10,11), col = c(0,"palevioletred"), border=c("blue","red"),outline=FALSE,
        xaxt="n",xlab="Sample Size")
axis(1,labels=c(25,50,75,100),at=c(1.5,4.5,7.5,10.5),cex.axis=1.2)
axis(2,labels=F,at=seq(0,2.5,0.5))
legend("topright", fill = c(0,"palevioletred"), border=c("blue","red"), legend = c("MNLMM","MNLMM-CM"), cex=1.2, bty='n')

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
boxplot(MSE ~ model+N, data = Bay.MCAR10,ylim=c(0,2.7),las=1,cex.lab=1.6,cex.axis=1.2,
        at = c(1,2,4,5,7,8,10,11), col = c(0,"palevioletred"), border=c("blue","red"),outline=FALSE,
        xaxt="n",xlab="Sample Size", ylab = expression(MSE~(hat(theta))))
axis(1,labels=c(25,50,75,100),at=c(1.5,4.5,7.5,10.5),cex.axis=1.2)
legend("topright", fill = c(0,"palevioletred"), border=c("blue","red"), legend = c("MNLMM","MNLMM-CM"), cex=1.2, bty='n')

par(mar=c(4, 0, 0, 0), cex.lab=2)
boxplot(MSE ~ model+N, data = Bay.MCAR30,ylim=c(0,2.7), yaxt='n',cex.lab=1.6,
        at = c(1,2,4,5,7,8,10,11), col = c(0,"palevioletred"), border=c("blue","red"),outline=FALSE,
        xaxt="n",xlab="Sample Size")
axis(1,labels=c(25,50,75,100),at=c(1.5,4.5,7.5,10.5),cex.axis=1.2)
axis(2,labels=F,at=seq(0,2.5,0.5))
legend("topright", fill = c(0,"palevioletred"), border=c("blue","red"), legend = c("MNLMM","MNLMM-CM"), cex=1.2, bty='n')

par(mar=c(4, 0, 0, 0.5), cex.lab=2)
boxplot(MSE ~ model+N, data = Bay.MCAR50,ylim=c(0,2.7), yaxt='n',cex.lab=1.6,
        at = c(1,2,4,5,7,8,10,11), col = c(0,"palevioletred"), border=c("blue","red"),outline=FALSE,
        xaxt="n",xlab="Sample Size")
axis(1,labels=c(25,50,75,100),at=c(1.5,4.5,7.5,10.5),cex.axis=1.2)
axis(2,labels=F,at=seq(0,2.5,0.5))
legend("topright", fill = c(0,"palevioletred"), border=c("blue","red"), legend = c("MNLMM","MNLMM-CM"), cex=1.2, bty='n')
dev.off()
