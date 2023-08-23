# Reproduct Supplementary Figure S.2
library(nlme)
A5055 = read.table(paste(PATH,'data/A5055data.txt',sep=""), header=T)
cen.rna = ifelse(A5055$rna<=50,1,0)
A5055 = data.frame(cbind(A5055, cen.rna, new.rna=A5055$rna))
A5055$new.rna[A5055$new.rna<=50]=50
Time = A5055$day/7
lrna = round(log10(A5055$new.rna),2)
cen.lrna = cen.rna
lcd4 = log10(A5055$cd4)
lcd8 = log10(A5055$cd8)
cd4.cd8 = A5055$cd4/A5055$cd8
na.cd4.cd8 = ifelse(is.na(cd4.cd8)==T,1,0)
A5055 = data.frame(cbind(A5055,Time=Time,lrna=lrna,cen.lrna=cen.lrna,lcd4=lcd4,
                         lcd8=lcd8,cd4.cd8=cd4.cd8,na.cd4.cd8))
r = 2
censored = 1.7

naid = unique(A5055$Subject[which(A5055$na.cd4.cd8 == 1)])

############################# IDV-RTV = 0 Group ################################
Trt1 = unique(A5055$Subject[which(A5055$arm == 1)]); n1 = length(Trt1); # n1 = 22
armA = A5055[which(A5055$arm == 1), ]
nj1 = numeric(n1)
for(j in 1:n1) nj1[j]=sum(armA$Subject==Trt1[j])
armAG = data.frame(cbind(A5055[which(A5055$arm==1),2], armA))
colnames(armAG) = c('ID','arm','Subject','patid','day','logrna','rna','cd4','cd8','cen.rna', 
                    'new.rna','Time','lrna','cen.lrna','lcd4','lcd8','cd4.cd8','na.cd4.cd8')

### censored = 1 ##
Trt1.C = unique(armAG$ID[which(armAG$cen.rna == 1)]); n1.C = length(Trt1.C); # n1 = 22
nj1.C = numeric(n1.C)
for(j in 1:n1.C) nj1.C[j]=sum(armA$Subject==Trt1.C[j])
armAG.C = armAG[which(armAG$ID==1 | armAG$ID==2 | armAG$ID==3 | armAG$ID==5 | 
                        armAG$ID==9 | armAG$ID==10 | armAG$ID==11 | armAG$ID==12 |
                        armAG$ID==13 | armAG$ID==14 | armAG$ID==15 | armAG$ID==16 | 
                        armAG$ID==17 | armAG$ID==19 | armAG$ID==20 | armAG$ID==22),]

### CD4/CD8 = NA ##
Trt1.NA = unique(armAG$ID[which(armAG$na.cd4.cd8 == 1)]); n1.NA = length(Trt1.NA); # n1 = 22
nj1.NA = numeric(n1.NA)
for(j in 1:n1.NA) nj1.NA[j]=sum(armAG$ID==Trt1.NA[j])
armAG.NA = armAG[-which(armAG$ID==6),]

############################# IDV-RTV = 1 Group ################################
Trt2 = unique(A5055$Subject[which(A5055$arm == 2)]); n2 = length(Trt2); # n2 = 22
armB = A5055[which(A5055$arm == 2), ]
nj2 = numeric(n2)
for(j in 1:n2) nj2[j]=sum(armB$Subject==Trt2[j])
armBG = data.frame(cbind(A5055[which(A5055$arm==2),2], armB))
colnames(armBG) = c('ID','arm','Subject','patid','day','logrna','rna','cd4','cd8','cen.rna', 
                    'new.rna','Time','lrna','cen.lrna','lcd4','lcd8','cd4.cd8','na.cd4.cd8')

Trt2.C = unique(armBG$ID[which(armBG$cen.rna == 1)]); n2.C = length(Trt2.C); # n2 = 22
nj2.C = numeric(n2.C)
for(j in 1:n2.C) nj2.C[j]=sum(armBG$ID==Trt2.C[j])
armBG.C = armBG[-which(armBG$ID==23 | armBG$ID==28 | armBG$ID==31 | armBG$ID==33 |
                         armBG$ID==38 | armBG$ID==41),]

### CD4/CD8 = NA ##
Trt2.NA = unique(armBG$ID[which(armBG$na.cd4.cd8 == 1)]); n2.NA = length(Trt2.NA); # n2 = 22
nj2.NA = numeric(n2.NA)
for(j in 1:n2.NA) nj2.NA[j]=sum(armBG$ID==Trt2.NA[j])
armBG.NA = armBG[-which(armBG$ID == 39),]

##################################### plot #####################################
postscript(paste(PATH, 'results/figS2.eps', sep=''), width=10, height=13)
layout(matrix(c(1:6), 3, 2), widths=c(6,5.2), heights=c(1.5, 8, 9.5))
############### IDV-RTV = 0 Group ###############
par(mar=c(0, 5.5, 0.5, 0), cex.lab=2)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', ylab='')
mtext('IDV-RTV = 0', 1, line=-2, cex=1.5, font=3)

### log(RNA)
par(mar=c(0, 5.5, 0, 0), cex.lab=2)
plot(1:10, 1:10, xlim=range(-4,207), ylim=c(0.5,6.2), type='n', xlab='', 
    ylab=expression(log[10](RNA)), xaxt='n', yaxt='n', las=1,cex.lab = 1.5)
abline(h=1:6,v=c(0,14,28,56,84,112,140,168,196), col = "gray", lty = 3, lwd=0.5)
abline(h=censored,col="blue",lty=2,lwd=2)
axis(1, at=c(0,14,28,56,84,112,140,168,196), labels = F)
axis(2, at=seq(0, 6, 1), labels = (seq(0, 6, 1)), las=1,cex.axis=1.2)

lines(armA$day[1:nj1[1]], armA[1:nj1[1],]$logrna, lwd=1, col="black", type="n", pch=16, cex=0.8)
text(armAG.C$day[nj1.C[1]]+5,armAG.C[nj1.C[1],]$logrna+0.05,armAG.C[nj1.C[1],]$ID,col="red")
for(i in 1:n1) {
  lines(armA$day[armA$Subject==i], armA$logrna[armA$Subject==i], lwd=1, col="black", type="o", pch=16, cex=0.8)     
}

for(i in c(2,3,5,7,9)){
  text(armAG.C$day[sum(nj1.C[0:(i+2)])]+7,armAG.C[sum(nj1.C[0:(i+2)]),]$logrna+0.05,armAG.C[sum(nj1.C[0:(i+2)]),]$ID, col="red",cex=0.9)
}
for(i in 13){
  text(armAG.C$day[sum(nj1.C[0:(i+2)])]+8,armAG.C[sum(nj1.C[0:(i+2)]),]$logrna-0.005,armAG.C[sum(nj1.C[0:(i+2)]),]$ID, col="red",cex=0.9)
}
for(i in 14){
  text(armAG.C$day[sum(nj1.C[0:(i+2)])]+8,armAG.C[sum(nj1.C[0:(i+2)]),]$logrna+0.05,armAG.C[sum(nj1.C[0:(i+2)]),]$ID, col="red",cex=0.9)
}
for(i in 8){
  text(armAG.C$day[sum(nj1.C[0:(i+2)])]+7,armAG.C[sum(nj1.C[0:(i+2)]),]$logrna,armAG.C[sum(nj1.C[0:(i+2)]),]$ID, col="red",cex=0.9)
}
for(i in 4){
  text(armAG.C$day[sum(nj1.C[0:(i+2)])]+6,armAG.C[sum(nj1.C[0:(i+2)]),]$logrna-0.02,armAG.C[sum(nj1.C[0:(i+2)]),]$ID, col="red",cex=0.9)
}
for(i in 12){
  text(armAG.C$day[sum(nj1.C[0:(i+2)])],armAG.C[sum(nj1.C[0:(i+2)]),]$logrna-0.15,armAG.C[sum(nj1.C[0:(i+2)]),]$ID, col="red",cex=0.9)
}
for(i in c(0,1)){
  text(armAG.C$day[(sum(nj1[0:(i+1)])+1)]-6,armAG.C[(sum(nj1[0:(i+1)])+1),]$logrna+0.05,armAG.C[sum(nj1.C[0:(i+2)]),]$ID, col="red",cex=0.9)
}
for(i in 6){
  text(armAG.C$day[(sum(nj1[0:(i+1)])+1)]+5,armAG.C[(sum(nj1[0:(i+1)])+1),]$logrna+0.05,armAG.C[sum(nj1.C[0:(i+2)]),]$ID, col="red",cex=0.9)
}
for(i in 10){
  text(armAG.C$day[(sum(nj1[0:(i+1)])+1)]+5,armAG.C[(sum(nj1[0:(i+1)])+1),]$logrna+0.05,armAG.C[sum(nj1.C[0:(i+2)]),]$ID, col="red",cex=0.9)
}
for(i in 11){
  text(armAG.C$day[(sum(nj1[0:(i+1)])+1)]-7,armAG.C[(sum(nj1[0:(i+1)])+1),]$logrna+0.05,armAG.C[sum(nj1.C[0:(i+2)]),]$ID, col="red",cex=0.9)
}

#### CD4/CD8
par(mar=c(5, 5.5, 0, 0), cex.lab=2)
plot(1:10, 1:10, xlim=range(-4,207), ylim=c(0, 2.1), type='n', xlab='Day', ylab="CD4/CD8",
     xaxt='n', yaxt='n', las=1,cex.lab = 1.5)
abline(h=seq(0, 2, 0.5),v=c(0,14,28,56,84,112,140,168,196), col = "gray", lty = 3, lwd=0.5)
axis(1, at=c(0,14,28,56,84,112,140,168,196), labels = (c(0,14,28,56,84,112,140,168,196)),cex.axis=1.2)
axis(2, at=seq(0, 2, 0.5), labels = (seq(0, 2, 0.5)), las=1,cex.axis=1.2)
lines(armA$day[1:nj1[1]],armA[1:nj1[1],]$cd4.cd8, lwd=1, col="black", type="n", pch=16, cex=0.8)
text(armAG.NA$day[nj1.NA[1]]+5,armAG.NA[nj1.NA[1],]$cd4.cd8+0.01,armAG.NA[nj1.NA[1],]$ID,col="red",cex=0.9)
for(i in 1:n1){
  lines(armA$day[armA$Subject==i], armA$cd4.cd8[armA$Subject==i], lwd=1, col="black", type="o", pch=16, cex=0.8)     
}
lines(armA$day[armA$Subject==6], armA$cd4.cd8[armA$Subject==6], lwd=1, col="black", type="o", pch=16, cex=0.8)   

for(i in c(0,2,5,12,17:19)){
  text(armAG.NA$day[sum(nj1.NA[0:(i+2)])]+7,armAG.NA[sum(nj1.NA[0:(i+2)]),]$cd4.cd8+0.02,armAG.NA[sum(nj1.NA[0:(i+2)]),]$ID, col="red",cex=0.9)
}
for(i in c(1,3,4,6,7)){
  text(armAG.NA$day[sum(nj1.NA[0:(i+2)])]+5,armAG.NA[sum(nj1.NA[0:(i+2)]),]$cd4.cd8+0.02,armAG.NA[sum(nj1.NA[0:(i+2)]),]$ID, col="red",cex=0.9)
}
for(i in 8){
  text(armAG.NA$day[sum(nj1.NA[0:(i+2)])]+7,armAG.NA[sum(nj1.NA[0:(i+2)]),]$cd4.cd8,armAG.NA[sum(nj1.NA[0:(i+2)]),]$ID, col="red",cex=0.9)
}
for(i in 11){
  text(armAG.NA$day[sum(nj1.NA[0:(i+2)])],armAG.NA[sum(nj1.NA[0:(i+2)]),]$cd4.cd8-0.1,armAG.NA[sum(nj1.NA[0:(i+2)]),]$ID, col="red",cex=0.9)
}
for(i in 15){
  text(armAG.NA$day[sum(nj1.NA[0:(i+2)])]+7,armAG.NA[sum(nj1.NA[0:(i+2)]),]$cd4.cd8-0.02,armAG.NA[sum(nj1.NA[0:(i+2)]),]$ID, col="red",cex=0.9)
}
for(i in 13){
  text(armAG.NA$day[sum(nj1.NA[0:(i+2)])],armAG.NA[sum(nj1.NA[0:(i+2)]),]$cd4.cd8+0.06,armAG.NA[sum(nj1.NA[0:(i+2)]),]$ID, col="red",cex=0.9)
}
for(i in 16){
  text(armAG.NA$day[sum(nj1.NA[0:(i+2)])],armAG.NA[sum(nj1.NA[0:(i+2)]),]$cd4.cd8+0.08,armAG.NA[sum(nj1.NA[0:(i+2)]),]$ID, col="red",cex=0.9)
}
for(i in 14){
  text(armAG.NA$day[(sum(nj1.NA[0:(i+1)])+1)]-7,armAG.NA[(sum(nj1.NA[0:(i+1)])+1),]$cd4.cd8+0.02,armAG.NA[sum(nj1.NA[0:(i+2)]),]$ID, col="red",cex=0.9)
}
for(i in 9){
  text(armAG.NA$day[(sum(nj1.NA[0:(i+1)])+1)]-7,armAG.NA[(sum(nj1.NA[0:(i+1)])+1),]$cd4.cd8+0.02,armAG.NA[sum(nj1.NA[0:(i+2)]),]$ID, col="red",cex=0.9)
}
for(i in 10){
  text(armAG.NA$day[(sum(nj1.NA[0:(i+1)])+1)]+5,armAG.NA[(sum(nj1.NA[0:(i+1)])+1),]$cd4.cd8+0.02,armAG.NA[sum(nj1.NA[0:(i+2)]),]$ID, col="red",cex=0.9)
}

############### IDV-RTV = 1 Group ###############
par(mar=c(0, 0, 0.5, 0.5), cex.lab=2)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', ylab='')
mtext('IDV-RTV = 1', 1, line=-2.1, cex=1.5, font=3)
### log(RNA)
par(mar=c(0, 0, 0, 0.5), cex.lab=2)
plot(1:10, 1:10, xlim=range(-4,207), ylim=c(0.5,6.2), type='n', xlab='', ylab='', xaxt='n', yaxt='n', las=1)
abline(h=1:6,v=c(0,14,28,56,84,112,140,168,196), col = "gray", lty = 3, lwd=0.5)
abline(h=censored,col="blue",lty=2,lwd=2)
axis(1, at=c(0,14,28,56,84,112,140,168,196), labels = F)
axis(2, at=seq(0, 6, 1), labels = F)
lines(armB$day[1:nj2[1]],armB[1:nj2[1],]$logrna, lwd=1, col="black", type="n", pch=16, cex=0.8)
text(armBG.C$day[nj2.C[1]]+7,armBG.C[nj2.C[1],]$logrna+0.05,armBG.C[nj2.C[1],]$ID,col="red")
for(i in n1+(1: n2)){
  lines(armB$day[armB$Subject==i], armB$logrna[armB$Subject==i], lwd=1, col="black", type="o", pch=16, cex=0.8)     
}

for(i in 8){
  text(armBG.C$day[sum(nj2.C[0:(i+2)])]+7,armBG.C[sum(nj2.C[0:(i+2)]),]$logrna+0.05,armBG.C[sum(nj2.C[0:(i+2)]),]$ID, col="red")
}
for(i in 13){
  text(armBG.C$day[sum(nj2.C[0:(i+2)])]+7,armBG.C[sum(nj2.C[0:(i+2)]),]$logrna,armBG.C[sum(nj2.C[0:(i+2)]),]$ID, col="red")
}
for(i in c(12,14)){
  text(armBG.C$day[sum(nj2.C[0:(i+2)])]+7,armBG.C[sum(nj2.C[0:(i+2)]),]$logrna+0.05,armBG.C[sum(nj2.C[0:(i+2)]),]$ID, col="red")
}
for(i in 11){
  text(armBG.C$day[sum(nj2.C[0:(i+2)])]+7,armBG.C[sum(nj2.C[0:(i+2)]),]$logrna,armBG.C[sum(nj2.C[0:(i+2)]),]$ID, col="red")
}
for(i in c(1,4)){
  text(armBG.C$day[sum(nj2.C[0:(i+2)])]+7,armBG.C[sum(nj2.C[0:(i+2)]),]$logrna+0.1,armBG.C[sum(nj2.C[0:(i+2)]),]$ID, col="red")
}
for(i in c(0,2)){
  text(armBG.C$day[sum(nj2.C[0:(i+2)])]+7,armBG.C[sum(nj2.C[0:(i+2)]),]$logrna+0.02,armBG.C[sum(nj2.C[0:(i+2)]),]$ID, col="red")
}
for(i in 9){
  text(armBG.C$day[sum(nj2.C[0:(i+2)])]-4,armBG.C[sum(nj2.C[0:(i+2)]),]$logrna+0.25,armBG.C[sum(nj2.C[0:(i+2)]),]$ID, col="red")
}
for(i in 7){
  text(armBG.C$day[sum(nj2.C[0:(i+2)])],armBG.C[sum(nj2.C[0:(i+2)]),]$logrna-0.1,armBG.C[sum(nj2.C[0:(i+2)]),]$ID, col="red")
}
for(i in c(3,5,6,10)){
  text(armBG.C$day[(sum(nj2.C[0:(i+1)])+1)]-7,armBG.C[(sum(nj2.C[0:(i+1)])+1),]$logrna+0.05,armBG.C[sum(nj2.C[0:(i+2)]),]$ID, col="red",cex=0.9)
}

#### CD4/CD8
par(mar=c(5, 0, 0, 0.5), cex.lab=2)
plot(1:10, 1:10, xlim=range(-4,207), ylim=c(0, 2.1), type='n', xlab='Day', ylab='', xaxt='n', yaxt='n', las=1,cex.lab = 1.5)
abline(h=seq(0, 2, 0.5),v=c(0,14,28,56,84,112,140,168,196), col = "gray", lty = 3, lwd=0.5)
axis(1, at=c(0,14,28,56,84,112,140,168,196), labels = (c(0,14,28,56,84,112,140,168,196)),cex.axis=1.2)
axis(2, at=seq(0, 2, 0.5), F)
lines(armB$day[1:nj2[1]],armB[1:nj2[1],]$cd4.cd8, lwd=1, col="black", type="o", pch=16, cex=0.8)
text(armBG.NA$day[nj2.NA[1]]+5,armBG.NA[nj2.NA[1],]$cd4.cd8+0.01,armBG.NA[nj2.NA[1],]$ID,col="red")
for(i in n1+(1:n2)) {
  lines(armB$day[armB$Subject==i], armB$cd4.cd8[armB$Subject==i], lwd=1, col="black", type="o", pch=16, cex=0.8)     
}
lines(armB$day[armB$Subject==39], armB$cd4.cd8[armB$Subject==39], lwd=1, col="black", type="o", pch=16, cex=0.8)  
  
for(i in c(1,2,6:7,9:19)){
  text(armBG.NA$day[sum(nj2.NA[0:(i+2)])]+7,armBG.NA[sum(nj2.NA[0:(i+2)]),]$cd4.cd8+0.01,armBG.NA[sum(nj2.NA[0:(i+2)]),]$ID, col="red")
}
for(i in c(0,8)){
  text(armBG.NA$day[sum(nj2.NA[0:(i+2)])],armBG.NA[sum(nj2.NA[0:(i+2)]),]$cd4.cd8-0.05,armBG.NA[sum(nj2.NA[0:(i+2)]),]$ID, col="red")
}
for(i in 5){
  text(armBG.NA$day[sum(nj2.NA[0:(i+2)])],armBG.NA[sum(nj2.NA[0:(i+2)]),]$cd4.cd8+0.06,armBG.NA[sum(nj2.NA[0:(i+2)]),]$ID, col="red")
}
for(i in 3){
  text(armBG.NA$day[(sum(nj2.NA[0:(i+1)])+1)]+7,armBG.NA[(sum(nj2.NA[0:(i+1)])+1),]$cd4.cd8+0.02,armBG.NA[sum(nj2.NA[0:(i+2)]),]$ID, col="red",cex=0.9)
}
for(i in 4){
  text(armBG.NA$day[(sum(nj2.NA[0:(i+1)])+1)]-7,armBG.NA[(sum(nj2.NA[0:(i+1)])+1),]$cd4.cd8+0.02,armBG.NA[sum(nj2.NA[0:(i+2)]),]$ID, col="red",cex=0.9)
}
dev.off()
