library(nlme)
A5055 = read.table(paste(PATH,'data/A5055data.txt', sep=""), header=T)
cen.rna = ifelse(A5055$rna<=50,1,0)
A5055 = data.frame(cbind(A5055, cen.rna, new.rna=A5055$rna))
A5055$new.rna[A5055$new.rna<=50]=50
Time = A5055$day/7
lrna = round(log10(A5055$new.rna),2)
cen.lrna = cen.rna
lcd4 = log10(A5055$cd4)
lcd8 = log10(A5055$cd8)
cd4.cd8=A5055$cd4/A5055$cd8
A5055 = data.frame(cbind(A5055,Time=Time,lrna=lrna,cen.lrna=cen.lrna,lcd4=lcd4,
                         lcd8=lcd8,cd4.cd8=cd4.cd8))
r = 2

################################## ARM-A Group #################################
Trt1 = unique(A5055$Subject[which(A5055$arm == 1)]); n1 = length(Trt1); 
armA = A5055[which(A5055$arm == 1), ]
nj1 = numeric(n1)
for(j in 1:n1) nj1[j]=sum(armA$Subject==Trt1[j])
armAG = data.frame(cbind(rep(1:n1, nj1), armA))
colnames(armAG) = c('ID','arm','Subject','patid','day','logrna','rna','cd4','cd8','cd4.cd8')

################################# ARM-B Group #################################
Trt2 = unique(A5055$Subject[which(A5055$arm == 2)]); n2 = length(Trt2); 
armB = A5055[which(A5055$arm == 2), ]
nj2 = numeric(n2)
for(j in 1:n2) nj2[j]=sum(armB$Subject==Trt2[j])
armBG = data.frame(cbind(rep(1:n2, nj2), armB))
colnames(armBG) = c('ID','arm','Subject','patid','day','logrna','rna','cd4','cd8','cd4.cd8')

################################################################################
# -------------------------- Supplementary Table S2-----------------------------
################################################################################

################################## IDV-RTV=0 ###################################
### Group weeks
W0.armA<-armA[which(armA$day==0),]
W14.armA<-armA[which(armA$day<=14 & armA$day>0),]
W28.armA<-armA[which(armA$day<=28 & armA$day>14),]
W56.armA<-armA[which(armA$day<=56 & armA$day>28),]
W84.armA<-armA[which(armA$day<=84 & armA$day>56),]
W112.armA<-armA[which(armA$day<=112 & armA$day>84),]
W140.armA<-armA[which(armA$day<=140 & armA$day>112),]
W168.armA<-armA[which(armA$day>140),]

### log(RNA) ###
# mean
armA.rna.mean<-c(mean(W0.armA$logrna),mean(W14.armA$logrna),mean(W28.armA$logrna),
                 mean(W56.armA$logrna),mean(W84.armA$logrna),mean(W112.armA$logrna),
                 mean(W140.armA$logrna),mean(W168.armA$logrna))
# sd
armA.rna.sd<-c(sd(W0.armA$logrna),sd(W14.armA$logrna),sd(W28.armA$logrna),sd(W56.armA$logrna),
               sd(W84.armA$logrna),sd(W112.armA$logrna),sd(W140.armA$logrna),sd(W168.armA$logrna))

# censored rate
armA.cen<-c(sum(W0.armA$cen.lrna)/length(W0.armA$cen.lrna),sum(W14.armA$cen.lrna)/length(W14.armA$cen.lrna),
            sum(W28.armA$cen.lrna)/length(W28.armA$cen.lrna),sum(W56.armA$cen.lrna)/length(W56.armA$cen.lrna),
            sum(W84.armA$cen.lrna)/length(W84.armA$cen.lrna),sum(W112.armA$cen.lrna)/length(W112.armA$cen.lrna),
            sum(W140.armA$cen.lrna)/length(W140.armA$cen.lrna),sum(W168.armA$cen.lrna)/length(W168.armA$cen.lrna))*100

### CD4/CD8 ###
# mean
armA.CD.mean<-c(mean(W0.armA$cd4.cd8,na.rm=T),mean(W14.armA$cd4.cd8,na.rm=T),mean(W28.armA$cd4.cd8,na.rm=T),
                mean(W56.armA$cd4.cd8,na.rm=T),mean(W84.armA$cd4.cd8,na.rm=T),mean(W112.armA$cd4.cd8,na.rm=T),
                mean(W140.armA$cd4.cd8,na.rm=T),mean(W168.armA$cd4.cd8,na.rm=T))
# sd
armA.CD.sd<-c(sd(W0.armA$cd4.cd8,na.rm=T),sd(W14.armA$cd4.cd8,na.rm=T),sd(W28.armA$cd4.cd8,na.rm=T),
              sd(W56.armA$cd4.cd8,na.rm=T),sd(W84.armA$cd4.cd8,na.rm=T),sd(W112.armA$cd4.cd8,na.rm=T),
              sd(W140.armA$cd4.cd8,na.rm=T),sd(W168.armA$cd4.cd8,na.rm=T))

# missing rate
armA.NA<-c(sum(is.na(W0.armA$cd4.cd8))/length(W0.armA$cd4.cd8),sum(is.na(W14.armA$cd4.cd8))/length(W14.armA$cd4.cd8),
           sum(is.na(W28.armA$cd4.cd8))/length(W28.armA$cd4.cd8),sum(is.na(W56.armA$cd4.cd8))/length(W56.armA$cd4.cd8),
           sum(is.na(W84.armA$cd4.cd8))/length(W84.armA$cd4.cd8),sum(is.na(W112.armA$cd4.cd8))/length(W112.armA$cd4.cd8),
           sum(is.na(W140.armA$cd4.cd8))/length(W140.armA$cd4.cd8),sum(is.na(W168.armA$cd4.cd8))/length(W168.armA$cd4.cd8))*100

################################## IDV-RTV=1 ###################################
### Group weeks
W0.armB<-armB[which(armB$day==0),]
W14.armB<-armB[which(armB$day<=14 & armB$day>0),]
W28.armB<-armB[which(armB$day<=28 & armB$day>14),]
W56.armB<-armB[which(armB$day<=56 & armB$day>28),]
W84.armB<-armB[which(armB$day<=84 & armB$day>56),]
W112.armB<-armB[which(armB$day<=112 & armB$day>84),]
W140.armB<-armB[which(armB$day<=140 & armB$day>112),]
W168.armB<-armB[which(armB$day>140),]

### log(RNA) ###
# mean
armB.rna.mean<-c(mean(W0.armB$logrna),mean(W14.armB$logrna),mean(W28.armB$logrna),
                 mean(W56.armB$logrna),mean(W84.armB$logrna),mean(W112.armB$logrna),
                 mean(W140.armB$logrna),mean(W168.armB$logrna))

# sd
armB.rna.sd<-c(sd(W0.armB$logrna),sd(W14.armB$logrna),sd(W28.armB$logrna),sd(W56.armB$logrna),
               sd(W84.armB$logrna),sd(W112.armB$logrna),sd(W140.armB$logrna),sd(W168.armB$logrna))

armB.CD.len<-c(sum(complete.cases(W0.armB$cd4.cd8)),sum(complete.cases(W14.armB$cd4.cd8)),
               sum(complete.cases(W28.armB$cd4.cd8)),sum(complete.cases(W56.armB$cd4.cd8)),
               sum(complete.cases(W84.armB$cd4.cd8)),sum(complete.cases(W112.armB$cd4.cd8)),
               sum(complete.cases(W140.armB$cd4.cd8)),sum(complete.cases(W168.armB$cd4.cd8)))

# censored rate
armB.cen<-c(sum(W0.armB$cen.lrna)/length(W0.armB$cen.lrna),sum(W14.armB$cen.lrna)/length(W14.armB$cen.lrna),
            sum(W28.armB$cen.lrna)/length(W28.armB$cen.lrna),sum(W56.armB$cen.lrna)/length(W56.armB$cen.lrna),
            sum(W84.armB$cen.lrna)/length(W84.armB$cen.lrna),sum(W112.armB$cen.lrna)/length(W112.armB$cen.lrna),
            sum(W140.armB$cen.lrna)/length(W140.armB$cen.lrna),sum(W168.armB$cen.lrna)/length(W168.armB$cen.lrna))*100

### CD4/CD8 ###
# mean
armB.CD.mean<-c(mean(W0.armB$cd4.cd8,na.rm=T),mean(W14.armB$cd4.cd8,na.rm=T),mean(W28.armB$cd4.cd8,na.rm=T),
                mean(W56.armB$cd4.cd8,na.rm=T),mean(W84.armB$cd4.cd8,na.rm=T),mean(W112.armB$cd4.cd8,na.rm=T),
                mean(W140.armB$cd4.cd8,na.rm=T),mean(W168.armB$cd4.cd8,na.rm=T))

# sd
armB.CD.sd<-c(sd(W0.armB$cd4.cd8,na.rm=T),sd(W14.armB$cd4.cd8,na.rm=T),sd(W28.armB$cd4.cd8,na.rm=T),
              sd(W56.armB$cd4.cd8,na.rm=T),sd(W84.armB$cd4.cd8,na.rm=T),sd(W112.armB$cd4.cd8,na.rm=T),
              sd(W140.armB$cd4.cd8,na.rm=T),sd(W168.armB$cd4.cd8,na.rm=T))

# missing rate
armB.NA<-c(sum(is.na(W0.armB$cd4.cd8))/length(W0.armB$cd4.cd8),sum(is.na(W14.armB$cd4.cd8))/length(W14.armB$cd4.cd8),
           sum(is.na(W28.armB$cd4.cd8))/length(W28.armB$cd4.cd8),sum(is.na(W56.armB$cd4.cd8))/length(W56.armB$cd4.cd8),
           sum(is.na(W84.armB$cd4.cd8))/length(W84.armB$cd4.cd8),sum(is.na(W112.armB$cd4.cd8))/length(W112.armB$cd4.cd8),
           sum(is.na(W140.armB$cd4.cd8))/length(W140.armB$cd4.cd8),sum(is.na(W168.armB$cd4.cd8))/length(W168.armB$cd4.cd8))*100

##################################### Table S2 #################################
TableS2 <- round(rbind(armA.rna.mean, armA.rna.sd, armA.cen,
                         armA.CD.mean, armA.CD.sd, armA.NA,
                         armB.rna.mean, armB.rna.sd, armB.cen,
                         armB.CD.mean, armB.CD.sd, armB.NA),4)

colnames(TableS2)<-c(0,14,28,56,84,112,140,168)
rownames(TableS2)<-c('A.RNA.Mean','A.RNA.SD','A.censor','A.CD4CD8.Mean','A.CD4CD8.SD','A.CD4CD8.na',
                     'B.RNA.Mean','B.RNA.SD','B.censor','B.CD4CD8.Mean','B.CD4CD8.SD','B.CD4CD8.na')

write.table(round(TableS2, 4), paste(PATH,'results/TableS2.csv', sep=""), sep = ",",na = "NA", row.names =T, col.names = T)
