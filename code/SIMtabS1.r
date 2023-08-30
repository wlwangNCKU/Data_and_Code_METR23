PATH1 = paste(PATH, 'results/N=25/', sep='')
PATH2 = paste(PATH, 'results/N=50/', sep='')
PATH3 = paste(PATH, 'results/N=75/', sep='')
PATH4 = paste(PATH, 'results/N=100/', sep='')

################################# CP & Length ##################################
CPLen.table = function(baycm, baym)
{
  bayM.CP = c(mean(baym[which(baym$para=='b1'),10],na.rm=T), mean(baym[which(baym$para=='b2'),10],na.rm=T),
              mean(baym[which(baym$para=='b3'),10],na.rm=T), mean(baym[which(baym$para=='b4'),10],na.rm=T),
              mean(baym[which(baym$para=='b5'),10],na.rm=T), mean(baym[which(baym$para=='d11'),10],na.rm=T),
              mean(baym[which(baym$para=='d21'),10],na.rm=T), mean(baym[which(baym$para=='d22'),10],na.rm=T),
              mean(baym[which(baym$para=='sigma11'),10],na.rm=T), mean(baym[which(baym$para=='sigma21'),10],na.rm=T),
              mean(baym[which(baym$para=='sigma22'),10],na.rm=T))

  bayCM.CP = c(mean(baycm[which(baycm$para=='b1'),10],na.rm=T), mean(baycm[which(baycm$para=='b2'),10],na.rm=T),
               mean(baycm[which(baycm$para=='b3'),10],na.rm=T), mean(baycm[which(baycm$para=='b4'),10],na.rm=T),
               mean(baycm[which(baycm$para=='b5'),10],na.rm=T), mean(baycm[which(baycm$para=='d11'),10],na.rm=T),
               mean(baycm[which(baycm$para=='d21'),10],na.rm=T), mean(baycm[which(baycm$para=='d22'),10],na.rm=T),
               mean(baycm[which(baycm$para=='sigma11'),10],na.rm=T), mean(baycm[which(baycm$para=='sigma21'),10],na.rm=T),
               mean(baycm[which(baycm$para=='sigma22'),10],na.rm=T))

  bayM.Len = c(mean(baym[which(baym$para=='b1'),9],na.rm=T), mean(baym[which(baym$para=='b2'),9],na.rm=T),
               mean(baym[which(baym$para=='b3'),9],na.rm=T), mean(baym[which(baym$para=='b4'),9],na.rm=T),
               mean(baym[which(baym$para=='b5'),9],na.rm=T), mean(baym[which(baym$para=='d11'),9],na.rm=T),
               mean(baym[which(baym$para=='d21'),9],na.rm=T), mean(baym[which(baym$para=='d22'),9],na.rm=T),
               mean(baym[which(baym$para=='sigma11'),9],na.rm=T), mean(baym[which(baym$para=='sigma21'),9],na.rm=T),
               mean(baym[which(baym$para=='sigma22'),9],na.rm=T))

  bayCM.Len = c(mean(baycm[which(baycm$para=='b1'),9],na.rm=T), mean(baycm[which(baycm$para=='b2'),9],na.rm=T),
                mean(baycm[which(baycm$para=='b3'),9],na.rm=T), mean(baycm[which(baycm$para=='b4'),9],na.rm=T),
                mean(baycm[which(baycm$para=='b5'),9],na.rm=T), mean(baycm[which(baycm$para=='d11'),9],na.rm=T),
                mean(baycm[which(baycm$para=='d21'),9],na.rm=T), mean(baycm[which(baycm$para=='d22'),9],na.rm=T),
                mean(baycm[which(baycm$para=='sigma11'),9],na.rm=T), mean(baycm[which(baycm$para=='sigma21'),9],na.rm=T),
                mean(baycm[which(baycm$para=='sigma22'),9],na.rm=T))
 
  bay = rbind(bayM.CP, bayCM.CP, bayM.Len, bayCM.Len)
  colnames(bay) = c("b1","b2","b3","b4","b5","d11","d21","d22","sigma11","sigma21","sigma22")
  rownames(bay) = c("CP.bayM", "CP.bayCM", "Len.bayM","Len.bayCM")
  return(bay)
}

##################################### N=25 #####################################
bayesian.m25<-read.table(paste(PATH1, "table(bayM).txt",sep=""))
colnames(bayesian.m25)<-c("loop","para","case","real data","Est","SE","LCI","UCI","Len","CP")
bayesian.cm25<-read.table(paste(PATH1, "table(bayCM).txt",sep=""))
colnames(bayesian.cm25)<-c("loop","para","case","real data","Est","SE","LCI","UCI","Len","CP")

bayesianm25.c1<-bayesian.m25[bayesian.m25$case=='1',]
bayesianm25.c2<-bayesian.m25[bayesian.m25$case=='2',]
bayesianm25.c3<-bayesian.m25[bayesian.m25$case=='3',]
bayesianm25.c4<-bayesian.m25[bayesian.m25$case=='4',]
bayesianm25.c5<-bayesian.m25[bayesian.m25$case=='5',]
bayesianm25.c6<-bayesian.m25[bayesian.m25$case=='6',]

bayesiancm25.c1<-bayesian.cm25[bayesian.cm25$case=='1',]
bayesiancm25.c2<-bayesian.cm25[bayesian.cm25$case=='2',]
bayesiancm25.c3<-bayesian.cm25[bayesian.cm25$case=='3',]
bayesiancm25.c4<-bayesian.cm25[bayesian.cm25$case=='4',]
bayesiancm25.c5<-bayesian.cm25[bayesian.cm25$case=='5',]
bayesiancm25.c6<-bayesian.cm25[bayesian.cm25$case=='6',]

case1.25<-CPLen.table(bayesiancm25.c1, bayesianm25.c1) 
case2.25<-CPLen.table(bayesiancm25.c2, bayesianm25.c2)
case3.25<-CPLen.table(bayesiancm25.c3, bayesianm25.c3) 
case4.25<-CPLen.table(bayesiancm25.c4, bayesianm25.c4)
case5.25<-CPLen.table(bayesiancm25.c5, bayesianm25.c5) 
case6.25<-CPLen.table(bayesiancm25.c6, bayesianm25.c6)

##################################### N=50 #####################################
bayesian.m50<-read.table(paste(PATH2, "table(bayM).txt",sep=""))
colnames(bayesian.m50)<-c("loop","para","case","real data","Est","SE","LCI","UCI","Len","CP")
bayesian.cm50<-read.table(paste(PATH2, "table(bayCM).txt",sep=""))
colnames(bayesian.cm50)<-c("loop","para","case","real data","Est","SE","LCI","UCI","Len","CP")

bayesianm50.c1<-bayesian.m50[bayesian.m50$case=='1',]
bayesianm50.c2<-bayesian.m50[bayesian.m50$case=='2',]
bayesianm50.c3<-bayesian.m50[bayesian.m50$case=='3',]
bayesianm50.c4<-bayesian.m50[bayesian.m50$case=='4',]
bayesianm50.c5<-bayesian.m50[bayesian.m50$case=='5',]
bayesianm50.c6<-bayesian.m50[bayesian.m50$case=='6',]

bayesiancm50.c1<-bayesian.cm50[bayesian.cm50$case=='1',]
bayesiancm50.c2<-bayesian.cm50[bayesian.cm50$case=='2',]
bayesiancm50.c3<-bayesian.cm50[bayesian.cm50$case=='3',]
bayesiancm50.c4<-bayesian.cm50[bayesian.cm50$case=='4',]
bayesiancm50.c5<-bayesian.cm50[bayesian.cm50$case=='5',]
bayesiancm50.c6<-bayesian.cm50[bayesian.cm50$case=='6',]

case1.50<-CPLen.table(bayesiancm50.c1, bayesianm50.c1) 
case2.50<-CPLen.table(bayesiancm50.c2, bayesianm50.c2)
case3.50<-CPLen.table(bayesiancm50.c3, bayesianm50.c3) 
case4.50<-CPLen.table(bayesiancm50.c4, bayesianm50.c4)
case5.50<-CPLen.table(bayesiancm50.c5, bayesianm50.c5) 
case6.50<-CPLen.table(bayesiancm50.c6, bayesianm50.c6)

##################################### N=75 #####################################
bayesian.m75<-read.table(paste(PATH3, "table(bayM).txt",sep=""))
colnames(bayesian.m75)<-c("loop","para","case","real data","Est","SE","LCI","UCI","Len","CP")
bayesian.cm75<-read.table(paste(PATH3, "table(bayCM).txt",sep=""))
colnames(bayesian.cm75)<-c("loop","para","case","real data","Est","SE","LCI","UCI","Len","CP")

bayesianm75.c1<-bayesian.m75[bayesian.m75$case=='1',]
bayesianm75.c2<-bayesian.m75[bayesian.m75$case=='2',]
bayesianm75.c3<-bayesian.m75[bayesian.m75$case=='3',]
bayesianm75.c4<-bayesian.m75[bayesian.m75$case=='4',]
bayesianm75.c5<-bayesian.m75[bayesian.m75$case=='5',]
bayesianm75.c6<-bayesian.m75[bayesian.m75$case=='6',]

bayesiancm75.c1<-bayesian.cm75[bayesian.cm75$case=='1',]
bayesiancm75.c2<-bayesian.cm75[bayesian.cm75$case=='2',]
bayesiancm75.c3<-bayesian.cm75[bayesian.cm75$case=='3',]
bayesiancm75.c4<-bayesian.cm75[bayesian.cm75$case=='4',]
bayesiancm75.c5<-bayesian.cm75[bayesian.cm75$case=='5',]
bayesiancm75.c6<-bayesian.cm75[bayesian.cm75$case=='6',]

case1.75<-CPLen.table(bayesiancm75.c1, bayesianm75.c1) 
case2.75<-CPLen.table(bayesiancm75.c2, bayesianm75.c2)
case3.75<-CPLen.table(bayesiancm75.c3, bayesianm75.c3) 
case4.75<-CPLen.table(bayesiancm75.c4, bayesianm75.c4)
case5.75<-CPLen.table(bayesiancm75.c5, bayesianm75.c5) 
case6.75<-CPLen.table(bayesiancm75.c6, bayesianm75.c6)

##################################### N=100 #####################################
bayesian.m100<-read.table(paste(PATH4, "table(bayM).txt",sep=""))
colnames(bayesian.m100)<-c("loop","para","case","real data","Est","SE","LCI","UCI","Len","CP")
bayesian.cm100<-read.table(paste(PATH4, "table(bayCM).txt",sep=""))
colnames(bayesian.cm100)<-c("loop","para","case","real data","Est","SE","LCI","UCI","Len","CP")

bayesianm100.c1<-bayesian.m100[bayesian.m100$case=='1',]
bayesianm100.c2<-bayesian.m100[bayesian.m100$case=='2',]
bayesianm100.c3<-bayesian.m100[bayesian.m100$case=='3',]
bayesianm100.c4<-bayesian.m100[bayesian.m100$case=='4',]
bayesianm100.c5<-bayesian.m100[bayesian.m100$case=='5',]
bayesianm100.c6<-bayesian.m100[bayesian.m100$case=='6',]

bayesiancm100.c1<-bayesian.cm100[bayesian.cm100$case=='1',]
bayesiancm100.c2<-bayesian.cm100[bayesian.cm100$case=='2',]
bayesiancm100.c3<-bayesian.cm100[bayesian.cm100$case=='3',]
bayesiancm100.c4<-bayesian.cm100[bayesian.cm100$case=='4',]
bayesiancm100.c5<-bayesian.cm100[bayesian.cm100$case=='5',]
bayesiancm100.c6<-bayesian.cm100[bayesian.cm100$case=='6',]

case1.100<-CPLen.table(bayesiancm100.c1, bayesianm100.c1) 
case2.100<-CPLen.table(bayesiancm100.c2, bayesianm100.c2)
case3.100<-CPLen.table(bayesiancm100.c3, bayesianm100.c3) 
case4.100<-CPLen.table(bayesiancm100.c4, bayesianm100.c4)
case5.100<-CPLen.table(bayesiancm100.c5, bayesianm100.c5) 
case6.100<-CPLen.table(bayesiancm100.c6, bayesianm100.c6)

######################## Supplementary Table S.1 ###############################
Tab.CPLen = data.frame(rep(c('MAR', 'MCAR'), each=48), 
                  c(rep('10%', 16), rep('30%', 16), rep('50%', 16)), 
                  rep(rep(c(25, 50, 75, 100), each=4), 6),
                  rep(rep(c('CP','Length'), each=2), 24),
                  rep(c('M', 'CM'), 48),  
                  round(rbind(case1.25, case1.50, case1.75, case1.100, 
                        case3.25, case3.50, case3.75, case3.100,
                        case5.25, case5.50, case5.75, case5.100, 
                        case2.25, case2.50, case2.75, case2.100,
                        case4.25, case4.50, case4.75, case4.100,
                        case6.25, case6.50, case6.75, case6.100), 4))

colnames(Tab.CPLen) = c('Missingness', 'CensorRate', 'SampleSize', 'Criterion', 'Model',
       "b1","b2","b3","b4","b5","d11","d21","d22","sigma11","sigma21","sigma22")

write.table(Tab.CPLen, paste(PATH, 'results/SIMTabS1.csv',sep=""), 
            sep = ",", na = "NA", row.names = F, col.names = T)





