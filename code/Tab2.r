# Reproduct Table 2
load(paste(PATH, 'data/A5055FitResult.RData',sep=""))
EST = cbind(t(IBM2$pos.summ$theta.out)[,-5], t(IBCM1$pos.inf$theta.out)[,-5])
EST[c(29:33), ] = EST[c(29:33), ] * 1e4
Table2 = round(EST, 3)
rownames(Table2) = c(paste('beta', 1:12, sep=''), 
                     c('d11','d21','d22','d31','d32','d33','d41','d42','d43','d44',
                       'd51','d52','d53','d54','d55',
                       'd61','d62','d63','d64','d65','d66'),
                     c('sig11', 'sig21', 'sig22'), 'omega1', 'omega2')
PATH = paste(getwd(),"/Data_and_Code_METR23/",sep="")
write.table(Table2, paste(PATH,'results/Table2.csv', sep=""), sep = ",",na = "NA", row.names = T, col.names = T)
