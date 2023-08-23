# Reproduct Table 1
load(paste(PATH, 'data/A5055FitResult.RData',sep=""))
Table1 = round(cbind(rbind(c(M0$model.inf$num.par, M1$model.inf$num.par, M2$model.inf$num.par),
                         c(IBM0$model.inf$EAIC, IBM1$model.inf$EAIC, IBM2$model.inf$EAIC),
                         c(IBM0$model.inf$EBIC, IBM1$model.inf$EBIC, IBM2$model.inf$EBIC),
                         c(IBM0$model.inf$DIC, IBM1$model.inf$DIC, IBM2$model.inf$DIC)),
                   rbind(c(CM0$model.inf$m, CM1$model.inf$m, CM2$model.inf$m),
                         c(IBCM0$model.inf$EAIC, IBCM1$model.inf$EAIC, IBCM2$model.inf$EAIC),
                         c(IBCM0$model.inf$EBIC, IBCM1$model.inf$EBIC, IBCM2$model.inf$EBIC),
                         c(IBCM0$model.inf$DIC, IBCM1$model.inf$DIC, IBCM2$model.inf$DIC))), 3)

rownames(Table1) = c('m', 'EAIC', 'EBIC', 'DIC')
colnames(Table1) = c('M-UNC', 'M-AR1', 'M-DEC', 'CM-UNC', 'CM-AR1', 'CM-DEC')
PATH = paste(getwd(),"/Data_and_Code_METR23/",sep="")
write.table(Table1, paste(PATH,'results/Table1.csv', sep=""), sep = ",",na = "NA", row.names =T, col.names = T)
