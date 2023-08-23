rm(list = ls())
PATH = paste(getwd(),"/Data_and_Code_METR23/",sep="")

###### Simulation ######
# Re-generate simulation results for sample size N=25
N = 25
PATH1 = paste(PATH, 'results/N=25/', sep='')
source(paste(PATH, 'code/sim.r', sep=''))

# Re-generate simulation results for sample size N=50
N = 50
PATH1 = paste(PATH, 'results/N=50/', sep='')
source(paste(PATH, 'code/sim.r', sep=''))

# Re-generate simulation results for sample size N=75
N = 75
PATH1 = paste(PATH, 'results/N=75/', sep='')
source(paste(PATH, 'code/sim.r', sep=''))

# Re-generate simulation results for sample size N=100
N = 100
PATH1 = paste(PATH, 'results/N=100/', sep='')
source(paste(PATH, 'code/sim.r', sep=''))

# Re-produce Figure 1
source(paste(PATH, 'code/SIMfig1.r', sep=''))

# Re-produce Figure 2
source(paste(PATH, 'code/SIMfig2.r', sep=''))

# Re-produce Supplementary Figure S.1
source(paste(PATH, 'code/SIMfigS1.r', sep=''))

# Re-produce Supplementary Table S.1
source(paste(PATH, 'code/SIMtabS1.r', sep=''))

###### Application to A5055 Data ######
# Re-produce Table 1
source(paste(PATH, 'code/Tab1.r', sep=''))

# Re-produce Table 2
source(paste(PATH, 'code/Tab2.r', sep=''))

# Re-produce Figure 3a & 3b
source(paste(PATH, 'code/Fig3.r', sep=''))

# Re-produce Figure 4
source(paste(PATH, 'code/Fig4.r', sep=''))

# Re-perform model fitting
source(paste(PATH, 'code/A5055Fit.r', sep=''))

# Re-produce Supplementary Table S.2
source(paste(PATH, 'code/TabS2.r', sep=''))

# Re-produce Supplementary Figure S.2
source(paste(PATH, 'code/FigS2.r', sep=''))
