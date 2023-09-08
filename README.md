# Data_and_Code_METR23
Supplement: "Bayesian multivariate nonlinear mixed models for censored longitudinal trajectories with non-monotone missing values"

README

#######################################################################################

Source code and data for the manuscript 
"Bayesian multivariate nonlinear mixed models for censored longitudinal trajectories with non-monotone missing values",
by Wan-Lun Wang, Luis M. Castro, and Tsung-I Lin*

#######################################################################################

# Configurations
The code was written/evaluated in R with the following software versions:
R version 4.1.1 (2021-08-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default

locale:
[1] LC_COLLATE=Chinese (Traditional)_Taiwan.950  LC_CTYPE=Chinese (Traditional)_Taiwan.950    LC_MONETARY=Chinese (Traditional)_Taiwan.950
[4] LC_NUMERIC=C                                 LC_TIME=Chinese (Traditional)_Taiwan.950    

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggplot2_3.4.3   lemon_0.4.6     tmvtnorm_1.4-10 gmm_1.6-6       sandwich_3.0-1  Matrix_1.3-4    MCMCpack_1.6-0  MASS_7.3-54     coda_0.19-4     mvtnorm_1.1-2  
[11] nlme_3.1-152   

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.10        plyr_1.8.8         compiler_4.1.1     pillar_1.9.0       tools_4.1.1        lifecycle_1.0.3    tibble_3.2.1       gtable_0.3.3      
 [9] lattice_0.20-44    pkgconfig_2.0.3    rlang_1.1.0        cli_3.6.1          SparseM_1.81       xfun_0.39          gridExtra_2.3      withr_2.5.0       
[17] knitr_1.36         dplyr_1.1.2        generics_0.1.2     vctrs_0.6.1        MatrixModels_0.5-0 tidyselect_1.2.0   grid_4.1.1         glue_1.6.2        
[25] R6_2.5.1           fansi_0.5.0        farver_2.1.1       conquer_1.0.2      magrittr_2.0.1     scales_1.2.1       matrixStats_0.61.0 mcmc_0.9-7        
[33] colorspace_2.1-0   labeling_0.4.2     quantreg_5.86      utf8_1.2.2         munsell_0.5.0      zoo_1.8-9         

# Descriptions of the codes
Please download the folder file "Data_and_Code_METR23" to the "current working directory" of the R package.
The getwd() function shall determine an absolute pathname of the "current working directory".

Before running the codes 'sim.r' and 'A5055Fit.r', one needs to install the following R packages:

    install.packages("mvtnorm")  Version: 1.1-3
    install.packages("tmvtnorm") Version: 1.5
    install.packages("MCMCpack") Version：1.6-3
    install.packages("nlme")     Version: 3.1-162
    install.packages("lemon")    Version：0.4.6
    install.packages("ggplot2")  Version：3.4.3

R codes for the implementation of our methodology are provided.

## Subfolder: ./function ##
./function
	contains the program (function) of
 
 	(1) 'mnlmm.na.fn.r' for carrying out EM-based maximum likelihood (ML) estimation of multivariate nonlinear mixed models (MNLMM) with missing data;
  	(2) 'mnlmmcm.fn.r' for carrying out EM-based ML estimation of multivariate nonlinear mixed models with censored and missing responses (MNLMM-CM);
   	(3) 'Bay.mnlmm.na.fn.r' for performing Bayesian estimation of MNLMM with missing data; and
    (4) 'Bay.mnlmm.na.fn.r' for performing Bayesian estimation of MNLMM-CM for simulation studies.
    (5) 'mnlmm.na.fn2.r' for carrying out EM-based ML estimation of MNLMM with missing data;
    (6) 'mnlmmcm.fn2.r' for carrying out EM-based ML estimation of MNLMM-CM;
    (7) 'Bay.mnlmm.na.fn2.r' for performing Bayesian estimation of MNLMM with missing data; and
    (8) 'Bay.mnlmm.na.fn2.r' for performing Bayesian estimation of MNLMM-CM for the A5055 data example.

## Subfolder: ./code ##
./code
       contains 
       	
	(1) 'sim.r' main script for re-generating part of intermediate results for simualtion (note: The cases of sample sizes 'N=25', 'N=50', 'N=75' and 'N=100' should be done separately.);
	(2) 'SIMfig1.r' main script for reproducting Fig. 1; 
	(3) 'SIMfig2.r' main script for reproducting Fig. 2;
	(4) 'SIMfigS1.r' main script for reproducting Fig. S.1;
	(5) 'SIMtabS1.r' main script for reproducting Table S.1;
 	(6) 'Tab1.r' main script for reproducting Table 1 (load 'A5055FitResult.RData' directly, and then run 'Tab1.r'); 
	(7) 'Tab2.r' main script for reproducting Table 2 (load 'A5055FitResult.RData' directly, and then run 'Tab2.r');
	(8) 'TabS2.r' main script for reproducting Table S.2;
	(9) 'Fig3.r' main script for reproducting Fig. 3;
	(10) 'Fig4.r' main script for reproducting Fig. 4 (load 'A5055FitResult.RData' directly, and then run 'Fig4.r');
	(11) 'FigS2.r' main script for reproducting Fig. S.2;
	(12) 'A5055Fit.r' main script for performing ML and Bayesian model fitting to the A5055 dataset.

### Note for Section 4 - Simulation

(1) R code 'sim.r' generates the intermediate results of Figs. 1-2 in the manuscript and Fig. S.1 and Table S.1 in the supplementary materials.

(2) To conduct simulation studies, please source the 'sim.r' script in subfolder "./code/", and then run the 'SIMfig1.r', 'SIMfig2.r', 'SIMfigS1' and 'SIMtabS1' scripts in subfolder "./code/".

(3) Because the code takes a huge amount of time to run, we record these intermediate results so that one can use the R codes 'SIMfig1.r', 'SIMfig2.r', 'SIMfigS1' and 'SIMtabS1' to obtain the final results based on files stored in "./results/N=25", "./results/N=50", "./results/N=75" and '"./results/N=100"' subfolders.

### Note for Section 5 - Application to A5055 data

(1) Because the 'A5055Fit.r' code takes a huge amount of time to run the MCMC sampling procedure for Bayesian model fitting, we record these intermediate results in 'A5055FitResults.RData' so that 
    one can use the R codes 'Tab1.r', 'Tab2.r' and 'Fig4' to obtain the final results immediately.

(2) To reproduce the results presented in Tables 1-2 and Figure 4, just load 'A5055FitResult.RData' file in the "./data/" 
    and then run the script 'Tab1.r', 'Tab2.r', and 'Fig4.r' in the subfolder "./code/". 

## Subfolder: ./data ##
./data
      contains
      
      (1) 'A5055data.txt': the dataset for the A5055 HIV-AIDS study;
      (2) 'A5055FitResult.RData': the fitting results for A5055 dataset.

## Subfolder: ./results ##
./results
      contains 
        
	(1) 'SIMfig1.eps': (Fig. 1) boxplots of the mean squared errors for the posterior estimates of entire model parameters obtained by fitting the MNLMM-CM and MNLMM;
	(2) 'SIMfig2.eps': (Fig. 2) split violin plots of the mean squared errors for the posterior mean of fitted responses;
 	(3) 'SIMfigS1.eps': (Supplementary Fig. S.1) mean squared errors for the posterior estimates of parameters obtained by fitting the MNLMM-CM and MNLMM;
  	(4) 'SIMtabS1.csv': (Supplementary Table S.1) coverage probabilities and average lengths of the 95% posterior credible intervals for model parameters;

	(5) 'Tab1.csv': (Table 1) Bayesian model selection criteria under the six competing models for the A5055 data;
 	(6) 'Tab2.csv': (Table 2) posterior estimates under the fitted MNLMM with DEC errors and MNLMM-CM with AR(1) errors for the A5055 data;
  	(7) 'TabS2.csv': (Supplementary Table S.2) calculating averages (Mean) and standard deviations (SD) of log10(RNA) and CD4/CD8 ratios at scheduled days after drug-regimen treatments;
   	(8) 'Fig3a.eps' & 'Fig3b.eps': (Fig. 3) showing censoring and missingness patterns for the A5055 data;
    (9) 'Fig4.eps': (Fig. 4) showing observations, fitted values and 95% posterior bands of responses, and imputed values and 95% posterior credible intervals for missing CD4/CD8 ratios;
    (10) 'FigS2.eps': (Supplementary Fig. S.2) drawing trajectory plot for the log10(RNA) and CD4/CD8 ratios of 44 patients in the A5055 dataset;

./results/N=25; ./results/N=50; ./results/N=75; ./results/N=100
      contain intermediately numerical results from the simulation studies under four different sample sizes.
      
### Additional Remark 
One can directly run each "source(.)" described in 'master.r' file in the seperate R session to obtain the results.
