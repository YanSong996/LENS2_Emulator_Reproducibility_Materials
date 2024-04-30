# Reproducibility Introduction for the Supplementary Materials
This file documents the artifacts associated with the Supplementary Materials and describes how to reproduce the findings. To reproduce all figures and results in the Supplementary Materials, please download and extract the provided "LENS2_Emulator_Reproducibility_Materials.zip", and then set your working directory to the folder "LENS2_Emulator_Reproducibility_Materials".

## Step 1. Download and process the data in sub-repository "LENS2_Data"
Please follow the README.md file in the sub-repository "LENS2_Data" to download the raw data, process and store them using file “Data_treatment.R” in the "LENS2_Data". Processing montly and annual data together for each ensemble takes about 63 seconds. Processing daily data for each ensemble takes about 32.5 seconds. 

## Step 2. Reproduce each section and figure sequentially
### Figure S1 in Section S2
Please refer to the reproduction of Figure 1 in Section 2.

### Figure S2 in Section S3.2.3
Please refer to the reproduction of Figure 2 in Section 3.2.1.

### Figure S3 in Section S3.2.4
Figure S3 illustrates the approximation performance of LatticeKrig and SHT, where Figures S3(d) and S3(b) are Figures 2(e) and 2(f), respectively. More details about choosing tuning parameters can be found in "Supplement/FigureS3/LatticeKrig_Annual.R" All intermediate outputs can be found in sub-repository "Supplement/FigureS3/Outputs". The total computational time for ploting Figure S3 (excluding Figures S3(d) and S3(b)) is about 2569.2 seconds.

1. Load 7 ensembles of annual data and necessary information (about 5.2 seconds)
2. Calculate one set of stochastic component $Z_9^{(1)}(L_i,l_j)$ (about 2.4 seconds)
3. Use LatticeKrig with nlevel=5 to approximate $Z_9^{(1)}(L_i,l_j)$ and plot Figure S3(a) (about 180.1 seconds)
4. Use SHT with Q=58 to approximate $Z_9^{(1)}(L_i,l_j)$ and plot Figure S3(c) (about 4.2 seconds)
5. Compare LatticeKrig and SHT using multiple stochastic components and plot Figures S3(e) and S3(f) (about 2377.3 seconds using 4 cores)

### Figure S4 in Section S3.4
Figure S4 validates the assumption of diagonal matrix $\Phi_p$. Here we take the annual case as an example to illustrate the reproduction. Assume that we have got the real-valued SHT coefficients "TDat.rsd.SHT" and their temporal dependence "TPhi.hat" by following lines 1-192 in code "Annual/Annual_SG.R". The total computational time is about 396.1 seconds.

1. Calculate the residuals of autoregressive model (about 0.2 seconds)
2. Calculate p-values of the first temporal lag of the cross-correlation and plot Figure S4(a) (about 395.9 seconds)

### Figure S5 in Section S4.1.1
Please refer to the reproduction of Figure 3 in Section 4.1.

### Figure S6 in Section S4.1.2
Figure S6 illustrates the inference results of annual data obtained by HCBG. All the intermediate results can be found in sub-repository "Supplement/Huang/Annual/Outputs". Here we use outputs from https://github.com/hhuang90/stochastic_emulator 
as initial values. The total computational time for ploting Figure S6 is about 113278.9 seconds without using the parallel. Beyond Figure S6, other procedures for generating annual emulations using HCBG and their computational time are provided in "Supplement/Huang/Annual/Huang_Annual.R".

1. Load 7 ensembles of annual data and necessary information for HCBG (about 4.6 seconds)
2. Model the mean trend and temporal dependence at each grid point and plot Figures S6(a) (about 2.0 seconds for each grid point $(L_i,l_j)$. We provide the intermediate result "Res_hatrhophi.csv".)
3. Evaluate mean trend and sd at each grid point and plot Figure S6(b) (about 723.9 seconds. We provide the intermediate result "Sig_Annual.csv".)
4. Model the longitude dependence at each latitude and plot Figure S6(c) (about 10.2 seconds for each latitude $L_i$, $i=1,\ldots,192$. We provide the intermediate result "Res_longLO.csv".)

### Figures S7 and S8 in Section S4.1.3
Please refer to the reproduction of Figure 5 in Section 4.1.

### Figure S9 in Section S4.2
Figure S9 illustrates inference results for deterministic components of the monthly data. All intermediate outputs can be found in sub-repository "Monthly/Outputs". The total computational time is about 774210.1 seconds without using the parallel. Note that Figures S9-S13 are all for the monthly data, so please do not clear the environment after reproducing Figure S9 if you also want to reproduce Figures S10-S13.           

1. Load 7 ensembles of monthly data and necessary information (about 66.1 secconds)
2. Model the deterministic components, evaluate mean trend $m_t(L_i,l_j)$ and sd $\sigma(L_i,l_j)$ for each grid point $(L_i,l_j)$, and plot Figure S9 (about 14.0 seconds for each grid point. We provide intermediate results "Res_Hatrho.csv", "Res_BetaAB.mat", and "Monthly_Sig.csv".)

### Figure S10 in Section S4.2
Figure S10 illustrates inference results for stochastic components of the monthly data. Assume that we keep all the intermediate results of Figure S9. The total computational time is about 237107.2 seconds without using the parallel. All intermediate outputs can be found in sub-repository "Monthly/Outputs".

1. Calculate stochastic components $Z_t^{(r)}(L_i,l_j)$  by detrending and rescaling (about 19.8 seconds)
2. Do SHT with $Q=144$ for the stochastic component at each ensemble $r$ and time point $t$ (about 4.2 seconds for each ensemble r and each time point t. Without using the parallel, it will take about $4.2\times 7\times (86\times 12)=30340.8 seconds.)
3. Calculate BIC values under different $Q$ values and plot Figure S10(a) (The major computation time is for calculating the inverse SHT. But this step takes time because it calculates the inverse SHT for each ensemble $r$, each time point $t$, and all candidates of $Q$. Take Q=90 as an example, which maximizes the computational time of inversing SHT, it takes about 0.9 seconds for each ensemble $t$ and time point $t$. If we run this step without doing parallel, we would take at most $0.9\times 7\times (86\times 12)\times (8+11+11)=195048 seconds. We provide intermediate results "BIC_land.csv", "BIC_ocean.csv", "BICd_land.csv", and "BICd_ocean.csv".)
4. Calculate v^2(L_i,l_j) under Ql=36 and Qo=70 and plot Figure S10(b) (The major computational time is to calculate the inverse SHT with $Q=36$ and $70$ for each $r$ and $t$. It will take about $0.8\times 7\times 86\times 12=5779.2$ seconds without the parallel. We provide the intermediate result "v2hat.csv".)
5. Do the real-valued transformation to SHT coefficients so that they are real values (about 2.8 seconds)
6. Test the Gaussianity of the real-valued coefficients (about 6.7 seconds)
7. Model the temporal dependence structure using a Tukey g-and-h autoregressive model with order $P=1$ and plot Figure S10(c) (about 5909.9 seconds without using the parallel. We provide intermediate results "bicp_noTukey.csv", "bicp_Tukey.csv", "Phihat_noTukey.csv", "Phihat_Tukey.csv", "rtc.csv", "Tukeyres.csv".)

### Figures S11 and S12 in Section S4.2
Figures S11 and S12 illustrates the performance of the generated monthly emulations. Assume that we keep all the intermediate results of Figures S9 and S10. The total computational time is about 4426.7 seconds with 4 cores. All intermediate outputs are in sub-repository "Monthly/Outputs".

1. Model the spatial dependence by evaluating the covariance matrix of (read-valued and Gaussianized) SHT coefficients (about 48.3 seconds)
2. Calculate the covaraince matrix $\Check{**U**}$ (about 0.3 seconds. We provide the intermediate result "U.mat")
3. Generate $R'=7$ ensembles of monthly emulations using 4 cores (about 2873.5 seconds)
4. Calculate $I_{uq}$ values using 4 cores and plot Figures S11(a), S11(b), and S11(e) (about 1241.9 seconds. We provide intermediate results "Iuq_Tukey.csv" and "Iuq_noTukey.csv".)
5. Calculate $WD_{S}$ values using 4 cores and plot Figures S11(c), S11(d), and S11(f) (about 262.7 seconds. We provide intermediate results "WD_time_Tukey.csv" and "WD_time_noTukey.csv".)
6. Plot Figure S12 (about 0 seconds)

### Figure S13 Section S4.2
Figure S13 illustrates the performance of the generated monthly emulations by aggregating them to be annual emulations. Assume that we keep all the intermediate results of Figures S9-S12, especially the generated monthly emulations. The total computational time is about 745.4 seconds. All intermediate outputs are in sub-repository "Monthly/Aggregate".

1. Aggregate monthly emulations to annual emulations (about 429.6 seconds)
2. Load the annual simulations (about 5.2 seconds)
3. Calculate $I_{uq}$ values and plot Figures S13(a) and S13(b) (about 139.2 seconds. We provide the intermediate result "Iuq_Monthly_to_Annual.csv".)
4. Calculate $WD_S$ values and plot Figures S13(b) and S13(c) (about 42.3 seconds. We provide the intermediate result "WD_Monthly_to_Annual.csv".)
5. Calculate $I_{fit}$ values and plot Figure S13(e) (about 129.1 seconds. We provide the intermediate result "Ifit_Monthly_to_Annual.csv".)
6. Plot Figure S13(f) using "BIC_KM.csv" in "Monthly/Outputs" (about 0 seconds)




