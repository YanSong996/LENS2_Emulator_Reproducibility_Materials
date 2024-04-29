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













