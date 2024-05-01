# Reproducibility Instruction
This file documents the artifacts associated with the article (i.e., the data and code supporting the computational findings) and describes how to reproduce the findings. To reproduce all figures and results in the main manuscript, please download and extract the provided "LENS2_Emulator_Reproducibility_Materials.zip", and then set your working directory to the folder "LENS2_Emulator_Reproducibility_Materials". 


## Step 1. Download and process the data in sub-repository "LENS2_Data"
Please follow the README.md file in the sub-repository "LENS2_Data" to download the raw data, process and store them using file “Data_treatment.R” in the "LENS2_Data". Processing montly and annual data together for each ensemble takes about 63 seconds. Processing daily data for each ensemble takes about 32.5 seconds.

## Step 2. Reproduce each section and figure sequentially
### Figure 1 in Section 2 (and Figure S1 in Section S2)
Figures 1 and S1 are used to illustrate the characteristics of surface temperature simulations in different scales. All the intermediate outputs can be found in sub-repository "Figure1". The total computational time for ploting Figures 1 and S1 is about 491.2 seconds.

1. Load 7 ensembles of annual data and necessary information (about 5.2 seconds)
2. Calculate annual ensemble mean to plot Figure 1(a) (about 1.9 seconds)
3. Calculate annual ensemble sd to plot Figure 1(b) (about 3.8 seconds)
4. Plot Figure 1(c) (about 0.6 seconds)
5. Save 7 annual time series at two grid points for ploting Figure 1(e) (about 0.4 seconds)
6. Calculate skewness and kurtosis for residuals of annaul data to plot Figures S1(a) and S1(d) (about 13.9 seconds)
7. Load 7 ensembles of monthly data and necessary information (about 71.6 seconds)
8. Plot Figure 1(d) (about 0.4 seconds)
9. Save 7 monthly time series at two grid points for ploting Figure 1(f) (about 0.5 seconds)
10. Calculate skewness and kurtosis for residuals of monthly data to plot Figures S1(b) and S1(e) (about 100.8 seconds)
11. Load 7 ensembles of daily data and necessary information (about 104.1 seconds)
12. Save 7 daily time series at two grid points for ploting Figure 1(g) (about 0.4 seconds)
13. Calculate skewness and kurtosis for residuals of daily data to plot Figures S1(c) and S1(f) (about 180.7 seconds)
14. Plot Figures 1(e)-1(g) (about 3.2 seconds)
15. Plot Figures S1(a)-S1(c) (about 1.8 seconds) and Figures S1(d)-S1(f) (about 1.7 seconds)

### Figure 2 in Section 3.2.1 (and Figure S2 in Section S3.2.3)
Figures 2 and S2 are used to illustrate the performance of spherical harmonics transformation (SHT). The total computational time for ploting Figures 2 and S2 is about 841.7 seconds. The intermediate outputs for Figure 2(f), which takes major time to reproduce, are saved in sub-repository "Figure2".

1. Load 7 ensembles of annual data and necessary information (about 5.2 seconds)
2. Calculate one set of stochastic component $Z_9^{(1)}(L_i,l_j)$ and plot Figure 2(a) (about 2.4 seconds)
3. Do SHT for $Z_9^{(1)}(L_i,l_j)$ and plot Figure 2(b) (about 4.2 seconds for SHT)
4. Do inverse SHT with Q=36 and plot Figure 2(c) (about 0.3 seconds for inverse SHT)
5. Do inverse SHT with Q=72 and plot Figure 2(d) (about 0.4 seconds for inverse SHT)
6. Do inverse SHT with Q=116 and plot Figure 2(e) (about 1.3 seconds for inverse SHT)
7. Use LatticeKrig to approximate $Z_9^{(1)}(L_i,l_j)$ and plot Figure 2(f) (about 814.5 seconds for LatticeKrig. We provide the intermediate results "LatticeKrig_level6_res.csv" in sub-repository "Figure2".)
8. Calculate one set of stochastic component $Z_9^{(3)}(L_i,l_j)$, do SHT and inverse SHT, and plot Figure S2(a) and S2(c) (about 6.0 seconds)
9. Calculate one set of stochastic component $Z_{69}^{(1)}(L_i,l_j)$, do SHT and inverse SHT, and plot Figure S2(b) and S2(d) (about 7.4 seconds)

### Figure 3 in Section 4.1 (and Figure S5 in Section S4.1.1)
Figure 3 illustrates some inference results of the annual data. Assume that we would use the intermediate results in item 2, then the total computational time without using parallel is about 63765.6 seconds. The detailed computational time is provided in each item below. We also give all the intermediate results in the sub-repository "Annual/Outputs". Note that Figures 3-5 are all for annual data, so please do not clear the environment after reproducing Figure 3 if you want to reproduce Figures 4 and 5 as well.

1. Load 10 ensembles of annual data and necessary information (about 11.0 seconds)
2. Calculate $I_{fit}$ values under $R=2,\ldots,10$ and plot Figure 3(a) and S5(a) (Here we give the maximum computational time, about 10917.9 seconds, which is obtained by calculating $I_{fit}$ values under $R=10$ and using 10 cores in parallel. Then, the whole step will take no more than 98261.2 seconds. We provide the intermediate result "IfitwithRs.csv" so that readers could use it to plot Figure 3(a) directly if necessary.)
3. Evaluate the deterministic components $m_t(L_i,l_j)$ and $\sigma(L_i,l_j)$ at each grid point using $R=7$ ensembles and plot Figure 3(b) and S5(b) (about 0.6 seconds for each grid point $(L_i,l_j)$, $i=1,\ldots,I=192$ and $j=1,\ldots,J=288$. Without parallel, this step takes 33177.6 seconds. We provide the intermediate results "Res_hatrho.csv", "Res_hatBeta.csv" and "Sig_Annual.csv".)
4. Calculate stochastic components $Z_t^{(r)}(L_i,l_j)$ by detrending $m_t(L_i,l_j)$ and rescaling $\sigma(L_i,l_j)$ (about 1.4 seconds)
5. Do SHT with $Q=144$ for the stochastic component $Z_t^{(r)}(L_i,l_j)$ at each ensemble $r$ and time point $t$ (about 4.2 seconds for each ensemble $r$ and each time point $t$, $r=1,\ldots,7$ and $t=1,\ldots,86$. Without parallel, this step takes about 2528.4 seconds.)
6. Calculate BIC values under different $Q$ values and plot Figure 3(c) (The major computation time is to calculate the inverse SHT with different values of $Q$. But this step takes time because it calculates the inverse SHT for each ensemble $r$, each time point $t$, and all candidates of $Q$. Take $Q=100$ as an example, which maximizes the computational time of inversing SHT in this step, it takes about 0.9 seconds for each $r$ and $t$. Without parallel, this step takes at most $0.9\times 7\times 86\times (9+21+21)=27631.8$ seconds. We provide the intermediate results "BIC_land.csv", "BIC_ocean.csv", "BICd_land.csv", and "BICd_ocean.csv".)
7. Calculate $v^2(L_i,l_j)$ under $Q_l=35$ and $Q_o=69$ and plot Figure 3(d) (The major computational time in this step is to calculate the inverse SHT with $Q=35$ and $69$ for each ensemble $r$, time point $t$. Without paralle, this step takes about $(0.27+0.42)\times 7\times 86=415.38$ seconds. We provide the intermediate outputs "v2hat.csv".)

### Figure 4 in Section 4.1
Figure 4 illustrates the temporal and spatial dependence structures of the annual data in spectral domain. Assume that we keep all intermediate results of reproducing Figure 3. The total computational time for ploting Figure 4 is about 439.4 seconds. All the intermediate outputs can be found in sub-repository "Annual/Outputs".

1. Do the real-valued transformation to SHT coefficients so that they are real-valued (about 1.4 seconds)
2. Model the temporal dependence structure using an autoregressive model with order $P=1$ and plot Figure 4(a) (about 0.357 seconds. We provide the intermediate results "Phihat.csv".)
3. Model the spatial dependence by evaluating the covariance matrix of (read-valued) SHT coefficients and plot Figures 4(b) and 4(c) (about 127.1 seconds. We provide the intermediate results "BIC_axial.csv".)
4. Assess spatial models by selecting auto-covariance at a selected latitide and plot Figures 4(d) (about 155.3 seconds)
5. As the same as the above item but use another latitude (about 155.3 seconds)

### Figure 5 in Section 4.1 (and Figures S7 and S8 in Section S4.1.3)
Figures 5, S7, and S8 illustrate the performance of the generated annual emulations. Assume that we keep all intermediate results of reproducing Figures 3 and 4. The total computational time for ploting Figures 5, S7, and S8 is about 369.0 seconds. All the intermediate outputs can be found in sub-repository "Annual/Outputs".

1. Caculate the covaraince matrix $\tilde{**U**}$ (about 0.4 seconds. We provide the intermediate result "U.mat".)
2. Do Cholescky decompostion on $\tilde{**U**}$ (about 0.6 seconds)
3. Generate $R'=7$ ensembles of annual emulations using 4 cores (about 207.8 seconds)
4. Calculate $I_{uq}$ values using 4 cores and plot Figures 5(a) and 5(b) (about 119.1 seconds. We provide the intermediate result "Iuq_axialnon.csv".)
5. Calculate $WD_S$ values using 4 cores and plot Figures 5(c) and 5(d) (about 31.1 seconds. We provide the intermediate result "WD_time.csv".)
6. Compare ensemble means and sds of annual emulations with those of simulations and plot Figures S7(a) and S7(c) (or Figures S7(b) and S7(d)) (about 2.4 seconds)
7. Compare time series of annual emulations with those of simulations and plot Figures S7(e) and S7(f) (0 seconds)
8. Compare periodograms of emulations with those of simulations and plot Figure S8(a) (or S8(b)) (about 1.4 seconds)

### Figures 6 and 7 in Section 4.2 (and Figures S14 and S15 in Section S4.3)
Figures 6 and 7 illustrate the performance of the generated daily emulations. These two figures take about 6199.6 seconds using the parallel with 4 cores. The inference process, which should be done before generating emulations, is illustrated in Figure S14. This figure takes about 2784269 seconds without using the parallel. Figure S15 just shows daily time series. All these four plots are for daily data. The order should be Figure S14, 6, 7 and S15. Therefore, we give the entire procedure here together. All intermediate outputs are in sub-repository "Daily/Outputs".  

1. Load 7 ensembles of daily data and necessary information (about 116.5 seconds)
2. Model the deterministic components, evaluate $m_t(L_i,l_j)$ and $\sigma(L_i,l_j)$ for each grid point $(L_i,l_j)$, and plot Figures S14(a) and S14(b) (about 39.6 seconds for each $(L_i,l_j)$. Without the parallel, this steop takes about $39.6\times 192\times 288=2191878$ seconds. We provide intermediate results "BICforM.csv", "Res_Hatrho.csv", "BetaAB.mat", and "Daily_Sig.csv".)
3. Calculate stochastic components $Z_t^{(r)}(L_i,l_j)$ by detrending and rescaling (about 32.2 seconds)
4. Do SHT with $Q=144$ for the stochastic component at each ensemble $r$ and time point $t$ (about 4.2 seconds for each $r$ and each $t$. Without the parallel, it takes 4.2*7*(5*365)=53655 seconds.)
5. Calculate BIC values under different Q values and plot Figure S14(c) (The major computation time is for calculating the inverse SHT. But this step takes time because it calculates the inverse SHT for $r$, each $t$, and all candidates of $Q$. Take $Q=90$ as an example, which maximizes the computational time of inversing SHT, it takes about 0.9 seconds for each $r$ and $t$. If we run this step without doing parallel, we would take at most $0.9\times 7\times (365\times5)\times(8+11+11)=344925$ seconds. We provide intermeidate results "BIC_land.csv", "BIC_ocean.csv", "BICd_land.csv", and "BICd_ocean.csv".)
6. Calculate $v^2(L_i,l_j)$ under $Q_l=36$ and $Q_o=68$ and plot Figure S14(d) (The major computational time in this step is to calculate the inverse SHT with $Q=36$ and $68$ for each $r$ and $t$. It will take about $0.8\times 7\times 365\times 5=10220$ seconds without the parallel. We provide the intermediate result "v2hat_do.csv".)
7. Do the real-valued transformation to SHT coefficients so that they are real values (about 4.2 seconds)
8. Test the normality of coefficients (about 11.5 seconds)
9. Choose the order of TGH autoregressive model for (non-Gaussian) time series at each coefficient $q$ using BIC and plot Figure S14(e) (about 166635.3 seconds without the parallel. We provide intermediate results "bicp_noTukey.csv" and "bicp_Tukey.csv".)
10. Model the temporal dependence structure using a Tukey g-and-h autoregressive model with order $P=1$ and plot Figure S14(f) (about 16791.3 seconds without the parallel. We provide the intermediate result "Tukeyres.csv".)
11. Model the spatial dependence by evaluating the covariance matrix of (read-valued and Gaussianized) SHT coefficients (about 80.4 seconds)
12. Calculate the covaraince matrix $\check{**U**}$  (about 0.3 seconds. We provide the intermediate result "U.mat".)
13. Generate $R'=7$ ensembles of daily emulations using 4 cores (about 3636.2 seconds)
14. Calculate $I_{uq}$ values using 4 cores and plot Figures 6(a), 6(b), and 6(e) (about 2352.8 seconds. We provide intermediate results "Iuq_noTukey.csv" and "Iuq_Tukey.csv".)
15. Calculate $WD_S$ values using 4 cores and plot Figures 6(c), 6(d), and 6(f) (about 129.9 seconds. We provide intermediate results "WD_time_noTukey.csv" and "WD_time_Tukey.csv")
16. Plot Figures 7 and S15 (about 0 second)


