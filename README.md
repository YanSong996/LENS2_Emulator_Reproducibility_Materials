# Reproducibility Instruction
This file documents the artifacts associated with the article (i.e., the data and code supporting the computational findings) and describes how to reproduce the findings. To reproduce all figures and results in the main manuscript and Supplementary Materials, please download and extract the provided "LENS2_Emulator_Reproducibility_Materials.zip", and then set your working directory to the folder "LENS2_Emulator_Reproducibility_Materials". 


## Step 1. Download and process the data in sub-repository "LENS2_Data"
Please follow the README.md file in the sub-repository "LENS2_Data" to download the raw data, process and store them using file “Data_treatment.R” in the "LENS2_Data". Processing montly and annual data together for each ensemble takes about 63 seconds. Processing daily data for each ensemble takes about 32.5 seconds.

## Step 2. Reproduce each section and figure sequentially
### Figure 1 in Section 2 (and Figure S1 in Section S2)
Figures 1 and S1 are used to illustrate the characteristics of surface temperature simulations in different scales.
1. Load the annual data and necessary information (about 5.2 seconds)
2. Calculate annual ensemble mean to plot Figure 1(a) (about 1.9 seconds)
3. Calculate annual ensemble sd to plot Figure 1(b) (about 3.8 seconds)
4. Plot Figure 1(c) (about 0.6 seconds)
5. Save R=7 annual time series at two grid points for ploting Figure 1(e) (about 0.4 seconds)
6. Calculate skewness and kurtosis for residuals of annaul data to plot Figures S1(a) and S1(d) (about 13.9 seconds)
7. Load the monthly data and necessary information (about 71.6 seconds)
8. Plot Figure 1(d) (about 0.4 seconds)
9. Save R=7 monthly time series at two grid points for ploting Figure 1(f) (about 0.5 seconds)
10. Calculate skewness and kurtosis for residuals of monthly data to plot Figures S1(b) and S1(e) (about 100.8 seconds)
11. Load the monthly data and necessary information (about 104.1 seconds)
12. Save R=7 daily time series at two grid points for ploting Figure 1(g) (about 0.4 seconds)
13. Calculate skewness and kurtosis for residuals of daily data to plot Figures S1(c) and S1(f) (about 180.7 seconds)
14. Plot Figures 1(e)-1(g) (about 3.2 seconds)
15. Plot Figures S1(a)-S1(c) (about 1.8 seconds) and Figures S1(d)-S1(f) (about 1.7 seconds)




