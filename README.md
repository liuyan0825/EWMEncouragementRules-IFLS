# EWMEncouragementRules-IFLS

This repository contains data and analysis code for the empirical application in "Policy Learning under Endogeneity Using Instrumental Variables".

**Dataset: IFLS2000_main.csv**

This data comes from [the third wave of the Indonesian Family Life Survey (IFLS3)](https://www.rand.org/well-being/social-and-behavioral-policy/data/FLS/IFLS/ifls3.html).

- **Variable “pidlink”:** tracking id

- **Variable “hhid00”:** household id in 2000

- **Variable “commid00”:** community id in 2000

- **Variable “lwages”:** log hourly wages constructed using data from Book 3A, Section TK

- **Variable “edu”:** years of education constructed using data from Book 3A, Section DL

- **Variable “upsec”:** indicator of attendance of upper secondary school or higher constructed from “edu”

- **Variable “exp”:** tuition fee of the cheapest secondary school using data from School Questionnaire, Section E

- **Variable “dist_sec”:** distance to the nearest secondary school using data from Service Availability Roster (SAR)

- **Variable “dist_health”:** distance to the nearest health post using data from Service Availability Roster (SAR)

- **Variable “ar09”:** age from IFLS3 Control Book (Book K), Section AR

- **Variable “ar15”:** religion from IFLS3 Control Book (Book K), Section AR

- **Variables “protestant”, “catholic”, “muslim”, “religion_other”:** religion dummies constructed from “ar15”

- **Variables “une_p”, “ele_p”, “sec_p”, “missing_p”, “une_m”, “ele_m”, “sec_m”, “missing_m”:** parental education constructed using data from Book 3B, Section BA

- **Variable “sc05”:** area (urban/rural) of residence from IFLS3 Control Book (Book K), Section SC

- **Variable “rural”:** indicator of rural residence constructed from “sc05”

- **Variable “sc01”:** province of residence from IFLS3 Control Book (Book K), Section SC

- **Variables “n_sumatra”, “w_sumatra”, “s_sumatra”, “lampung”, “jakarta”, “c_java”, “yogyakarta”, “e_java”, “bali”, “w_nussa_tengara”, “s_kalimanthan”, “s_sulawesi”:** province dummies constructed from “sc01”

**Stata Code: IFLS_construction.do**

This code generates the dataset “IFLS2000_main.csv” from the original data files “b3a_dl1.dta”, “b3a_tk1.dta”, “b3a_tk2.dta”, “b3b_ba0.dta”, “bk_ar1.dta”, “bk_sc.dta”, “sar.dta”, “schl_e.dta”.

**Stata code: samplestat.do**

This code calculates sample averages for the main variables used in the analysis.

**MATLAB Code: propensity.m**

This code estimates the propensity score p(X,Z) and trims observations for p(X,Z) to have common support for treated and untreated. The coefficients from the logistic regression are saved in "propensity_coefs.mat". The trimmed data is saved in "IFLS2000_main_trim.mat".

**MATLAB Code: EWM.m**

This code calculates the feasible EWM encouragement rule and the budget-constrained EWM encouragement rule for two tuition subsidy levels (median and maximum tuition fee) and replicate Figure 1 and Table 1. The optimization problem is solved with [IBM ILOG CPLEX Optimization Studio](https://www.ibm.com/products/ilog-cplex-optimization-studio). Part of the codes are adapted from [replication files](https://www.econometricsociety.org/publications/econometrica/2018/03/01/who-should-be-treated-empirical-welfare-maximization-methods/supp/13288_Data_and_Programs.zip) for [Kitagawa and Tetenov (2018)](https://onlinelibrary.wiley.com/doi/abs/10.3982/ECTA13288).

**MATLAB Code: decomposition.m**

This code plots the level sets of changes in average treatment take-up and average treatment effects for those induced to switch treatment status conditional on (Z1,Z2) when going from the status quo to a full tuition waiver (Figure 2). Function [contourLegend.m](https://www.mathworks.com/matlabcentral/fileexchange/115120-legend-for-contour-plots) is required.

**MATLAB Code: specification.m**

This code compares the polynomial and local quadratic specifications for the MTE.

**MATLAB Code: cv.m**

This code selects the order for the polynomial specification using leave-one-out cross-validation.
