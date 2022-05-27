# EWMEncouragementRules-IFLS

This repository contains data and Matlab codes for the second empirical application in "Policy Learning under Endogeneity Using Instrumental Variables". It calculates propensity scores and hybrid EWM encouragement rules, using data from the third wave of the Indonesia Family Life Survey (IFLS).

**Data source**: https://www.rand.org/well-being/social-and-behavioral-policy/data/FLS/IFLS/ifls3.html

**File description**:
1. propensity.m: Estimate the propensity score p(X,Z) and trim observations for p(X,Z) to have common support for treated and untreated. The coefficients from the logistic regression are saved in "propensity_coefs.mat". The trimmed data is saved in "IFLS2000_main_trim.mat".
2. EWM.m: Calculate the hybrid EWM encouragement rule for two tuition subsidy levels (median and maximum tuition fee) and replicate Figure 5.1 and Table 5.1. The optimization problem is solved with [IBM ILOG CPLEX Optimization Studio](https://www.ibm.com/products/ilog-cplex-optimization-studio). Part of the codes are adapted from [replication files](https://www.econometricsociety.org/sites/default/files/13288_Data_and_Programs.zip) for [Kitagawa and Tetenov (2018)](https://onlinelibrary.wiley.com/doi/abs/10.3982/ECTA13288).
3. propensity_level.m: Plot the level sets of changes in estimated propensity scores under a full tuition waiver (Panel A of Figure 5.2).
4. MTE_level.m: Plot the level sets of the integral of the estimated MTE (Panel B of Figure 5.2).
