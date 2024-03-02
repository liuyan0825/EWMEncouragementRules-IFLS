% 1/20/2024 Yan Liu
% Estimate the propensity score p(X,Z)

% Data input and preparation
clear all
data = readtable('IFLS2000_main.csv');
n = size(data,1);
D = data.upsec;
Z1 = data.exp/1000;
Z2 = data.dist_sec;
Z12 = Z1.*Z2;
X = [data.ar09 data.ar09.^2 data.rural data.dist_health ...
    data.protestant data.catholic data.religion_other ...
    data.ele_p data.sec_p data.missing_p data.ele_m data.sec_m data.missing_m...
    data.n_sumatra data.w_sumatra data.s_sumatra data.lampung ...
    data.c_java data.yogyakarta data.e_java data.bali ...
    data.w_nussa_tengara data.s_kalimanthan data.s_sulawesi];
XZ1 = X.*Z1;
XZ2 = X.*Z2;

% Estimate propensity score by logistic regression
gamma = glmfit([X Z1 XZ1 Z2 XZ2 Z12],D,'binomial','link','logit');
save('propensity_coefs.mat','gamma')

% Calculate predicted propensity scores
t = [ones(n,1) X Z1 XZ1 Z2 XZ2 Z12]*gamma;
phat = exp(t)./(1+exp(t));

data.phat = phat;
IFLS2000_main = data;
save('IFLS2000_main.mat','IFLS2000_main')