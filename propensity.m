% 05/27/2022 Yan Liu
% Estimate the propensity score p(X,Z)

% Data input and preparation
clear all
data = readtable('IFLS2000_main.csv');
n = size(data,1);
Y = data.lwages;
D = data.upsec;
Z1 = data.exp/1000;
Z2 = data.dist_sec;
Z12 = Z1.*Z2;
X = [data.dist_health data.ar09 data.ar09.^2 data.une_p data.ele_p data.sec_p ...
    data.une_m data.ele_m data.sec_m data.rural data.n_sumatra data.w_sumatra ...
    data.s_sumatra data.lampung data.c_java data.yogyakarta data.e_java ...
    data.bali data.w_nussa_tengara data.s_kalimanthan data.s_sulawesi];
XZ1 = X.*Z1;
XZ2 = X.*Z2;

% Estimate propensity score by logistic regression
gamma = glmfit([X Z1 XZ1 Z2 XZ2 Z12],D,'binomial','link','logit');
save('propensity_coefs.mat','gamma')

% Calculate predicted propensity scores
t = [ones(n,1) X Z1 XZ1 Z2 XZ2 Z12]*gamma;
phat = exp(t)./(1+exp(t));
select1 = (D==1);
select0 = (D==0);
p1 = phat(select1);
p0 = phat(select0);

% Trim observations to have common support for treated and untreated
select = (phat>=max(min(p0),min(p1)) & phat<=min(max(p0),max(p1)));
data.phat = phat;
IFLS2000_main_trim = data(select,:);
save('IFLS2000_main_trim.mat','IFLS2000_main_trim')
[min(IFLS2000_main_trim.phat) max(IFLS2000_main_trim.phat)]