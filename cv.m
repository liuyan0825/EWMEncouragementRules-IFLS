% 12/12/2022 Yan Liu
% Select the order for global polynomial specification using leave-one-out
% cross-validation

% Data input and preparation
warning('off','all')
clear all
load IFLS2000_main_trim.mat

data = IFLS2000_main_trim;
n = size(data,1);

Y = data.lwages;
D = data.upsec;
Z1 = data.exp/1000;
Z2 = data.dist_sec;
X = [data.ar09 data.ar09.^2 data.rural data.dist_health ...
    data.protestant data.catholic data.religion_other ...
    data.ele_p data.sec_p data.missing_p data.ele_m data.sec_m data.missing_m...
    data.n_sumatra data.w_sumatra data.s_sumatra data.lampung ...
    data.c_java data.yogyakarta data.e_java data.bali ...
    data.w_nussa_tengara data.s_kalimanthan data.s_sulawesi];
X = [ones(n,1) X];
p = data.phat;

X1 = X.*p;
X0 = X.*(1-p);
Z21 = Z2.*p;
Z20 = Z2.*(1-p);

% Compute mean squared error for polynomials of orders 2-15
MSE = zeros(14,1);
for k = 2:15
    poly = zeros(n,k-1);
    for j = 2:k
        poly(:,j-1) = p.^j-p;
    end
    W = [X0 Z20 X1 Z21 poly];
    MSE(k-1) = crossval('mse',W,Y,'Predfun',@regf,'Leaveout',1);
end