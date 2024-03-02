% 1/20/2024 Yan Liu
% Estimate the coefficients in the partially linear model for E[Y|X,Z2,p(X,Z)]
% using the double residual regression procedure of Robinson (1988)

% Data input and preparation
clear all
load IFLS2000_main.mat

data = IFLS2000_main;
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
p = data.phat;
dX = size(X,2);
Xp = X.*p;
Z2p = Z2.*p;
h = 0.06; %bandwidth

% Step 1: fit a local linear regression of Y and each regressor on p
Xhat = zeros(n,dX+1);
Xphat = zeros(n,dX+1);
Yhat = zeros(n,1);
for i = 1:n
    W = diag(normpdf((p-p(i))/h));
    P = [ones(n,1) p-p(i)];
    for j = 1:dX
        thetaX = (P'*W*P)\(P'*W*X(:,j));
        Xhat(i,j) = thetaX(1);
        thetaXp = (P'*W*P)\(P'*W*Xp(:,j));
        Xphat(i,j) = thetaXp(1);
    end
    thetaZ2 = (P'*W*P)\(P'*W*Z2);
    Xhat(i,dX+1) = thetaZ2(1);
    thetaZ2p = (P'*W*P)\(P'*W*Z2p);
    Xphat(i,dX+1) = thetaZ2p(1);
    thetaY = (P'*W*P)\(P'*W*Y);
    Yhat(i) = thetaY(1);
end
            
% Step 2: generate residual for Y and each regressor
eY = Y-Yhat;
eX = [X Z2]-Xhat;
eXp = [Xp Z2p]-Xphat;
eW = [eX eXp];
            
% Step 3: regress e_Y on e_X
betae = (eW'*eW)\(eW'*eY);
save('ddr_coefs.mat','betae')