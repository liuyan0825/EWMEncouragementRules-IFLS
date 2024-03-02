% 1/20/2024 Yan Liu
% Select the bandwidth using leave-one-out cross-validation

function hcv(k)

% Data input and preparation
load IFLS2000_main.mat

data = IFLS2000_main;
n = size(data,1);
ind = 1:n;

Y = data.lwages;
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

hs = 0.01:0.01:1;
h = hs(k);
llr = zeros(n,1);

for i = 1:n
    Xhat = zeros(n-1,dX+1);
    Xphat = zeros(n-1,dX+1);
    Yhat = zeros(n-1,1);
    pi = p(ind~=i);
    
    % Step 1: fit a local linear regression of Y and each regressor on p
    for ii = 1:(n-1)
        W = diag(normpdf((pi-pi(ii))/h));
        P = [ones(n-1,1) pi-pi(ii)];
        for j = 1:dX
            thetaX = (P'*W*P)\(P'*W*X(ind~=i,j));
            Xhat(ii,j) = thetaX(1);
            thetaXp = (P'*W*P)\(P'*W*Xp(ind~=i,j));
            Xphat(ii,j) = thetaXp(1);
        end
        thetaZ2 = (P'*W*P)\(P'*W*Z2(ind~=i));
        Xhat(ii,dX+1) = thetaZ2(1);
        thetaZ2p = (P'*W*P)\(P'*W*Z2p(ind~=i));
        Xphat(ii,dX+1) = thetaZ2p(1);
        thetaY = (P'*W*P)\(P'*W*Y(ind~=i));
        Yhat(ii) = thetaY(1);
    end
            
    % Step 2: generate residual for Y and each regressor
    eY = Y(ind~=i)-Yhat;
    eX = [X(ind~=i,:) Z2(ind~=i)]-Xhat;
    eXp = [Xp(ind~=i,:) Z2p(ind~=i)]-Xphat;
    eW = [eX eXp];
            
    % Step 3: regress e_Y on e_X
    betae = (eW'*eW)\(eW'*eY);
    V = Y(ind~=i)-[X(ind~=i,:) Z2(ind~=i) Xp(ind~=i,:) Z2p(ind~=i)]*betae;
            
    % Step 4: fit a local linear regression of V on p
    W = diag(normpdf((pi-p(i))/h));
    P = [ones(n-1,1) pi-p(i)];
    theta = (P'*W*P)\(P'*W*V);
    llr(i) = [X(i,:) Z2(i) Xp(i,:) Z2p(i)]*betae+theta(1);
end
MSE = sum((Y-llr).^2)/n;

save(['MSE_h',num2str(h),'.mat'],'MSE')