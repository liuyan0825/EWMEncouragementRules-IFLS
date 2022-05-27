% 05/27/2022 Yan Liu
% Plot level sets of changes in estimated propensity scores as a function of (Z1, Z2)

% Data input and preparation
clear all
load IFLS2000_main_trim.mat
load propensity_coefs.mat

data = IFLS2000_main_trim;
n = size(data,1);
Z1 = data.exp/1000;
Z2 = data.dist_sec;
X = [data.dist_health data.ar09 data.ar09.^2 data.une_p data.ele_p data.sec_p ...
    data.une_m data.ele_m data.sec_m data.rural data.n_sumatra data.w_sumatra ...
    data.s_sumatra data.lampung data.c_java data.yogyakarta data.e_java ...
    data.bali data.w_nussa_tengara data.s_kalimanthan data.s_sulawesi];

dX = size(X,2);
gamma0 = gamma(1);
gammaX = gamma(2:dX+1);
gammaZ1 = gamma(dX+2);
gammaXZ1 = gamma(dX+3:end-dX-2);
gammaZ2 = gamma(end-dX-1);
gammaXZ2 = gamma(end-dX:end-1);
gammaZ12 = gamma(end);

% Create grids for contour plot
z1all = linspace(min(Z1),max(Z1),100);
z2all = linspace(min(Z2),max(Z2),100);
[Z1all,Z2all] = meshgrid(z1all,z2all);

% Calculate changes in estimated propensity scores evaluated at mean values of X
Xbar = mean(X);
tX = gamma0+Xbar*gammaX+gammaZ1*Z1all+Xbar*gammaXZ1*Z1all+gammaZ2*Z2all...
    +Xbar*gammaXZ2*Z2all+gammaZ12*Z1all.*Z2all;
p = exp(tX)./(1+exp(tX)); % propensity score before manipulation
tXM = gamma0+Xbar*gammaX+gammaZ1*0+Xbar*gammaXZ1*0+gammaZ2*Z2all...
    +Xbar*gammaXZ2*Z2all+gammaZ12*0.*Z2all;
pM = exp(tXM)./(1+exp(tXM)); % propensity score after manipulation (full tuition waiver)

% Prepare the density plot of the covariates
Zr = [Z1 Z2];
[Zr, Ind] = sortrows(Zr);
Zu = unique(Zr,'rows');
exp_u = Zu(:,1);
dist_sec_u = Zu(:,2);
nu = size(Zu,1);
nw = zeros(nu,1);
jj = 1;
for j = 1:n
    if ~(sum(Zr(j,:)~=Zu(jj,:)))
        nw(jj) = nw(jj) + 1;
    else
        jj = jj+1;
        nw(jj) = nw(jj) + 1;
    end
end

h = figure('Color','white');
contour(Z2all,Z1all,pM-p,'ShowText','on')
hold on
title('Panel A: Level Sets of Propensity Score Changes');
xlabel('Distance to School (in km)');
ylabel('Fees per Continuing Student (in 1000 Rupiah)');

% Overlay the density plot of the covariates
scatter(dist_sec_u,exp_u,nw,'MarkerEdgeColor','none','MarkerFaceColor','black');
saveas(h,'propensity_level.png');