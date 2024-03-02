% 1/20/2024 Yan Liu
% Plot level sets of changes in treatment take-up and PRTE when going 
% from the status quo to a full tuition waiver for individuals with 
% different values of (Z1,Z2) and the median value of X 
% Replicate Figure 2

% Data input and preparation
clear all
load IFLS2000_main.mat
load propensity_coefs.mat
load 'ddr_coefs.mat'

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
Xbar = median(X);
h = 0.06; %bandwidth

% Create grids for contour plot
z1all = linspace(min(Z1),max(Z1),100);
z2all = linspace(min(Z2),max(Z2),100);
[Z1all,Z2all] = meshgrid(z1all,z2all);

% Calculate changes in treatment take-up evaluated at (Xbar,Z1all,Z2all)
dX = size(X,2);
gamma0 = gamma(1);
gammaX = gamma(2:dX+1);
gammaZ1 = gamma(dX+2);
gammaXZ1 = gamma(dX+3:end-dX-2);
gammaZ2 = gamma(end-dX-1);
gammaXZ2 = gamma(end-dX:end-1);
gammaZ12 = gamma(end);
tX = gamma0+Xbar*gammaX*ones(100)+gammaZ1*Z1all+Xbar*gammaXZ1*Z1all...
    +gammaZ2*Z2all+Xbar*gammaXZ2*Z2all+gammaZ12*Z1all.*Z2all;
p0 = exp(tX)./(1+exp(tX)); % Propensity score before manipulation
tXalpha = gamma0+Xbar*gammaX*ones(100)+gammaZ2*Z2all+Xbar*gammaXZ2*Z2all;
palpha = exp(tXalpha)./(1+exp(tXalpha)); % Propensity score after manipulation (full tuition waiver)
TT = palpha-p0;

% Calculate PRTE evaluated at (Xbar,Z1all,Z2all)
Xp = X.*p;
Z2p = Z2.*p;
Vhat = Y-[X Z2 Xp Z2p]*betae;

llr = zeros(100);
for i = 1:100
    for j = 1:100
        W = diag(normpdf((p-p0(i,j))/h));
        P = [ones(n,1) p-p0(i,j)];
        theta = (P'*W*P)\(P'*W*Vhat);
        llr(i,j) = theta(1);
    end
end
llr_alpha = zeros(100);
for i = 1:100
    for j = 1:100
        W = diag(normpdf((p-palpha(i,j))/h));
        P = [ones(n,1) p-palpha(i,j)];
        theta = (P'*W*P)\(P'*W*Vhat);
        llr_alpha(i,j) = theta(1);
    end
end

betaX = betae(26:49);
betaZ2 = betae(50);
TE = Xbar*betaX*ones(100)+Z2all*betaZ2+(llr_alpha-llr)./TT;

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
set(gcf,'Position',[0 0 585 225]);

subplot(1,2,1)
[C1,h1] = contour(Z2all,Z1all,TT,'ShowText','on','LineWidth',1);
clabel(C1,h1,'FontSize',6)
hold on
title('Panel A: Changes in Treatment Take-up')
xlabel('Distance to School (in km)');
ylabel('Fees per Continuing Student (in 1000 Rupiah)')
set(gca,'FontSize',6)
set(gca,'Position',[0.04 0.12 0.44 0.83])
set(gca,'YTick',[0,5,10,15,20,25])
set(gca,'YTickLabel',{'0','5','10','15','20','25'})

% Overlay the density plot of the covariates
scatter(dist_sec_u,exp_u,nw*0.4,'MarkerEdgeColor','none','MarkerFaceColor','black');

subplot(1,2,2)
[C2,h2] = contour(Z2all,Z1all,TE,'LineWidth',1);
h2.LevelStep = 0.4;
hold on
title('Panel B: PRTE')
xlabel('Distance to School (in km)')
ylabel('Fees per Continuing Student (in 1000 Rupiah)')
set(gca,'FontSize',6)
set(gca,'Position',[0.53 0.12 0.44 0.83])
set(gca,'YTick',[0,5,10,15,20,25])
set(gca,'YTickLabel',{'0','5','10','15','20','25'})

% Overlay the density plot of the covariates
scatter(dist_sec_u,exp_u,nw*0.4,'MarkerEdgeColor','none','MarkerFaceColor','black');
contourLegend(h2)
saveas(h,'decompose','epsc');