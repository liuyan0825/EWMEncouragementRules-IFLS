% 08/06/2022 Yan Liu
% Plot level sets of changes in average treatment take-up and changes in 
% average treatment effects for those induced to switch treatment status
% conditional on (X,Z) evaluated at the mean value of X 
% when going from the status quo to a full tuition waiver   
% Replicate Figure 2

% Data input and preparation
clear all
load IFLS2000_main_trim.mat
load propensity_coefs.mat

data = IFLS2000_main_trim;
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
Z = [Z1 XZ1 Z2 XZ2 Z12];
p = data.phat;

% Create grids for contour plot
z1all = linspace(min(Z1),max(Z1),100);
z2all = linspace(min(Z2),max(Z2),100);
[Z1all,Z2all] = meshgrid(z1all,z2all);

% Calculate changes in estimated propensity scores evaluated at mean value of X
dX = size(X,2);
gamma0 = gamma(1);
gammaX = gamma(2:dX+1);
gammaZ1 = gamma(dX+2);
gammaXZ1 = gamma(dX+3:end-dX-2);
gammaZ2 = gamma(end-dX-1);
gammaXZ2 = gamma(end-dX:end-1);
gammaZ12 = gamma(end);
Xbar = mean(X);
tX = gamma0+Xbar*gammaX+gammaZ1*Z1all+Xbar*gammaXZ1*Z1all+gammaZ2*Z2all...
    +Xbar*gammaXZ2*Z2all+gammaZ12*Z1all.*Z2all;
p0 = exp(tX)./(1+exp(tX)); % propensity score before manipulation
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

% Parametric estimation of MTE
X1 = [ones(n,1) X].*p;
X0 = [ones(n,1) X].*(1-p);
Z21 = Z2.*p;
Z20 = Z2.*(1-p);
W = [X0 Z20 X1 Z21 p.^2-p]; % Second-order polynomial in propensity score
theta = (W'*W)\(W'*Y);
betaX = theta(24:45)-theta(1:22);
betaZ2 = theta(46)-theta(23);
alpha = theta(47);

% Calculate integral of estimated MTE evaluated at mean value of X
PRTE = ([1 Xbar]*betaX*(pM-p0)+Z2all*betaZ2.*(pM-p0)+(pM.^2-pM-(p0.^2-p0))*alpha)./(pM-p0);

h = figure('Color','white');
set(gcf,'Position',[0 0 585 225]);

subplot(1,2,1)
[C1,h1] = contour(Z2all,Z1all,pM-p0,'ShowText','on','LineWidth',1);
clabel(C1,h1,'FontSize',6)
hold on
title('Panel A: Changes in Average Treatment Take-up')
xlabel('Distance to School (in km)');
ylabel('Fees per Continuing Student (in 1000 Rupiah)')
set(gca,'FontSize',6)
set(gca,'Position',[0.04 0.12 0.44 0.83])
set(gca,'YTick',[0,5,10,15,20,25])
set(gca,'YTickLabel',{'0','5','10','15','20','25'})

% Overlay the density plot of the covariates
scatter(dist_sec_u,exp_u,nw*0.4,'MarkerEdgeColor','none','MarkerFaceColor','black');

subplot(1,2,2)
[C2,h2] = contour(Z2all,Z1all,PRTE,'ShowText','on','LineWidth',1);
clabel(C2,h2,'FontSize',6)
hold on
title('Panel B: Changes in Average Treatment Effects')
xlabel('Distance to School (in km)')
ylabel('Fees per Continuing Student (in 1000 Rupiah)')
set(gca,'FontSize',6)
set(gca,'Position',[0.53 0.12 0.44 0.83])
set(gca,'YTick',[0,5,10,15,20,25])
set(gca,'YTickLabel',{'0','5','10','15','20','25'})

% Overlay the density plot of the covariates
scatter(dist_sec_u,exp_u,nw*0.4,'MarkerEdgeColor','none','MarkerFaceColor','black');

saveas(h,'decompose','epsc');