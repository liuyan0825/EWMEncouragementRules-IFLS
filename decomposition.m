% 12/12/2022 Yan Liu
% Plot level sets of changes in average treatment take-up and average treatment 
% effects for those induced to switch treatment status conditional on (Z1,Z2)
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
X = [data.ar09 data.ar09.^2 data.rural data.dist_health ...
    data.protestant data.catholic data.religion_other ...
    data.ele_p data.sec_p data.missing_p data.ele_m data.sec_m data.missing_m...
    data.n_sumatra data.w_sumatra data.s_sumatra data.lampung ...
    data.c_java data.yogyakarta data.e_java data.bali ...
    data.w_nussa_tengara data.s_kalimanthan data.s_sulawesi];
p = data.phat;

% Create grids for contour plot
z1all = linspace(min(Z1),max(Z1),100);
z2all = linspace(min(Z2),max(Z2),100);
[Z1all,Z2all] = meshgrid(z1all,z2all);

% Calculate changes in average treatment take-up conditional on (Z1,Z2)
dX = size(X,2);
gamma0 = gamma(1);
gammaX = gamma(2:dX+1);
gammaZ1 = gamma(dX+2);
gammaXZ1 = gamma(dX+3:end-dX-2);
gammaZ2 = gamma(end-dX-1);
gammaXZ2 = gamma(end-dX:end-1);
gammaZ12 = gamma(end);
tX = gamma0+kron(X*gammaX,ones(100))+kron(ones(n,1),gammaZ1*Z1all)...
    +kron(X*gammaXZ1,Z1all)+kron(ones(n,1),gammaZ2*Z2all)...
    +kron(X*gammaXZ2,Z2all)+kron(ones(n,1),gammaZ12*Z1all.*Z2all);
p0 = exp(tX)./(1+exp(tX)); % Propensity score before manipulation
tXM = gamma0+kron(X*gammaX,ones(100))+kron(ones(n,1),gammaZ1*zeros(100))...
    +kron(X*gammaXZ1,zeros(100))+kron(ones(n,1),gammaZ2*Z2all)...
    +kron(X*gammaXZ2,Z2all)+kron(ones(n,1),gammaZ12*zeros(100).*Z2all);
pM = exp(tXM)./(1+exp(tXM)); % Propensity score after manipulation (full tuition waiver)
J = kron(ones(1,n),eye(100));
TT = J*(pM-p0)/n;

% Parametric estimation of MTE
X1 = [ones(n,1) X].*p;
X0 = [ones(n,1) X].*(1-p);
Z21 = Z2.*p;
Z20 = Z2.*(1-p);
W = [X0 Z20 X1 Z21 p.^2-p]; % Second-order polynomial in propensity score
theta = (W'*W)\(W'*Y);
betaX = theta(27:51)-theta(1:25);
betaZ2 = theta(52)-theta(26);
alpha = theta(53);

% Calculate average treatment effects for those induced to switch treatment
% status conditional on (Z1,Z2)
numerator = kron([ones(n,1) X]*betaX,ones(100)).*(pM-p0)+...
      +kron(ones(n,1),Z2all*betaZ2).*(pM-p0)+(pM.^2-pM-(p0.^2-p0))*alpha;
TE = J*(numerator)/n./TT;
TE = min(max(TE,quantile(TE,0.025,'all')),quantile(TE,0.975,'all')); %Winsorize at 2.5% level

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
hold on
title('Panel B: Treatment Effects')
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