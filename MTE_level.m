% 05/27/2022 Yan Liu
% Plot level sets of integral of estimated MTE as a function of (Z1, Z2)

% Data input and preparation
clear all
load IFLS2000_main_trim.mat

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
pbar = mean(p);

X1 = [ones(n,1) X].*p;
X0 = [ones(n,1) X].*(1-p);
Z21 = Z2.*p;
Z20 = Z2.*(1-p);

% Parametric estimation of MTE
W = [X0 Z20 X1 Z21 p.^2-p]; % Second-order polynomial in propensity score
theta = (W'*W)\(W'*Y);
betaX = theta(24:45)-theta(1:22);
betaZ2 = theta(46)-theta(23);
alpha2 = theta(47);

% Create grids for contour plot
u = linspace(0,1,100);
z2all = linspace(min(Z2),max(Z2),100);
[U,Z2all] = meshgrid(u,z2all);

% Calculate integral of estimated MTE evaluated at mean values of X
Xbar = mean(X);
intMTE = [1 Xbar]*betaX*(U-pbar)+Z2all*betaZ2.*(U-pbar)+(U.^2-pbar^2-(U-pbar))*alpha2;

h = figure('Color','white');
contour(Z2all,U-pbar,intMTE,'ShowText','on')
title('Panel B: Level Sets of Outcome Changes');
xlabel('Distance to School (in km)');
ylabel('Change in Propensity Score');
saveas(h,'MTE_level','epsc');