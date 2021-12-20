% 12/20/2021 Yan Liu
% Estimate the propensity score p(X,Z) and plot its conditional and
% unconditional supports

% Data input and preparation
clear all
load IFLS2000_main.mat

data = IFLS2000main;
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

% Create grids for heatmap
P = 0:0.05:1;
dX = size(X,2);
gammaX = gamma(2:dX+1);
gammaZ2 = gamma(end-dX-1);
gammaXZ2 = gamma(end-dX:end-1);
tX = X*gammaX+Z2*gammaZ2+XZ2*gammaXZ2;
T = linspace(min(tX),max(tX),21);

% Compute conditional probability density estimates of phat
f = zeros(20,20);
for i = 1:n
    for j = 1:20
        for k = 1:19
           if (phat(i)>=P(j)) && (phat(i)<P(j+1)) && (tX(i)>=T(k)) && (tX(i)<T(k+1))
              f(j,k) =  f(j,k)+1/n;
           end
        end
        if (phat(i)>=P(j)) && (phat(i)<P(j+1)) && (tX(i)>=T(20)) && (tX(i)<=T(21))
            f(j,20) =  f(j,20)+1/n;
        end
    end
end

f = f./sum(f);

% Plot conditional probability density estimates of phat in a heatmap
Pmid = 0.025:0.05:0.975;
Tstep = (max(tX)-min(tX))/20;
Tmid = min(tX)+Tstep/2:Tstep:max(tX)-Tstep/2;
[x,p] = meshgrid(Tmid,Pmid);

h1 = figure('Color','white');
surface(x,p,f)
colorbar
xlabel('$X$','interpreter','latex')
ylabel('$P$','interpreter','latex')
saveas(h1,'propensity_con_IFLS.png');

% Plot unconditional probability density estimates of phat in a histogram
select1 = (D==1);
select0 = (D==0);
p1 = phat(select1);
p0 = phat(select0);
Bins = linspace(0,1,21);
y1 = hist(p1,Bins)/n;   
y0 = hist(p0,Bins)/n;

h2 = figure('Color','white');
bar(Bins, [y1;y0]');
xlabel('$P$','interpreter','latex')
ylabel('$f(P)$','interpreter','latex')
legend('$D=1$ (Upper Secondary or Higher)','$D=0$ (Less than Upper Secondary)','interpreter','latex')
saveas(h2,'propensity_unc_IFLS.png');

% Trim observations to have common support for treated and untreated
select = (phat>=max(min(p0),min(p1)) & phat<=min(max(p0),max(p1)));
data.phat = phat;
IFLS2000_main_trim = data(select,:);
save('IFLS2000_main_trim.mat','IFLS2000_main_trim')
[min(IFLS2000_main_trim.phat) max(IFLS2000_main_trim.phat)]