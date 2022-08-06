% 08/06/2022 Yan Liu
% Calculate the feasible EWM encouragement rule
% Replicate Table 1 and Figure 1

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
X = [ones(n,1) X];
Z = [Z1 XZ1 Z2 XZ2 Z12];
p = data.phat;

X1 = X.*p;
X0 = X.*(1-p);
Z21 = Z2.*p;
Z20 = Z2.*(1-p);

% Parametric estimation of MTE
p2 = p.^2-p; % Second-order polynomial in propensity score
W = [X0 Z20 X1 Z21 p2];
theta = (W'*W)\(W'*Y);
beta0 = theta(1:23);
beta1 = theta(24:46);
alpha2 = theta(47);

% Calculate propensity scores after manipulation
% Case 1: alpha=2.5 (median tuition fee)
M1 = (Z1-2.5).*(Z1>=2.5);
XM1 = X(:,2:end).*M1;
M1Z2 = M1.*Z2;
ZM1 = [M1 XM1 Z2 XZ2 M1Z2];
p_M1 = predictp(X,ZM1,gamma);
p2_M1 = p_M1.^2-p_M1;

% Case 2: alpha=22.25 (maximum tuition fee)
M2 = (Z1-22.25).*(Z1>=22.25);
XM2 = X(:,2:end).*M2;
M2Z2 = M2.*Z2;
ZM2 = [M2 XM2 Z2 XZ2 M2Z2];
p_M2 = predictp(X,ZM2,gamma);
p2_M2 = p_M2.^2-p_M2;

% Calculte operator kernel in social welfare criterion
g_M1 = [X Z2]*(beta1-beta0).*(p_M1-p)+(p2_M1-p2)*alpha2;
g_M2 = [X Z2]*(beta1-beta0).*(p_M2-p)+(p2_M2-p2)*alpha2;

% Specify color and transparency for patches
color = true;
if (color)
    patch1_edge  = [0.6 0.4 0.1];
    patch1_color = [0.93 0.79 0.57];
    patch2_edge  = [0.3 0.4 0.3];
    patch2_color = [0.6 0.7 0.6];
    Face2Alpha1    = 0.7;
    Face2Alpha2    = 0.3;
    savefilename ='EWMrule_IFLS_color';
else
    patch1_edge  = [0.2 0.2 0.2];
    patch1_color = [0.5 0.5 0.5];
    patch2_edge  = [0.5 0.5 0.5];
    patch2_color = [0.8 0.8 0.8];
    Face2Alpha1    = 0.7;
    Face2Alpha2    = 0.3;
    savefilename = 'EWMrule_IFLS';
end

% Setup CPLEX optimization options
opt = cplexoptimset('cplex');
opt.parallel = 1;
opt.threads = 4;
opt.simplex.tolerances.feasibility = 1e-08;
opt.mip.tolerances.integrality = 1e-08;
opt.mip.strategy.variableselect = 3;
opt.mip.strategy.nodeselect = 2;
opt.mip.strategy.lbheur = 1;
opt.mip.limits.cutpasses = -1;
opt.display = 'off';

% create regressor matrix V from the data
V  = [Z1 Z2];

% Rescale covariates to [-1,1]
Vscale = ones(n,1)*max(abs(V));
V = [ones(n,1) V./Vscale];
k = size(V,2);

[V, Ind] = sortrows(V);
g_M1 = g_M1(Ind);
g_M2 = g_M2(Ind);
Vu = unique(V,'rows');
nu  = size(Vu,1); % number of unique covariate vectors
gu_M1 = zeros(nu,1);
gu_M2 = zeros(nu,1);
jj = 1;
for j = 1:n
    if ~(sum(V(j,:)~=Vu(jj,:)))
        gu_M1(jj) = gu_M1(jj) + g_M1(j);
        gu_M2(jj) = gu_M2(jj) + g_M2(j);
    else
        jj = jj+1;
        gu_M1(jj) = gu_M1(jj) + g_M1(j);
        gu_M2(jj) = gu_M2(jj) + g_M2(j);
    end
end

% add explicit monotonicity constraints
% Vu is ordered by increasing fees per continuing student (Vu(:,2))
% and then by increasing distance to school (Vu(:,3))
% For each number of children we impose treatment set inclusion
samefee = (Vu(1:end-1,2)==Vu(2:end,2));
% decreasing in distance to school
Mineq_d = [diag(samefee) zeros(nu-1,1)] + [zeros(nu-1,1) diag(-samefee)];
% increasing in distance to school
Mineq_i = [diag(-samefee) zeros(nu-1,1)] + [zeros(nu-1,1) diag(samefee)];

f_M1 = [zeros(k,1); -gu_M1]; % objective function coefficients for case 1
f_M2 = [zeros(k,1); -gu_M2]; % objective function coefficients for case 2
B = 1; % bounds on coefficients
C = B*sum(abs(Vu),2); % maximum values of v'beta
minmargin = max(1,C)*(1e-8); % prevent non-integer numbers the integrality 
% constraint of integers from being counted as integers
Aineq_d = [[-Vu diag(C)]; [Vu -diag(C)]; [zeros(nu-1,k) Mineq_d]];
Aineq_i = [[-Vu diag(C)]; [Vu -diag(C)]; [zeros(nu-1,k) Mineq_i]];
bineq = [[C-minmargin];[-minmargin];minmargin(1:nu-1,:)];
lb = [-B*ones(k,1); zeros(nu,1)];
ub = [ B*ones(k,1);  ones(nu,1)];

% Variable type string
ctype = strcat(repmat('C',1,k),repmat('B',1,nu));

% Welfare maximization with treatment decreasing in distance to school
[sol_pd, v_pd] = cplexmilp(f_M1,Aineq_d,bineq,[],[],[],[],[],lb,ub,ctype,[],opt);
%                      with treatment increasing in distance to school
[sol_pi, v_pi] = cplexmilp(f_M1,Aineq_i,bineq,[],[],[],[],[],lb,ub,ctype,[],opt);
    
if (v_pd < v_pi)
    lambda_M1 = sol_pd(1:k,:);
    v_M1    = -v_pd;
else
    lambda_M1 = sol_pi(1:k,:);
    v_M1    = -v_pi;
end

% Calculate Panel A of Table 1
in_Ghat_M1 = (V*lambda_M1>0);
fprintf('Share of Eligible Population (FEWM):\n%.4f\n',mean(in_Ghat_M1));
fprintf('Est. Welfare Gain (FEWM):\n%.4f\n',mean(g_M1.*in_Ghat_M1));
fprintf('Avg. Change in Treatment Take-up (FEWM):\n%.4f\n',mean((p_M1(Ind)-p(Ind)).*in_Ghat_M1));
fprintf('PRTE (FEWM):\n%.4f\n',mean(g_M1.*in_Ghat_M1)/mean((p_M1(Ind)-p(Ind)).*in_Ghat_M1));
fprintf('Est. Welfare Gain (encourage everyone):\n%.4f\n',mean(g_M1));
fprintf('Avg. Change in Treatment Take-up (encourage everyone):\n%.4f\n',mean(p_M1-p));
fprintf('PRTE (encourage everyone):\n%.4f\n',mean(g_M1)/mean(p_M1-p));

% Welfare maximization with encouragement decreasing in distance to school
[sol_pd, v_pd] = cplexmilp(f_M2,Aineq_d,bineq,[],[],[],[],[],lb,ub,ctype,[],opt);
%                      with encouragement increasing in distance to school
[sol_pi, v_pi] = cplexmilp(f_M2,Aineq_i,bineq,[],[],[],[],[],lb,ub,ctype,[],opt);
    
if (v_pd < v_pi)
    lambda_M2 = sol_pd(1:k,:);
    v_M2    = -v_pd;
else
    lambda_M2 = sol_pi(1:k,:);
    v_M2    = -v_pi;
end

% Calculate Panel B of Table 1
in_Ghat_M2 = (V*lambda_M2>0);
fprintf('Share of Eligible Population (FEWM):\n%.4f\n',mean(in_Ghat_M2));
fprintf('Est. Welfare Gain (FEWM):\n%.4f\n',mean(g_M2.*in_Ghat_M2));
fprintf('Avg. Change in Treatment Take-up (FEWM):\n%.4f\n',mean((p_M2(Ind)-p(Ind)).*in_Ghat_M2));
fprintf('PRTE (FEWM):\n%.4f\n',mean(g_M2.*in_Ghat_M2)/mean((p_M2(Ind)-p(Ind)).*in_Ghat_M2));
fprintf('Est. Welfare Gain (encourage everyone):\n%.4f\n',mean(g_M2));
fprintf('Avg. Change in Treatment Take-up (encourage everyone):\n%.4f\n',mean(p_M2-p));
fprintf('PRTE (encourage everyone):\n%.4f\n',mean(g_M2)/mean(p_M2-p));

h = figure('Color','white');

% Calculate the cutoff of fees at different levels of distance to school
% for case 1
line_dist_sec = [min(Z2)-0.5:0.01:max(Z2)+0.5]';
line_exp_M1 = Vscale(1,1)*(-(lambda_M1(1) + line_dist_sec*lambda_M1(3)./Vscale(1,2))./lambda_M1(2));

if (lambda_M1(2)>0) % encouragement rule is increasing in fees
    select_M1 = (line_exp_M1 <= 25);
    patchM1_X = [line_dist_sec(select_M1); flipud(line_dist_sec(select_M1))];
    patchM1_Y = [max(line_exp_M1(select_M1),0); 25*ones(length(line_exp_M1(select_M1)),1)];
else               % encouragement rule is decreasing in fees   
    select_M1 = (line_exp_M1 >= 0);
    patchM1_X = [line_dist_sec(select_M1); flipud(line_dist_sec(select_M1))];
    patchM1_Y = [min(line_exp_M1(select_M1),25); zeros(length(line_exp_M1(select_M1)),1)];
end
patch(patchM1_X,patchM1_Y,patch2_color,'LineStyle','-','FaceAlpha',Face2Alpha2,'EdgeColor',patch2_edge);
hold on

% Calculate the cutoff of fees at different levels of distance to school
% for case 2
line_exp_M2 = Vscale(1,1)*(-(lambda_M2(1) + line_dist_sec*lambda_M2(3)./Vscale(1,2))./lambda_M2(2));

if (lambda_M2(2)>0) % encouragement rule is increasing in fees
    select_M2 = (line_exp_M2 <= 25);
    patchM2_X = [line_dist_sec(select_M2); flipud(line_dist_sec(select_M2))];
    patchM2_Y = [max(line_exp_M2(select_M2),0); 25*ones(length(line_exp_M2(select_M2)),1)];
else               % encouragement rule is decreasing in fees   
    select_M2 = (line_exp_M2 >= 0);
    patchM2_X = [line_dist_sec(select_M2); flipud(line_dist_sec(select_M2))];
    patchM2_Y = [min(line_exp_M2(select_M2),25); zeros(length(line_exp_M2(select_M2)),1)];
end
patch(patchM2_X,patchM2_Y,patch1_color,'LineStyle','-','FaceAlpha',Face2Alpha1,'EdgeColor',patch1_edge);
hold on

% Prepare the density plot of instruments
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

% Overlay the density plot of instruments
scatter(dist_sec_u,exp_u,nw,'MarkerEdgeColor','none','MarkerFaceColor','black');

% Plot labeling
axis([min(Z2)-0.5 max(Z2)+0.5 0 25]);
xlabel('Distance to School (in km)');
ylabel('Fees per Continuing Student (in 1000 Rupiah)');
ax = gca;
ax.YTick = [0,5,10,15,20,25];
ax.YTickLabel = {'0','5','10','15','20','25'};
legend('2500 Rupiah subsidy','22250 Rupiah subsidy','Population Density','Location','northeast');
saveas(h,savefilename,'epsc');