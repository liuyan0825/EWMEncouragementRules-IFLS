% 12/20/2021 Yan Liu
% Calculate the hybrid EWM encouragement rule
% Also calculate the hybrid EWM treatment rule by Sasaki and Ura (2020)

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
W = [X0 Z20 X1 Z21 p.^2];
theta = (W'*W)\(W'*Y);
beta0 = theta(1:23);
beta1 = theta(24:46);
alpha2 = theta(47);

% Calculate propensity scores after manipulation
M = (Z1-2.5).*(Z1>=2.5);
XM = X(:,2:end).*M;
MZ2 = M.*Z2;
ZM = [M XM Z2 XZ2 MZ2];
p_M = predictp(X,ZM,gamma);

% Calculte operator kernel in social welfare criterion
g_M = [X Z2]*(beta1-beta0).*(p_M-p)+alpha2*(p_M.^2-p.^2);
g = [X Z2]*(beta1-beta0)+alpha2;

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

% create regressor matrix X from the data
V  = [Z1 Z2];

% Rescale covariates to [-1,1]
Vscale = ones(n,1)*max(abs(V));
V = [ones(n,1) V./Vscale];
k = size(V,2);

[V, Ind] = sortrows(V);
g_M = g_M(Ind);
g = g(Ind);
Vu = unique(V,'rows');
nu  = size(Vu,1); % number of unique covariate vectors
gu_M = zeros(nu,1);
gu = zeros(nu,1);
jj = 1;
for j = 1:n
    if ~(sum(V(j,:)~=Vu(jj,:)))
        gu_M(jj) = gu_M(jj) + g_M(j);
        gu(jj) = gu(jj) + g(j);
    else
        jj = jj+1;
        gu_M(jj) = gu_M(jj) + g_M(j);
        gu(jj) = gu(jj) + g(j);
    end
end

% add explicit monotonicity constraints
% Xu is ordered by increasing fees per continuing student (Xu(:,2))
% and then by increasing distance to school (Xu(:,3))
% For each number of children we impose treatment set inclusion
samefee = (Vu(1:end-1,2)==Vu(2:end,2));
% decreasing in distance to school
Mineq_d = [diag(samefee) zeros(nu-1,1)] + [zeros(nu-1,1) diag(-samefee)];
% increasing in distance to school
Mineq_i = [diag(-samefee) zeros(nu-1,1)] + [zeros(nu-1,1) diag(samefee)];

f_M = [zeros(k,1); -gu_M]; % objective function coefficients
f = [zeros(k,1); -gu];
B = 1; % bounds on coefficients
C = B*sum(abs(Vu),2); % maximum values of x'beta
minmargin = max(1,C)*(1e-8); % prevent non-integer numbers the integrality constraint of integers from being counted as integers
Aineq_d = [[-Vu diag(C)]; [Vu -diag(C)]; [zeros(nu-1,k) Mineq_d]];
Aineq_i = [[-Vu diag(C)]; [Vu -diag(C)]; [zeros(nu-1,k) Mineq_i]];
bineq = [[C-minmargin];[-minmargin];minmargin(1:nu-1,:)];
lb = [-B*ones(k,1); zeros(nu,1)];
ub = [ B*ones(k,1);  ones(nu,1)];

% Variable type string
ctype = strcat(repmat('C',1,k),repmat('B',1,nu));

% Welfare maximization with treatment decreasing in distance to school
[sol_pd, v_pd] = cplexmilp(f_M,Aineq_d,bineq,[],[],[],[],[],lb,ub,ctype,[],opt);
%                      with treatment increasing in distance to school
[sol_pi, v_pi] = cplexmilp(f_M,Aineq_i,bineq,[],[],[],[],[],lb,ub,ctype,[],opt);
    
if (v_pd < v_pi)
    lambda_M = sol_pd(1:k,:);
    v_M    = -v_pd;
else
    lambda_M = sol_pi(1:k,:);
    v_M    = -v_pi;
end

in_Ghat_M = (V*lambda_M>0);
fprintf('Proportion treated:\n%.4f\n',mean(in_Ghat_M));
fprintf('Est. welfare gain:\n%.4f\n',mean(g_M.*in_Ghat_M));
fprintf('Est. welfare gain (treat everyone):\n%.4f\n',mean(g_M));

% Welfare maximization with treatment decreasing in distance to school
[sol_pd, v_pd] = cplexmilp(f,Aineq_d,bineq,[],[],[],[],[],lb,ub,ctype,[],opt);
%                      with treatment increasing in distance to school
[sol_pi, v_pi] = cplexmilp(f,Aineq_i,bineq,[],[],[],[],[],lb,ub,ctype,[],opt);
    
if (v_pd < v_pi)
    lambda = sol_pd(1:k,:);
    v    = -v_pd;
else
    lambda = sol_pi(1:k,:);
    v    = -v_pi;
end

in_Ghat = (V*lambda>0);
fprintf('Proportion treated:\n%.4f\n',mean(in_Ghat));
fprintf('Est. welfare gain:\n%.4f\n',mean(g.*in_Ghat));
fprintf('Est. welfare gain (treat everyone):\n%.4f\n',mean(g));

patch1_edge  = [0.6 0.4 0.1];
patch1_color = [0.93 0.79 0.57];
patch2_edge  = [0.3 0.4 0.3];
patch2_color = [0.6 0.7 0.6];
Face2Alpha1    = 0.7;
Face2Alpha2    = 0.3;

h = figure('Color','white');

% Calculate the cutoff of fees at different levels of distance to school
% for encouragement rules
line_dist_sec = [min(Z2)-0.5:0.01:max(Z2)+0.5]';
line_exp_M = Vscale(1,1)*(-(lambda_M(1) + line_dist_sec*lambda_M(3)./Vscale(1,2))./lambda_M(2));

if (lambda_M(2)>0) % treatment rule is increasing in fees
    select_M = (line_exp_M <= 25);
    patchM_X = [line_dist_sec(select_M); flipud(line_dist_sec(select_M))];
    patchM_Y = [max(line_exp_M(select_M),0); 25*ones(length(line_exp_M(select_M)),1)];
else               % treatment rule is decreasing in fees   
    select_M = (line_exp_M >= 0);
    patchM_X = [line_dist_sec(select_M); flipud(line_dist_sec(select_M))];
    patchM_Y = [min(line_exp_M(select_M),25); zeros(length(line_exp_M(select_M)),1)];
end
patch(patchM_X,patchM_Y,patch1_color,'LineStyle','-','FaceAlpha',Face2Alpha1,'EdgeColor',patch1_edge);
hold on

% Calculate the cutoff of fees at different levels of distance to school
% for treatment rules
line_exp = Vscale(1,1)*(-(lambda(1) + line_dist_sec*lambda(3)./Vscale(1,2))./lambda(2));

if (lambda(2)>0) % treatment rule is increasing in fees
    select = (line_exp <= 25);
    patch_X = [line_dist_sec(select); flipud(line_dist_sec(select))];
    patch_Y = [max(line_exp(select),0); 25*ones(length(line_exp(select)),1)];
else             % treatment rule is decreasing in fees            
    select = (line_exp >= 0);
    patch_X = [line_dist_sec(select); flipud(line_dist_sec(select))];
    patch_Y = [min(line_exp(select),25); zeros(length(line_exp(select)),1)];
end
patch(patch_X,patch_Y,patch2_color,'LineStyle','-','FaceAlpha',Face2Alpha2,'EdgeColor',patch2_edge);
hold on

% Prepare the density plot of the covariates
Xr = [Z1 Z2];
[Xr, Ind] = sortrows(Xr);
Xu = unique(Xr,'rows');
exp_u = Xu(:,1);
dist_sec_u = Xu(:,2);
nu = size(Xu,1);
nw = zeros(nu,1);
jj = 1;
for j = 1:n
    if ~(sum(Xr(j,:)~=Xu(jj,:)))
        nw(jj) = nw(jj) + 1;
    else
        jj = jj+1;
        nw(jj) = nw(jj) + 1;
    end
end

% overlay the density plot of the covariates
scatter(dist_sec_u,exp_u,nw,'MarkerEdgeColor','none','MarkerFaceColor','black');

% plot labeling
axis([min(Z2)-0.5 max(Z2)+0.5 0 25]);
xlabel('Distance to School (in km)');
ylabel('Fees per Continuing Student (in 1000 Rupiah)');
ax = gca;
ax.YTick = [0,5,10,15,20,25];
ax.YTickLabel = {'0','5','10','15','20','25'};
legend('Encouragement Rule','Treatment Rule','Population Density','Location','northeast');
saveas(h,'EWMrule_IFLS.png');