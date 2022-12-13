% 12/12/2022 Yan Liu
% Calculate the feasible EWM encouragement rule and the budget-constrained 
% EWM encouragement rule
% Replicate Figure 1

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
X = [data.ar09 data.ar09.^2 data.rural data.dist_health ...
    data.protestant data.catholic data.religion_other ...
    data.ele_p data.sec_p data.missing_p data.ele_m data.sec_m data.missing_m...
    data.n_sumatra data.w_sumatra data.s_sumatra data.lampung ...
    data.c_java data.yogyakarta data.e_java data.bali ...
    data.w_nussa_tengara data.s_kalimanthan data.s_sulawesi];
XZ1 = X.*Z1;
XZ2 = X.*Z2;
Z = [Z1 XZ1 Z2 XZ2 Z12];
p = data.phat;

X1 = [ones(n,1) X].*p;
X0 = [ones(n,1) X].*(1-p);
Z21 = Z2.*p;
Z20 = Z2.*(1-p);

% Parametric estimation of MTE
p2 = p.^2-p; % Second-order polynomial in propensity score
W = [X0 Z20 X1 Z21 p2];
theta = (W'*W)\(W'*Y);
beta0 = theta(1:26);
beta1 = theta(27:52);
alpha2 = theta(53);

% Calculate propensity scores after manipulation
% Case 1: alpha=2.5 (median tuition fee)
M1 = (Z1-2.5).*(Z1>=2.5);
XM1 = X.*M1;
M1Z2 = M1.*Z2;
ZM1 = [M1 XM1 Z2 XZ2 M1Z2];
p_M1 = predictp([ones(n,1) X],ZM1,gamma);
p2_M1 = p_M1.^2-p_M1;

% Case 2: alpha=22.25 (maximum tuition fee)
M2 = (Z1-22.25).*(Z1>=22.25);
XM2 = X.*M2;
M2Z2 = M2.*Z2;
ZM2 = [M2 XM2 Z2 XZ2 M2Z2];
p_M2 = predictp([ones(n,1) X],ZM2,gamma);
p2_M2 = p_M2.^2-p_M2;

% Calculate operator kernel in social welfare criterion
g_M1 = [ones(n,1) X Z2]*(beta1-beta0).*(p_M1-p)+(p2_M1-p2)*alpha2;
g_M2 = [ones(n,1) X Z2]*(beta1-beta0).*(p_M2-p)+(p2_M2-p2)*alpha2;

% Calculate budget functions
B_M1 = (Z1-M1).*p_M1;
B_M2 = (Z1-M2).*p_M2;

% Report results for policies that encourage everyone
disp('Encourage everyone: alpha=2.5 (median tuition fee)')
fprintf('Est. welfare gain:\n%.4f\n',mean(g_M1));
fprintf('Average Propensity Score Changes:\n%.4f\n',mean(p_M1-p));
fprintf('PRTE:\n%.4f\n',mean(g_M1)/mean(p_M1-p));
fprintf('Est. total costs:\n%.4f\n',mean((Z1-M1).*p_M1));

disp('Encourage everyone: alpha=22.25 (maximum tuition fee)')
fprintf('Est. welfare gain:\n%.4f\n',mean(g_M2));
fprintf('Average Propensity Score Changes:\n%.4f\n',mean(p_M2-p));
fprintf('PRTE:\n%.4f\n',mean(g_M2)/mean(p_M2-p));
fprintf('Est. total costs:\n%.4f\n',mean((Z1-M2).*p_M2));

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

% Create regressor matrix V from the data
V  = [Z1 Z2];

% Rescale covariates to [-1,1]
Vscale = ones(n,1)*max(abs(V));
V = [ones(n,1) V./Vscale];
k = size(V,2);

f_M1 = [zeros(k,1); -g_M1]; % objective function coefficients for case 1
f_M2 = [zeros(k,1); -g_M2]; % objective function coefficients for case 2
B = 1; % bounds on coefficients
C = B*sum(abs(V),2); % maximum values of v'beta
minmargin = max(1,C)*(1e-8); % prevent non-integer numbers the integrality 
% constraint of integers from being counted as integers
Aineq = [[-V diag(C)]; [V -diag(C)]];
bineq = [[C-minmargin];[-minmargin]];
Aineq_M1 = [[-V diag(C)]; [V -diag(C)]; [zeros(1,k) B_M1']];
Aineq_M2 = [[-V diag(C)]; [V -diag(C)]; [zeros(1,k) B_M2']];
bineq_bc = [[C-minmargin];[-minmargin];0.28*n];
lb = [-B*ones(k,1); zeros(n,1)];
ub = [ B*ones(k,1);  ones(n,1)];

% Variable type string
ctype = strcat(repmat('C',1,k),repmat('B',1,n));

% Calculate feasible EWM encouragement rule in case 1
[sol, v] = cplexmilp(f_M1,Aineq,bineq,[],[],[],[],[],lb,ub,ctype,[],opt);
lambda_M1 = sol(1:k,:);
v_M1 = -v;
in_Ghat_M1 = (V*lambda_M1>0);
disp('Feasible EWM encouragement rule: alpha=2.5 (median tuition fee)')
fprintf('Proportion treated:\n%.4f\n',mean(in_Ghat_M1));
fprintf('Est. welfare gain:\n%.4f\n',mean(g_M1.*in_Ghat_M1));
fprintf('Average Propensity Score Changes:\n%.4f\n',mean((p_M1-p).*in_Ghat_M1));
fprintf('PRTE:\n%.4f\n',mean(g_M1.*in_Ghat_M1)/mean((p_M1-p).*in_Ghat_M1));
fprintf('Est. total costs:\n%.4f\n',mean((Z1-M1).*p_M1.*in_Ghat_M1));

% Calculate feasible EWM encouragement rule in case 2
[sol, v] = cplexmilp(f_M2,Aineq,bineq,[],[],[],[],[],lb,ub,ctype,[],opt);
lambda_M2 = sol(1:k,:);
v_M2 = -v;
in_Ghat_M2 = (V*lambda_M2>0);
disp('Feasible EWM encouragement rule: alpha=22.25 (maximum tuition fee)')
fprintf('Proportion treated:\n%.4f\n',mean(in_Ghat_M2));
fprintf('Est. welfare gain:\n%.4f\n',mean(g_M2.*in_Ghat_M2));
fprintf('Average Propensity Score Changes:\n%.4f\n',mean((p_M2-p).*in_Ghat_M2));
fprintf('PRTE:\n%.4f\n',mean(g_M2.*in_Ghat_M2)/mean((p_M2-p).*in_Ghat_M2));
fprintf('Est. total costs:\n%.4f\n',mean((Z1-M2).*p_M2.*in_Ghat_M2));

% Calculate budget-constrained EWM encouragement rule in case 1
[sol, v] = cplexmilp(f_M1,Aineq_M1,bineq_bc,[],[],[],[],[],lb,ub,ctype,[],opt);
lambda_M1_bc = sol(1:k,:);
v_M1_bc = -v;
in_Ghat_M1_bc = (V*lambda_M1_bc>0);
disp('Budget-constrained EWM encouragement rule: alpha=2.5 (median tuition fee)')
fprintf('Proportion treated:\n%.4f\n',mean(in_Ghat_M1_bc));
fprintf('Est. welfare gain:\n%.4f\n',mean(g_M1.*in_Ghat_M1_bc));
fprintf('Average Propensity Score Changes:\n%.4f\n',mean((p_M1-p).*in_Ghat_M1_bc));
fprintf('PRTE:\n%.4f\n',mean(g_M1.*in_Ghat_M1_bc)/mean((p_M1-p).*in_Ghat_M1_bc));
fprintf('Est. total costs:\n%.4f\n',mean((Z1-M1).*p_M1.*in_Ghat_M1_bc));

% Calculate budget-constrained EWM encouragement rule in case 2
[sol, v] = cplexmilp(f_M2,Aineq_M2,bineq_bc,[],[],[],[],[],lb,ub,ctype,[],opt);
lambda_M2_bc = sol(1:k,:);
v_M2_bc = -v;
in_Ghat_M2_bc = (V*lambda_M2_bc>0);
disp('Budget-constrained EWM encouragement rule: alpha=22.25 (maximum tuition fee)')
fprintf('Proportion treated:\n%.4f\n',mean(in_Ghat_M2_bc));
fprintf('Est. welfare gain:\n%.4f\n',mean(g_M2.*in_Ghat_M2_bc));
fprintf('Average Propensity Score Changes:\n%.4f\n',mean((p_M2-p).*in_Ghat_M2_bc));
fprintf('PRTE:\n%.4f\n',mean(g_M2.*in_Ghat_M2_bc)/mean((p_M2-p).*in_Ghat_M2_bc));
fprintf('Est. total costs:\n%.4f\n',mean((Z1-M2).*p_M2.*in_Ghat_M2_bc));

% Prepare patch for feasible EWM encouragement rule in case 1
% (Calculate the cutoff of fees at different levels of distance to school)
line_dist_sec = [min(Z2):0.01:max(Z2)]';
line_exp_M1 = Vscale(1,1)*(-(lambda_M1(1) + line_dist_sec*lambda_M1(3)./Vscale(1,2))./lambda_M1(2));

if (lambda_M1(2)>0) % encouragement rule is increasing in fees
    select_M1 = (line_exp_M1 <= 22.25);
    patchM1_X = [line_dist_sec(select_M1); flipud(line_dist_sec(select_M1))];
    patchM1_Y = [max(line_exp_M1(select_M1),0); 22.25*ones(length(line_exp_M1(select_M1)),1)];
else               % encouragement rule is decreasing in fees   
    select_M1 = (line_exp_M1 >= 0);
    patchM1_X = [line_dist_sec(select_M1); flipud(line_dist_sec(select_M1))];
    patchM1_Y = [min(line_exp_M1(select_M1),22.25); zeros(length(line_exp_M1(select_M1)),1)];
end

% Prepare patch for feasible EWM encouragement rule in case 2
line_exp_M2 = Vscale(1,1)*(-(lambda_M2(1) + line_dist_sec*lambda_M2(3)./Vscale(1,2))./lambda_M2(2));

if (lambda_M2(2)>0) % encouragement rule is increasing in fees
    select_M2 = (line_exp_M2 <= 22.25);
    patchM2_X = [line_dist_sec(select_M2); flipud(line_dist_sec(select_M2))];
    patchM2_Y = [max(line_exp_M2(select_M2),0); 22.25*ones(length(line_exp_M2(select_M2)),1)];
else               % encouragement rule is decreasing in fees   
    select_M2 = (line_exp_M2 >= 0);
    patchM2_X = [line_dist_sec(select_M2); flipud(line_dist_sec(select_M2))];
    patchM2_Y = [min(line_exp_M2(select_M2),22.25); zeros(length(line_exp_M2(select_M2)),1)];
end

% Prepare patch for budget-constrained EWM encouragement rule in case 1
line_exp_M1_bc = Vscale(1,1)*(-(lambda_M1_bc(1) + line_dist_sec*lambda_M1_bc(3)./Vscale(1,2))./lambda_M1_bc(2));

if (lambda_M1_bc(2)>0) % encouragement rule is increasing in fees
    select_M1_bc = (line_exp_M1_bc <= 22.25);
    patchM1_X_bc = [line_dist_sec(select_M1_bc); flipud(line_dist_sec(select_M1_bc))];
    patchM1_Y_bc = [max(line_exp_M1_bc(select_M1_bc),0); 22.25*ones(length(line_exp_M1_bc(select_M1_bc)),1)];
else               % encouragement rule is decreasing in fees   
    select_M1_bc = (line_exp_M1_bc >= 0);
    patchM1_X_bc = [line_dist_sec(select_M1_bc); flipud(line_dist_sec(select_M1_bc))];
    patchM1_Y_bc = [min(line_exp_M1_bc(select_M1_bc),22.25); zeros(length(line_exp_M1_bc(select_M1_bc)),1)];
end

% Prepare patch for budget-constrained EWM encouragement rule in case 2
line_exp_M2_bc = Vscale(1,1)*(-(lambda_M2_bc(1) + line_dist_sec*lambda_M2_bc(3)./Vscale(1,2))./lambda_M2_bc(2));

if (lambda_M2_bc(2)>0) % encouragement rule is increasing in fees
    select_M2_bc = (line_exp_M2_bc <= 22.25);
    patchM2_X_bc = [line_dist_sec(select_M2_bc); flipud(line_dist_sec(select_M2_bc))];
    patchM2_Y_bc = [max(line_exp_M2_bc(select_M2_bc),0); 22.25*ones(length(line_exp_M2_bc(select_M2_bc)),1)];
else               % encouragement rule is decreasing in fees   
    select_M2_bc = (line_exp_M2_bc >= 0);
    patchM2_X_bc = [line_dist_sec(select_M2_bc); flipud(line_dist_sec(select_M2_bc))];
    patchM2_Y_bc = [min(line_exp_M2_bc(select_M2_bc),22.25); zeros(length(line_exp_M2_bc(select_M2_bc)),1)];
end

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

h = figure('Color','white');
set(gcf,'Position',[0 0 585 225]);
subplot(1,2,1)
patch(patchM1_X,patchM1_Y,patch2_color,'LineStyle','-','FaceAlpha',Face2Alpha2,'EdgeColor',patch2_edge);
hold on
patch(patchM2_X,patchM2_Y,patch1_color,'LineStyle','-','FaceAlpha',Face2Alpha1,'EdgeColor',patch1_edge);
scatter(dist_sec_u,exp_u,nw*0.4,'MarkerEdgeColor','none','MarkerFaceColor','black');
title('Panel A: Feasible EWM Encouragement Rule')
set(gca,'FontSize',6)
set(gca,'Position',[0.04 0.12 0.44 0.83])
axis([min(Z2)-0.5 max(Z2)+0.5 -1 22.25]);
xlabel('Distance to School (in km)');
ylabel('Fees per Continuing Student (in 1000 Rupiah)');
ax = gca;
ax.YTick = [0,5,10,15,20];
ax.YTickLabel = {'0','5','10','15','20'};
legend('Up to 2500 Rupiah subsidy','Up to 22250 Rupiah subsidy','Population Density','Location','northeast');

subplot(1,2,2)
patch(patchM1_X_bc,patchM1_Y_bc,patch2_color,'LineStyle','-','FaceAlpha',Face2Alpha2,'EdgeColor',patch2_edge);
hold on
patch(patchM2_X_bc,patchM2_Y_bc,patch1_color,'LineStyle','-','FaceAlpha',Face2Alpha1,'EdgeColor',patch1_edge);
scatter(dist_sec_u,exp_u,nw*0.4,'MarkerEdgeColor','none','MarkerFaceColor','black');
title('Panel B: Budget-Constrained EWM Encouragement Rule')
set(gca,'FontSize',6)
set(gca,'Position',[0.53 0.12 0.44 0.83])
axis([min(Z2)-0.5 max(Z2)+0.5 -1 22.25]);
xlabel('Distance to School (in km)');
ylabel('Fees per Continuing Student (in 1000 Rupiah)');
ax = gca;
ax.YTick = [0,5,10,15,20];
ax.YTickLabel = {'0','5','10','15','20'};
legend('Up to 2500 Rupiah subsidy','Up to 22250 Rupiah subsidy','Population Density','Location','northeast');

saveas(h,savefilename,'epsc');