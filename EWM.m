% 1/20/2024 Yan Liu
% Calculate the feasible EWM encouragement rule and the budget-constrained 
% EWM encouragement rule
% Replicate Figure 1

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
Z12 = Z1.*Z2;
X = [data.ar09 data.ar09.^2 data.rural data.dist_health ...
    data.protestant data.catholic data.religion_other ...
    data.ele_p data.sec_p data.missing_p data.ele_m data.sec_m data.missing_m...
    data.n_sumatra data.w_sumatra data.s_sumatra data.lampung ...
    data.c_java data.yogyakarta data.e_java data.bali ...
    data.w_nussa_tengara data.s_kalimanthan data.s_sulawesi];
XZ1 = X.*Z1;
XZ2 = X.*Z2;
p = data.phat;
Xp = X.*p;
Z2p = Z2.*p;
Vhat = Y-[X Z2 Xp Z2p]*betae;
as = [2.5 22.25]; %tuition subsidy levels
h = 0.06; %bandwidth

% Specify color and transparency for patches
color = true;
Face2Alpha    = [0.3 0.7];
if (color)
    patch_edge  = [0.3 0.4 0.3;0.6 0.4 0.1];
    patch_color = [0.6 0.7 0.6;0.93 0.79 0.57];
    savefilename ='EWMrule_IFLS_color';
else
    patch_edge  = [0.5 0.5 0.5;0.2 0.2 0.2];
    patch_color = [0.8 0.8 0.8;0.5 0.5 0.5];
    savefilename = 'EWMrule_IFLS';
end

% Create regressor matrix V from the data
V  = [Z1 Z2];

% Rescale covariates to [-1,1]
Vscale = ones(n,1)*max(abs(V));
V = [ones(n,1) V./Vscale];
k = size(V,2);

C = sum(abs(V),2); % maximum values of v'beta
minmargin = max(1,C)*(1e-8); % prevent non-integer numbers the integrality 
% constraint of integers from being counted as integers
Aineq = [[-V diag(C)];[V -diag(C)]];
bineq = [[C-minmargin];[-minmargin]];
bineq_bc = [[C-minmargin];[-minmargin];0.28*n];
lb = [-ones(k,1);zeros(n,1)];
ub = [ones(k,1);ones(n,1)];

% Variable type string
ctype = strcat(repmat('C',1,k),repmat('B',1,n));

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

patch_X = [];
patch_Y = [];
patch_X_bc = [];
patch_Y_bc = [];
Npatch = zeros(2,2);

for d = 1:2
    a = as(d);
    
    % Calculate propensity scores after manipulation
    alphaZ = (Z1-a).*(Z1>=a);
    Zalpha = [alphaZ X.*alphaZ Z2 XZ2 alphaZ.*Z2];
    palpha = predictp([ones(n,1) X],Zalpha,gamma);
    
    % Calculate operator kernel in social welfare criterion
    llr = zeros(n,2);
    for i = 1:n
        W = diag(normpdf((p-p(i))/h));
        P = [ones(n,1) p-p(i)];
        theta = (P'*W*P)\(P'*W*Vhat);
        llr(i,1) = theta(1);
    end
    for i = 1:n
        W = diag(normpdf((p-palpha(i))/h));
        P = [ones(n,1) p-palpha(i)];
        theta = (P'*W*P)\(P'*W*Vhat);
        llr(i,2) = theta(1);
    end
    g = [X Z2]*betae(26:50).*(palpha-p)+llr(:,2)-llr(:,1);
    
    % Calculate budget functions
    B = (Z1-alphaZ).*palpha;
    
    % Report results for policies that encourage everyone
    if a == 2.5
        disp('Encourage everyone: alpha=2.5 (median tuition fee)')
    elseif a == 22.25
        disp('Encourage everyone: alpha=22.25 (maximum tuition fee)')
    end
    fprintf('Est. welfare gain:\n%.4f\n',mean(g));
    fprintf('Average Propensity Score Changes:\n%.4f\n',mean(palpha-p));
    fprintf('PRTE:\n%.4f\n',mean(g)/mean(palpha-p));
    fprintf('Est. total costs:\n%.4f\n',mean(B));
    
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
    
    f = [zeros(k,1); -g]; % objective function coefficients
    % Calculate feasible EWM encouragement rule
    [sol,] = cplexmilp(f,Aineq,bineq,[],[],[],[],[],lb,ub,ctype,[],opt);
    lambda = sol(1:k,:);
    in_Ghat = (V*lambda>0);
    if a == 2.5
        disp('Feasible EWM encouragement rule: alpha=2.5 (median tuition fee)')
    elseif a == 22.25    
        disp('Feasible EWM encouragement rule: alpha=22.25 (maximum tuition fee)')
    end
    fprintf('Proportion treated:\n%.4f\n',mean(in_Ghat));
    fprintf('Est. welfare gain:\n%.4f\n',mean(g.*in_Ghat));
    fprintf('Average Propensity Score Changes:\n%.4f\n',mean((palpha-p).*in_Ghat));
    fprintf('PRTE:\n%.4f\n',mean(g.*in_Ghat)/mean((palpha-p).*in_Ghat));
    fprintf('Est. total costs:\n%.4f\n',mean(B.*in_Ghat));
    
    % Prepare patch for feasible EWM encouragement rule
    % (Calculate the cutoff of fees at different levels of distance to school)
    line_dist_sec = [min(Z2):0.01:max(Z2)]';
    line_exp = Vscale(1,1)*(-(lambda(1) + line_dist_sec*lambda(3)./Vscale(1,2))./lambda(2));
    
    if (lambda(2)>0) % encouragement rule is increasing in fees
        select = (line_exp <= 22.25);
        temp1 = [line_dist_sec(select); flipud(line_dist_sec(select))];
        temp2 = [max(line_exp(select),0); 22.25*ones(length(line_exp(select)),1)];
    else               % encouragement rule is decreasing in fees   
        select = (line_exp >= 0);
        temp1 = [line_dist_sec(select); flipud(line_dist_sec(select))];
        temp2 = [min(line_exp(select),22.25); zeros(length(line_exp(select)),1)];
    end
    
    patch_X = [patch_X;temp1];
    patch_Y = [patch_Y;temp2];
    Npatch(d,1) = length(temp1);
    
    Aineq_bc = [[-V diag(C)]; [V -diag(C)]; [zeros(1,k) B']];
    % Calculate budget-constrained EWM encouragement rule
    [sol_bc,] = cplexmilp(f,Aineq_bc,bineq_bc,[],[],[],[],[],lb,ub,ctype,[],opt);
    lambda_bc = sol_bc(1:k,:);
    in_Ghat_bc = (V*lambda_bc>0);
    if a == 2.5
        disp('Budget-constrained EWM encouragement rule: alpha=2.5 (median tuition fee)')
    elseif a == 22.25
        disp('Budget-constrained EWM encouragement rule: alpha=22.25 (maximum tuition fee)')
    end
    fprintf('Proportion treated:\n%.4f\n',mean(in_Ghat_bc));
    fprintf('Est. welfare gain:\n%.4f\n',mean(g.*in_Ghat_bc));
    fprintf('Average Propensity Score Changes:\n%.4f\n',mean((palpha-p).*in_Ghat_bc));
    fprintf('PRTE:\n%.4f\n',mean(g.*in_Ghat_bc)/mean((palpha-p).*in_Ghat_bc));
    fprintf('Est. total costs:\n%.4f\n',mean(B.*in_Ghat_bc));
    
    % Prepare patch for budget-constrained EWM encouragement rule
    line_exp_bc = Vscale(1,1)*(-(lambda_bc(1) + line_dist_sec*lambda_bc(3)./Vscale(1,2))./lambda_bc(2));
    
    if (lambda_bc(2)>0) % encouragement rule is increasing in fees
        select_bc = (line_exp_bc <= 22.25);
        temp3 = [line_dist_sec(select_bc); flipud(line_dist_sec(select_bc))];
        temp4 = [max(line_exp_bc(select_bc),0); 22.25*ones(length(line_exp_bc(select_bc)),1)];
    else               % encouragement rule is decreasing in fees   
        select_bc = (line_exp_bc >= 0);
        temp3 = [line_dist_sec(select_bc); flipud(line_dist_sec(select_bc))];
        temp4 = [min(line_exp_bc(select_bc),22.25); zeros(length(line_exp_bc(select_bc)),1)];
    end
    
    patch_X_bc = [patch_X_bc;temp3];
    patch_Y_bc = [patch_Y_bc;temp4];
    Npatch(d,2) = length(temp3);
end

h = figure('Color','white');
set(gcf,'Position',[0 0 585 225]);
subplot(1,2,1)
patch(patch_X(1:Npatch(1,1)),patch_Y(1:Npatch(1,1)),patch_color(1,:),...
    'LineStyle','-','FaceAlpha',Face2Alpha(1),'EdgeColor',patch_edge(1,:));
hold on
patch(patch_X((Npatch(1,1)+1):(Npatch(1,1)+Npatch(2,1))),...
    patch_Y((Npatch(1,1)+1):(Npatch(1,1)+Npatch(2,1))),patch_color(2,:),...
    'LineStyle','-','FaceAlpha',Face2Alpha(2),'EdgeColor',patch_edge(2,:));
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
patch(patch_X_bc(1:Npatch(1,2)),patch_Y_bc(1:Npatch(1,2)),patch_color(1,:),...
    'LineStyle','-','FaceAlpha',Face2Alpha(1),'EdgeColor',patch_edge(1,:));
hold on
patch(patch_X_bc((Npatch(1,2)+1):(Npatch(1,2)+Npatch(2,2))),...
    patch_Y_bc((Npatch(1,2)+1):(Npatch(1,2)+Npatch(2,2))),patch_color(2,:),...
    'LineStyle','-','FaceAlpha',Face2Alpha(2),'EdgeColor',patch_edge(2,:));
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