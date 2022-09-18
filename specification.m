% 09/18/2022 Yan Liu
% Compare polynomial and local quadratic specifications for MTE

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

X1 = [ones(n,1) X].*p;
X0 = [ones(n,1) X].*(1-p);
Z21 = Z2.*p;
Z20 = Z2.*(1-p);

% 2-th order polynomial
W2 = [X0 Z20 X1 Z21 p.^2-p];
theta2 = (W2'*W2)\(W2'*Y);
beta2 = theta2(24:46)-theta2(1:23);
alpha2 = theta2(47);

% 3-th order polynomial
W3 = [X0 Z20 X1 Z21 p.^2-p p.^3-p];
theta3 = (W3'*W3)\(W3'*Y);
beta3 = theta3(24:46)-theta3(1:23);
alpha3 = theta3(47:48);

% 4-th order polynomial
W4 = [X0 Z20 X1 Z21 p.^2-p p.^3-p p.^4-p];
theta4 = (W4'*W4)\(W4'*Y);
beta4 = theta4(24:46)-theta3(1:23);
alpha4 = theta4(47:49);

% Local quadratic regression following Heckman, Urzua, and Vytlacil (2006)
h = 0.27; %bandwidth
% Steps 1&3: fit a local linear regression of Y and each regressor on p
dX = size(X,2);
Xhat = zeros(n,dX+1);
for i = 1:n
    W = diag(normpdf((p-p(i))/h));
    P = [ones(n,1) p-p(i)];
    for j = 1:dX
        theta = (P'*W*P)\(P'*W*X(:,j));
        Xhat(i,j) = theta(1);
    end
    theta = (P'*W*P)\(P'*W*Z2);
    Xhat(i,dX+1) = theta(1);
end
Xphat = zeros(n,dX+1);
Xp = X.*p;
Z2p = Z2.*p;
for i = 1:n
    W = diag(normpdf((p-p(i))/h));
    P = [ones(n,1) p-p(i)];
    for j = 1:dX
        theta = (P'*W*P)\(P'*W*Xp(:,j));
        Xphat(i,j) = theta(1);
    end
    theta = (P'*W*P)\(P'*W*Z2p);
    Xphat(i,dX+1) = theta(1);
end
Yhat = zeros(n,1);
for i = 1:n
    W = diag(normpdf((p-p(i))/h));
    P = [ones(n,1) p-p(i)];
    theta = (P'*W*P)\(P'*W*Y);
    Yhat(i) = theta(1);
end
% Step 2&4: generate residual for Y and each regressor
eY = Y-Yhat;
eX = [X Z2]-Xhat;
eXp = [Xp Z2p]-Xphat;
eW = [eX eXp];
% Step 5: regress e_Y on e_X
betae = (eW'*eW)\(eW'*eY);
Ytilde = Y-eW*betae;
% Step 6: fit a local quadratic regression of Ytilde on p
u = linspace(min(p),max(p),100).';
lqr = zeros(100,1);
for i = 1:100
   W = diag(normpdf((p-u(i))/h));
   P = [ones(n,1) p-u(i) (p-u(i)).^2];
   theta = (P'*W*P)\(P'*W*Ytilde);
   lqr(i) = theta(2);
end

% Calculate estimated MTE evaluated at mean values of X and Z2
Xbar = mean(X).*ones(100,1);
Z2bar = mean(Z2)*ones(100,1);
u2 = 2*u-1;
MTE2 = [ones(100,1) Xbar Z2bar]*beta2+u2*alpha2;
u3 = [2*u 3*u.^2]-1;
MTE3 = [ones(100,1) Xbar Z2bar]*beta3+u3*alpha3;
u4 = [2*u 3*u.^2 4*u.^3]-1;
MTE4 = [ones(100,1) Xbar Z2bar]*beta4+u4*alpha4;
beta = betae(23:44)-betae(1:22);
MTElqr = [Xbar Z2bar]*beta+lqr;

h = figure('Color','white');
plot(u,MTE2,'-','LineWidth',1.8)
hold on
plot(u,MTE3,'--','LineWidth',1.8)
plot(u,MTE4,':','LineWidth',1.8)
plot(u,MTElqr,'-.','LineWidth',1.8)
legend('2-th order','3-th order','4-th order','local quadratic','Location','northeast')
xlabel('$u_1$','interpreter','latex');
ylabel('$\widehat{MTE}_1(u_1,X,Z_2)$','interpreter','latex');
saveas(h,'MTEspecification','epsc');