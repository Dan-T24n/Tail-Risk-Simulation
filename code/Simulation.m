
clear; clc;

%% global params

n = 500;    % #firms
t = 2000;   % #time horizon (20*months)
df_m = 5;   % degree of freedom of marginal t-dist
nu = 1000;     % degree of freedom of t-copula
rho = 0;  % pairwise correlation

% rho_block *todo*  

%% generate data

data = genData(n,t,rho,nu,df_m);

%% test 1 day cross-section
day = randi(t);
X = data (day,:);
Var = 1/length(X)*((X-mean(X))*(X-mean(X))');
plot(X)
title(['cross-section return of day ',num2str(day)]);
xlabel(['Mean = ',num2str(mean(X)),', Variance = ',num2str(Var)]);
set(gca,'FontSize',15)

%% test month sampling - compute tail risk

%pick 1 random month
month = randi(t/20);
idx = month*20-linspace(19,0,20);

%sample daily returns
X = data (idx,:);

%compute Kelly measure
Kelly = CSTR(X);

%compute Smooth measure
[Smooth, kVec] = SmoothCSTR(X);

%% compare TR for long-run

[Kelly, Smooth]=ComputeTR(data);

z = [Kelly; Smooth];
month = linspace(1,100,100);
line(month,z)

mean(Kelly)
mean(Smooth)

kHat_Kelly = 1/mean(Kelly)
kHat_Smooth = 1/mean(Smooth)


%%

% function GPareto

% function Kelly

% function Smooth