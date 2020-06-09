
clear; clc;

%% global params

n = 500;    % #firms
t = 2000;   % #time horizon (20*months)

% loop params
df_m = 5;   % degree of freedom of marginal t-dist
nu = 3;     % degree of freedom of t-copula
rho = 0;  % pairwise correlation

%% run simulation
[delta,p_value,pMSE] = Simulate(n,t,rho,nu,df_m);
