
clear; clc;

%% global params

n = 500;    % #firms
t = 2000;   % #time horizon (20*months)
df_m = 3;   % degree of freedom of marginal t-dist
nu = 1000;     % degree of freedom of t-copula
rho = 0;  % pairwise correlation

%%
[delta,pMSE,p_KS] = Simulate(rho,nu,df_m);