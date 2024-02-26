function [Kelly_Month,Smooth_Month,Kelly_q,Smooth_q] = Cumul_Simulate(n,t,rho,nu,df_m)

% Simulate 1 full scenario and output: 
% delta = mean(pMSE_Smooth) - mean(pMSE_Kelly)
% p-value = p-value of Wilcoxon test for delta
% pMSE = vector of avg pMSE

%   ==params==
% n         % #firms
% t         % #time horizon (20*months)
% rho       % pairwise correlation 
% df_m      % degree of freedom of marginal t-dist
% nu        % degree of freedom of t-copula

%% Generate data

data = genData(n,t,rho,nu,df_m);

%% Compute Tail index measures for long-run

t = size(data,1);
M = t/20;

%pre-allocate empty vectors to speed up
Kelly = zeros(1,M);
Smooth = zeros(1,M);
MonthTR = zeros(1,M); 
q = zeros(1,M);


%Loop over months
for i = 1 : M
    
%split data into months
idx = i*20-linspace(19,0,20);
X = data (idx,:);

%compute Tal index
Kelly(i) = CSTR(X);
Smooth(i) = SmoothCSTR(X);
MonthTR(i) = TR_MonthCm(X);

%compute q
% x = reshape(X,1,[]); 
% q(i) = quantile(x,.05); 

end

%% Compare correlation

% q_ratio = q/(-1.645)-1;

% correlation tests
z = [Kelly; Smooth; MonthTR];
C = corr(z');

Kelly_Month = C(3,1);
Smooth_Month = C(3,2);
% Kelly_q = C(4,1);
% Smooth_q = C(4,2);

end

