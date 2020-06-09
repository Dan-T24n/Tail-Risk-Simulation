function [delta,p_value,pMSE] = Simulate(n,t,rho,nu,df_m)

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
GP_k = zeros(1,M);
GP_sigma = zeros(1,M);
pMSE_GP = zeros(1,M);
pMSE_Kelly = zeros(1,M);
pMSE_smooth = zeros(1,M);

%Loop over months
for i = 1 : M
    
%split data into months
idx = i*20-linspace(19,0,20);
X = data (idx,:);

%compute Tal index
Kelly(i) = CSTR(X);
Smooth(i) = SmoothCSTR(X);
[GP_k(i),GP_sigma(i)] = GP_Pool(X);

%fit cdf function
x = reshape(X,1,[]);    
q = quantile(x,.05);    
y = x(x<q)/q;
[F,yi] = ecdf(y);

z = q - x(x<q);
[F1,zi] = ecdf(z);

cdf_Kelly = CDF_Tail(yi,Kelly(i));
cdf_smooth = CDF_Tail(yi,Smooth(i));
cdf_GP = gpcdf(zi,GP_k(i),GP_sigma(i));

%compute MSE
pMSE_GP(i) = Fitness(cdf_GP,F);
pMSE_Kelly(i) = Fitness(cdf_Kelly,F);
pMSE_smooth(i) = Fitness(cdf_smooth,F);

end

%% Compare pMSE & check statistical difference

Z = [pMSE_GP; pMSE_Kelly; pMSE_smooth];

pMSE = mean(Z,2);

delta = pMSE(3) - pMSE(2);

% Wilcoxon ranksum test pair-wise: location test
[P_Wilcoxon,H_Wilcoxon] = Wilcoxon_matrix(Z);


% Kolmogorov-Smirnov 2-sample test pair-wise: distribution test
%[P_KS,H_KS] = KS_matrix(Z);

p_value = P_Wilcoxon(2,3);

end

