
clear; clc;

%% global params

n = 500;    % #firms
t = 2000;   % #time horizon (20*months)
df_m = 3;   % degree of freedom of marginal t-dist
nu = 5;     % degree of freedom of t-copula
rho = 0.1;  % pairwise correlation

%% generate data

data = genData(n,t,rho,nu,df_m);

%% test month sampling

% pick 1 random month
month = randi(t/20);
idx = month*20-linspace(19,0,20);

% sample daily returns and plot check
X = data (idx,:);
X1 = reshape(X',1,[]);
figure();
plot(X1)
title(['Month ',num2str(month),' selected at random']);
xlabel(['Each day has 500 data points ']);
set(gca,'FontSize',15)

%% compute Kelly measure
X = data(idx,:);

Kelly = CSTR(X);

[F,yi] = ecdf(y);


%plot(yi,gpcdf(yi,GP_kHat,sigmaHat),'-');


%% compute smooth 
X = data(idx,:);

Smooth = SmoothCSTR(X);


%% GP_poolfit  test

x = reshape(X,1,[]);
q = quantile(x,.05);
y = q - x(x<q);
histogram(y,b,'Normalization','probability'); 

[paramEsts,paramCI] = gpfit(y);
GP_kHat      = paramEsts(1)   % "Tail index" parameter - or shape
sigmaHat  = paramEsts(2);    % scale parameter
%[nll,acov] = gplike(paramEsts, y);
%stdErr = sqrt(diag(acov));
kGP_CI  = paramCI(:,1);

%fit CDF function of Pool_GP
[F,yi] = ecdf(y);
plot(yi,gpcdf(yi,GP_kHat,sigmaHat),'-');
hold on;
stairs(yi,F,'r');
hold off;
legend('Fitted Generalized Pareto CDF','Empirical CDF')

%% averaging GP


days = size(X,1);
temp1 = zeros(days,1);
temp2 = zeros(days,1);


for i = 1 : days
    Xtemp = X(i,:);

    q = quantile(Xtemp,.05);    %compute q left tail threshold at 5%

    y = q - Xtemp(Xtemp<q);      %extract the tail

    paramEsts = gpfit(y);
    
    temp1(i) = paramEsts(1);
    temp2(i) = paramEsts(2);
end

GP_kVec = temp1;
GP_kHat_smooth = mean(GP_kVec);
GP_sigma_smooth = mean(temp2);

%fit cdf of Smooth_GP
[F,yi] = ecdf(y);
plot(yi,gpcdf(yi,GP_kHat_smooth,GP_sigma_smooth),'-');
hold on;
stairs(yi,F,'r');
hold off;
legend('Fitted Generalized Pareto CDF','Empirical CDF')

