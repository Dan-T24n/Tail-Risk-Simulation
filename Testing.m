
clear; clc;

%% global params

n = 500;      % #firms
t = 2000;   % #time horizon (20*months)
df_m = 3;   % degree of freedom of marginal t-dist
nu = 5;     % degree of freedom of t-copula
rho = 0.1;  % pairwise correlation

% rho_block *todo*  

%% generate data

data = genData(n,t,rho,nu,df_m);

%% test month sampling - compute tail risk

%pick 1 random month
month = randi(t/20);
idx = month*20-linspace(19,0,20);

%sample daily returns
X = data (idx,:);

%% compute Kelly measure
X = reshape(X,1,[]);    %reshape data into 1 vector

q = quantile(X,.05);    %compute q left tail threshold at 5%

b = round(sqrt(t));
histogram(X,b,'Normalization','probability')
hold on
plot([q q],[0 2/b],'--r')
hold off

%Y = q - X(X<q);
%histogram(Y,b,'Normalization','probability');

Y = X(X<q); 
histogram(Y,b,'Normalization','probability');

Kelly = mean(log(Y/q));

%% compute smooth 
X = data (idx,:);
days = size(X,1);
temp = zeros(days,1);


for i = 1 : days
    Xtemp = X(i,:);

    q = quantile(Xtemp,.05);    %compute q left tail threshold at 5%

    Y = Xtemp(Xtemp<q);      %extract the tail

    temp(i) = mean(log(Y/q));
end

kVec = temp;    %careful complex number!

Smooth = mean(kVec);

day = 7;
x = X(day,:);
q = quantile(x,.05);
Y = x(x<q);
histogram(Y,b,'Normalization','probability');


%% pool_GP fit  test

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

%fit CDF function of GP
[F,yi] = ecdf(y);
plot(yi,gpcdf(yi,kHatGP,sigmaHat),'-');
hold on;
stairs(yi,F,'r');
hold off;
legend('Fitted Generalized Pareto CDF','Empirical CDF')

%% averaging GP

X = data (idx,:);
days = size(X,1);
temp = zeros(days,1);


for i = 1 : days
    Xtemp = X(i,:);

    q = quantile(Xtemp,.05);    %compute q left tail threshold at 5%

    y = q - Xtemp(Xtemp<q);      %extract the tail

    paramEsts = gpfit(y);
    
    temp(i) = paramEsts(1);
end

GP_kVec = temp;
GP_kHat_smooth = mean(temp); 

%% extract prob and curve fitting - in-house estimation

plot(ysort);

%extract frequency
hh = histogram(y,100,'Normalization','probability');
counts = hh.BinCounts;
freq = counts'/sum(counts);
T = length(counts)/length(y);
pos = (1:T:length(y));
y_bin = y(pos);
