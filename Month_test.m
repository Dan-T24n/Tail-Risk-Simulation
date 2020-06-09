
clear; clc;

%% global params

n = 500;      % #firms
t = 2000;   % #time horizon (20*months)
df_m = 3;   % degree of freedom of marginal t-dist
nu = 5;     % degree of freedom of t-copula
rho = 0.5;  % pairwise correlation

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

% sample daily returns and plot check
X = data (idx,:);
X1 = reshape(X',1,[]);
figure();
plot(X1)
title(['Month ',num2str(month),' selected at random']);
xlabel(['Each day has 500 data points ']);
set(gca,'FontSize',15)


%% compute tail risk

%compute Kelly measure
Kelly = CSTR(X);

%compute Smooth measure
Smooth = SmoothCSTR(X);

%compute Generalized Pareto params
[GP_k,GP_sigma] = GP_Pool(X);

%compute GP Smooth params
%[GPSmooth_k,GPSmooth_sigma] = GP_Smooth(X);


%% Fit CDF function

%fit relative exceedances to power law
x = reshape(X,1,[]);    %reshape data into 1 vector
q = quantile(x,.05);    %compute q left tail threshold at 5%
y = x(x<q)/q;

[F,yi] = ecdf(y); %empirical cumulative function
cdf_Kelly = CDF_Tail(yi,Kelly);
cdf_Smooth = CDF_Tail(yi,Smooth);


%fit exceedances to GP
z = q - x(x<q);

[F1,zi] = ecdf(z);
cdf_GP = gpcdf(zi,GP_k,GP_sigma);
%cdf_GP_Smooth = gpcdf(zi,GPSmooth_k,GPSmooth_sigma);

figure();
stairs(yi,F,'r'); %yi is just sorted values of y
hold on;
plot(yi,cdf_Kelly,'b-');
plot(yi,cdf_Smooth,'k-');
plot(yi,cdf_GP,'c-');   %plotting vs. yi to have same axis (cdf does not changes)
%plot(yi,cdf_GP_Smooth,'y-');  GP_Smooth is way off
hold off;
title(['Fitted cumulative functions']);
legend('Empirical CDF','Fitted Kelly CDF','Fitted Smooth CDF','Fitted GP CDF','Location','best');
set(gca,'FontSize',15)


%% GP_pool fit test separetly

x = reshape(X,1,[]);
q = quantile(x,.05);
y = q - x(x<q);

[paramEsts,paramCI] = gpfit(y);
GP_kHat      = paramEsts(1);   % "Tail index" parameter - or shape
sigmaHat  = paramEsts(2);    % scale parameter

%fit CDF function of GP
[F,yi] = ecdf(y);
figure();
plot(yi,gpcdf(yi,GP_kHat,sigmaHat),'-');
hold on;
stairs(yi,F,'r');
hold off;
legend('Fitted Generalized Pareto CDF','Empirical CDF','Location','best')
set(gca,'FontSize',15)


%% compare fitness - pseudo R2


[pMSE_GP,pR2_GP] = Fitness(cdf_GP,F);

[pMSE_Kelly,pR2_Kelly] = Fitness(cdf_Kelly,F);

[pMSE_smooth,pR2_smooth] = Fitness(cdf_Smooth,F);


