
clear; clc;

%% =====EQC MODEL global params=========

n = 500;      % #firms
t = 2000;   % #time horizon (20*months)
df_m = 5;   % degree of freedom of marginal t-dist
nu = 5;     % degree of freedom of t-copula
rho = 0.4;  % pairwise correlation

% generate data

data = genData(n,t,rho,nu,df_m);

%% =====BLOCK_MODEL global params========

n = 500;      % #firms
t = 2000;   % #time horizon (20*months)
df_m = 3;   % degree of freedom of marginal t-dist
nu = 5;     % degree of freedom of t-copula

b = 4; 
rvec = [0.1 0.5 0.3 0.4];  % pairwise correlation

%test block size
remain = rem(n,b);

if not(remain == 0)
    fprintf('ERRORR***change number of blocks! \n')
else
    fprintf('OK***Proceed \n')
end

clear remain

%generate data

data = genData_block(n,b,t,rvec,nu,df_m);

%% one month sampling
%pick 1 random month
month = randi(t/20);
idx = month*20-linspace(19,0,20);

% sample daily returns 
X = data(idx,:);

% plot month check
X1 = reshape(X',1,[]);
figure();
plot(X1)
title(['Month ',num2str(month),' selected at random']);
xlabel(['Each day has 500 data points ',', Avg = ',num2str(mean(X1)),', Std = ',num2str(std(X1))]);
set(gca,'FontSize',15)

% Compute cumul return
Xc = CmReturn(X);
Xc = Xc/std(Xc);
x = reshape(Xc,1,[]);    %reshape data into 1 vector
q = quantile(x,.05);     %compute q left tail threshold at 5%

figure();
plot(Xc)
title(['Month ',num2str(month),' standardized cumulative returns']);
xlabel(['Avg = ',num2str(mean(Xc)), ', Avg = ',num2str(mean(Xc)),', q = ',num2str(q)]);
set(gca,'FontSize',15)

% check extreme cumul returns
[Max,Imax] = max(Xc);
Xmax = X(:,Imax);
[Min,Imin] = min(Xc);
Xmin = X(:,Imin);

%% compute tail risk with cumulative returns


%estimator Hill measure
MonthTR = TR_MonthCm(X); 

%compute Generalized Pareto params
[GP_k,GP_sigma] = GP_Pool(X);

%% Fit CDF function

%fit relative exceedances to power law
x = reshape(Xc,1,[]);    %reshape data into 1 vector
q = quantile(x,.05);     %compute q left tail threshold at 5%
y = x(x<q)/q;

[F,yi] = ecdf(y);   %empirical cumulative function

cdf_MonthTR = CDF_Tail(yi,MonthTR);     %MonthTR cumulative function

%fit exceedances to GP
z = q - x(x<q);

[F1,zi] = ecdf(z);
cdf_GP = gpcdf(zi,GP_k,GP_sigma);
%cdf_GP_Smooth = gpcdf(zi,GPSmooth_k,GPSmooth_sigma);

figure();
stairs(yi,F,'r'); %yi is just sorted values of y
hold on;
plot(yi,cdf_MonthTR,'b-');
plot(yi,cdf_GP,'c-');   %plotting vs. yi to have same axis (cdf does not changes)
%plot(yi,cdf_GP_Smooth,'y-');  GP_Smooth is way off
hold off;
title(['Fitted cumulative functions']);
legend('Empirical CDF','Fitted MonthTR CDF','Fitted GP CDF','Location','best');
set(gca,'FontSize',15)


%Hill perform better when number of obs is very low.

%% compare fitness - pseudo R2


[pMSE_GP] = Fitness(cdf_GP,F)

[pMSE_MonthTR] = Fitness(cdf_MonthTR,F)



