
clear; clc;

%% global params 

n = 500;    % #firms
t = 2000;   % #time horizon (20*months)
df_m = 3;   % degree of freedom of marginal t-dist
nu = 5;     % degree of freedom of t-copula
rho = 0.5;  % pairwise correlation

% to check theoretical case, set rho = 0 and nu > 20 (approach Gaussian)
% exponent estimate should approach df_m as sample approach t(df_m)

%% generate data

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
[pMSE_GP(i),pR2_GP(i)] = Fitness(cdf_GP,F);
[pMSE_Kelly(i),pR2_Kelly(i)] = Fitness(cdf_Kelly,F);
[pMSE_smooth(i),pR2_smooth(i)] = Fitness(cdf_smooth,F);

end

%% Compare pMSE & check statistical difference

%pMSE
z = [pMSE_GP; pMSE_Kelly; pMSE_smooth];
month = linspace(1,100,100);

figure();
line(month,z);
title(['Long-run MSE Compare']);
xlabel(['Mean GP = ',num2str(mean(pMSE_GP)),', Mean Kelly = ',num2str(mean(pMSE_Kelly)),', Mean Smooth = ',num2str(mean(pMSE_smooth))]);
legend('Fitted GP','Fitted Kelly','Fitted Smooth','Location','best');
set(gca,'FontSize',15)

% statistical test
Z = [pMSE_GP; pMSE_Kelly; pMSE_smooth];

pMSE = mean(Z,2);

% Wilcoxon ranksum test pair-wise: location test
[P_Wilcoxon,H_Wilcoxon] = Wilcoxon_matrix(Z);

% Kolmogorov-Smirnov 2-sample test pair-wise: distribution test
[P_KS,H_KS] = KS_matrix(Z);

p_KS = P_KS(2,3);

%% Compare pR2

% %pR2
% z = [pR2_GP; pR2_Kelly; pR2_smooth];
% month = linspace(1,100,100);
% 
% figure();
% line(month,z);
% title(['Long-run Pseudo-R^2 Compare']);
% xlabel(['Mean GP = ',num2str(mean(pR2_GP)),', Mean Kelly = ',num2str(mean(pR2_Kelly)),', Mean Smooth = ',num2str(mean(pR2_smooth))]);
% legend('Fitted GP','Fitted Kelly','Fitted Smooth','Location','best');
% set(gca,'FontSize',15)

%% Plotting Kelly vs. Smooth
z = [Kelly; Smooth];
month = linspace(1,100,100);

figure();
line(month,z);
title(['Long-run Tail Index measures']);
xlabel(['Mean Kelly = ',num2str(mean(Kelly)),', Mean Smooth = ',num2str(mean(Smooth))]);
set(gca,'FontSize',15)

exponent_Kelly = 1/mean(Kelly);
exponent_Smooth = 1/mean(Smooth);

%% check month with largest difference in TR estimate - Warning: Black swan ahead 

% find the month and extract data
z = Kelly - Smooth;
[M,I] = max(z);
month = I;
idx = month*20-linspace(19,0,20);

% check daily returns & detect the swan
X = data (idx,:);
X1 = reshape(X',1,[]);

X2 = X.^2;  %square all returns
w = sum(X2,2); %sum by row: find day with largest variance
[M,I] = max(w);
day = round(I);

q = quantile(X(I,:),.05);  
daymin = min(X(I,:));

figure();
plot(X1)
title(['Month ',num2str(month),' may contain black swans']);
xlabel(['The biggest swan is day ',num2str(day)]);
set(gca,'FontSize',15)

figure();
plot(X(I,:))
title(['The biggest (black) swan is day ',num2str(day)]);
xlabel(['Average return = ',num2str(mean(X(I,:))),', Lowest return = ',num2str(daymin),', Threshold return = ',num2str(q)]);
set(gca,'FontSize',15)

% smooth will consistently under-estimate black swans
% avoid market downturns

%% Bring in month_test to check fitness

%compute params
Kelly = CSTR(X);
Smooth = SmoothCSTR(X);
[GP_k,GP_sigma] = GP_Pool(X);

% Fit CDF function
x = reshape(X,1,[]);    %reshape data into 1 vector
q = quantile(x,.05);    %compute q left tail threshold at 5%
y = x(x<q)/q;

[F,yi] = ecdf(y); %empirical cumulative function
cdf_Kelly = CDF_Tail(yi,Kelly);
cdf_Smooth = CDF_Tail(yi,Smooth);

z = q - x(x<q);
[F1,zi] = ecdf(z);
cdf_GP = gpcdf(zi,GP_k,GP_sigma);

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

%compare fitness
pMSE_GP = Fitness(cdf_GP,F);
pMSE_Kelly = Fitness(cdf_Kelly,F);
pMSE_smooth = Fitness(cdf_Smooth,F);



%% check month with smallest difference

% find the month and extract data
z = Kelly - Smooth;
[M,I] = min(z);
month = I;
idx = month*20-linspace(19,0,20);

% check daily returns
X = data (idx,:);
X1 = reshape(X',1,[]);

X2 = X.^2;
w = sum(X2,2);
[M,I] = max(w);
day = round(I);
q = quantile(X(I,:),.05);  
daymax = max(X(I,:));
daymin = min(X(I,:));

figure();
plot(X1)    
title(['Month ',num2str(month),' with lowest difference']);
xlabel(['The highest variance is day ',num2str(day)]);
set(gca,'FontSize',15)

figure();
plot(X(I,:))
title(['The highest variance day ',num2str(day)]);
xlabel(['Average return = ',num2str(mean(X(I,:))),', Lowest return = ',num2str(daymin),', Threshold return = ',num2str(q)]);
set(gca,'FontSize',15)

%% Bring in month_test to check fitness

%compute params
Kelly = CSTR(X);
Smooth = SmoothCSTR(X);
[GP_k,GP_sigma] = GP_Pool(X);

% Fit CDF function
x = reshape(X,1,[]);    %reshape data into 1 vector
q = quantile(x,.05);    %compute q left tail threshold at 5%
y = x(x<q)/q;

[F,yi] = ecdf(y); %empirical cumulative function
cdf_Kelly = CDF_Tail(yi,Kelly);
cdf_Smooth = CDF_Tail(yi,Smooth);

z = q - x(x<q);
[F1,zi] = ecdf(z);
cdf_GP = gpcdf(zi,GP_k,GP_sigma);

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

%compare fitness
pMSE_GP = Fitness(cdf_GP,F);
pMSE_Kelly = Fitness(cdf_Kelly,F);
pMSE_smooth = Fitness(cdf_Smooth,F);
