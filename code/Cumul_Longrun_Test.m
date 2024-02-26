
clear; clc;

%% =====EQC MODEL global params=========

n = 500;    % #firms
t = 2000;   % #time horizon (20*months)
df_m = 5;  % degree of freedom of marginal t-dist
nu = 5 ;    % degree of freedom of t-copula
rho = 0.7;  % pairwise correlation

% to check theoretical case, set rho = 0 and nu > 20 (approach Gaussian)
% exponent estimate should approach df_m as sample approach t(df_m)

% generate data

data = genData(n,t,rho,nu,df_m);

% add negative trend for crisis
 data = data - 2;

%% =====BLOCK_MODEL global params========

n = 500;    % #firms
t = 2000;   % #time horizon (20*months)
df_m = 5;   % degree of freedom of marginal t-dist
nu = 4;     % degree of freedom of t-copula

b = 10;      %# blocks
rvec = [0.5 0.5 0.1 0.1 0.1 0.1 0.5 0.1 0.1 0.1 ];  % pairwise correlation

%test block size
remain = rem(n,b);

if not(remain == 0)
    fprintf('***ERRORR***change number of blocks! \n')
else
    fprintf('***OK***Proceed \n')
end

clear remain

%generate data

data = genData_block(n,b,t,rvec,nu,df_m);

%% Compute Tail index measures for long run

t = size(data,1);
M = t/20;

%pre-allocate empty vectors to speed up
Kelly = zeros(1,M);
Smooth = zeros(1,M);
MonthTR = zeros(1,M);
qMonth = zeros(1,M); 
GP_k = zeros(1,M);
GP_sigma = zeros(1,M);
pMSE_GP = zeros(1,M);
pMSE_Kelly = zeros(1,M);
pMSE_smooth = zeros(1,M);
q = zeros(1,M);


for i = 1 : M
    
%split data into months
idx = i*20-linspace(19,0,20);
X = data(idx,:);

%compute Tal index
Kelly(i) = CSTR(X);
[Smooth(i),Smooth_flag(i)] = SmoothCSTR(X);
[MonthTR(i), qMonth(i), flagMonth(i)] = TR_MonthCm(X); 
[GP_k(i),GP_sigma(i)] = GP_Pool(X);

%fit cdf function
x = reshape(X,1,[]);    


q(i) = quantile(x,.05);    
y = x(x<q(i))/q(i);
[F,yi] = ecdf(y);   %build empirical cdf from ratio excedances (Hill)

z = q(i) - x(x<q(i));
[F1,zi] = ecdf(z);      %build empirical cdf from value excedances (GP)


cdf_Kelly = CDF_Tail(yi,Kelly(i));
cdf_smooth = CDF_Tail(yi,Smooth(i));
cdf_GP = gpcdf(zi,GP_k(i),GP_sigma(i));

%compute MSE
[pMSE_GP(i)] = Fitness(cdf_GP,F);
[pMSE_Kelly(i)] = Fitness(cdf_Kelly,F);
[pMSE_smooth(i)] = Fitness(cdf_smooth,F);

end

flagsDaily = find(Smooth_flag);
flagsMonthly = find(flagMonth);



%% Compare pMSE & check statistical difference

%pMSE
pMSE = [pMSE_Kelly; pMSE_smooth];
month = linspace(1,100,100);

figure();
line(month,pMSE);
title(['Long-run MSE Compare']);
xlabel(['Mean Kelly = ',num2str(mean(pMSE_Kelly)),', Mean Smooth = ',num2str(mean(pMSE_smooth))]);
legend('Fitted Kelly','Fitted Smooth','Location','best');
set(gca,'FontSize',15)

% statistical test
Z = [pMSE_GP; pMSE_Kelly; pMSE_smooth];

% Wilcoxon ranksum test pair-wise: location test
[P_Wilcoxon,H_Wilcoxon] = Wilcoxon_matrix(Z);

% Kolmogorov-Smirnov 2-sample test pair-wise: distribution test
%[P_KS,H_KS] = KS_matrix(Z);

p_value = P_Wilcoxon(2,3)


%% Plotting Kelly vs. Smooth vs. MonthTR
z = [Kelly; Smooth; MonthTR];
month = linspace(1,100,100);


% figure();
line(month,z);
title(['Long-run Tail Index measures']);
xlabel(['Mean Kelly = ',num2str(mean(Kelly)),', Mean Smooth = ',num2str(mean(Smooth)),', Mean q = ',num2str(mean(q))]);
hold on
legend('Fitted Kelly','Fitted Smooth','Fitted MonthTR','Location','best');
set(gca,'FontSize',15)
hold off

%with q ratio (q/mean(q) - 1)
% figure();
% line(month,z);
% title(['Long-run Tail Index measures']);
% xlabel(['Mean Kelly = ',num2str(mean(Kelly)),', Mean Smooth = ',num2str(mean(Smooth)),', Mean q = ',num2str(mean(q))]);
% hold on
% plot(q/mean(q)-1,'Color','black');
% legend('Fitted Kelly','Fitted Smooth','Fitted MonthTR','q ratio','Location','best');
% set(gca,'FontSize',15)
% hold off

%with GP
% figure();
% line(month,z);
% title(['Long-run Tail Index measures']);
% xlabel(['Mean Kelly = ',num2str(mean(Kelly)),', Mean Smooth = ',num2str(mean(Smooth)),', Mean q = ',num2str(mean(q))]);
% hold on
% plot(GP_k,'Color','black');
% legend('Fitted Kelly','Fitted Smooth','Fitted MonthTR','GP','Location','best');
% set(gca,'FontSize',15)
% hold off

%problem with GP: does not converge when rho is high, broken for rho >0.5
%need to use q ratio instead

% correlation tests
q_ratio = q/mean(q)-1;
z = [Kelly; Smooth; MonthTR;q_ratio];
corr_test = corr(z')

% rho(MonthTR,q_ratio)<0: stocks with strong swings (high q) may experience
% mean-reversing dynamics => weak CmReturn, while stocks with no swing (low q) have stable
% trajectories => strong CmReturn


%% check month with largest difference: MonthTR - Smooth 

% find the month and extract data
z = MonthTR - Smooth ;
[M,I] = max(z);
month = I;
idx = month*20-linspace(19,0,20);

% check daily returns & detect the swan
X = data (idx,:);
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
title(['Month ',num2str(month),' with lowest difference in estimates']);
xlabel(['The highest variance is day ',num2str(day)]);
set(gca,'FontSize',15)

figure();
plot(X(I,:))
title(['The highest variance day ',num2str(day)]);
xlabel(['Average return = ',num2str(mean(X(I,:))),', Lowest return = ',num2str(daymin),', Threshold return = ',num2str(q)]);
set(gca,'FontSize',15)

% month_test to check fitness
plotFitness(X)

%% when do smooth outperform? 

% find the month and extract data
z = (pMSE_smooth - pMSE_Kelly)';
[Min_z,I] = min(z);
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
title(['Month ',num2str(month),' with lowest -delta in pMSE']);
xlabel(['The highest variance is day ',num2str(day)]);
set(gca,'FontSize',15)

figure();
plot(X(I,:))
title(['The highest variance day ',num2str(day)]);
xlabel(['Average return = ',num2str(mean(X(I,:))),', Lowest return = ',num2str(daymin),', Threshold return = ',num2str(q)]);
set(gca,'FontSize',15)

% month_test to check fitness
plotFitness(X)

%% Smooth captures cluster of small swans, Kelly capture single swan!