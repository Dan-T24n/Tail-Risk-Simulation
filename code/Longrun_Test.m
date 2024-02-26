
clear; clc;

%% =====EQC MODEL global params=========

n = 1000;    % #firms
t = 2000;   % #time horizon (20*months)
df_m = 3;   % degree of freedom of marginal t-dist
nu = 5;   % degree of freedom of t-copula
rho = 0.2;  % pairwise correlation
df_m = 20;  % degree of freedom of marginal t-dist
nu = 20 ;    % degree of freedom of t-copula
rho = 0.1;  % pairwise correlation

% to check theoretical case, set rho = 0 and nu > 20 (approach Gaussian)
% exponent estimate should approach df_m as sample approach t(df_m)

% generate data

data = genData(n,t,rho,nu,df_m);

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

%% compare TR measures for long-run

[Kelly, Smooth]=ComputeTR(data);

z = [Kelly; Smooth];

%pre-allocate empty vectors to speed up
M = t/20;
Kelly = zeros(1,M);
Smooth = zeros(1,M);
MonthTR = zeros(1,M); 
GP_k = zeros(1,M);
GP_sigma = zeros(1,M);
pMSE_GP = zeros(1,M);
pMSE_Kelly = zeros(1,M);
pMSE_smooth = zeros(1,M);
q = zeros(1,M);


for i = 1 : M
    
%split data into months
idx = i*20-linspace(19,0,20);
X = data (idx,:);

%compute Tal index
Kelly(i) = CSTR(X);
[Smooth(i),Smooth_flag(i)] = SmoothCSTR(X);
MonthTR(i) = TR_MonthCm(X); 
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

flags = find(Smooth_flag);
flag_count = sum(Smooth_flag());

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
>>>>>>> Stashed changes
month = linspace(1,100,100);


figure();
line(month,z);
<<<<<<< Updated upstream
title(['Long-run Tail Risk measures']);
xlabel(['Mean Kelly = ',num2str(mean(Kelly)),', Mean Smooth = ',num2str(mean(Smooth))]);
=======
title(['Long-run Tail Index measures']);
xlabel(['Mean Kelly = ',num2str(mean(Kelly)),', Mean Smooth = ',num2str(mean(Smooth)),', Mean q = ',num2str(mean(q))]);
hold on
legend('Fitted Kelly','Fitted Smooth','Fitted MonthTR','Location','best');
>>>>>>> Stashed changes
set(gca,'FontSize',15)
hold off

%with q ratio (q/mean(q) - 1)
q_ratio = q/mean(q)-1;

figure();
line(month,z);
title(['Long-run Tail Index measures']);
xlabel(['Mean Kelly = ',num2str(mean(Kelly)),', Mean Smooth = ',num2str(mean(Smooth)),', Mean q = ',num2str(mean(q))]);
hold on
plot(q/mean(q)-1,'Color','black');
legend('Fitted Kelly','Fitted Smooth','Fitted MonthTR','q ratio','Location','best');
set(gca,'FontSize',15)
hold off

%with GP
% figure();
% line(month,z);
% title(['Long-run Tail Index measures']);
% xlabel(['Mean Kelly = ',num2str(mean(Kelly)),', Mean Smooth = ',num2str(mean(Smooth)),', Mean q = ',num2str(mean(q))]);
% hold on
% plot(GP_k,'Color','black');
% legend('Fitted Kelly','Fitted Smooth','Fitted MonthTR','GP_k','Location','best');
% set(gca,'FontSize',15)
% hold off


% correlation tests
z = [Kelly; Smooth; MonthTR; q_ratio];
C = corr(z');
C(3,1)
C(3,2)

% rho(MonthTR,q_ratio)<0: stocks with strong swings (high q) may experience
% mean-reversing dynamics => weak CmReturn, while stocks with no swing (low q) have stable
% trajectories => strong CmReturn


%% check month with largest difference - Warning: Black swan ahead 

% find the month and extract data
z = Kelly - Smooth;
[M,I] = max(z)
month = I;
idx = month*20-linspace(19,0,20);

% check daily returns & detect the swan
X = data (idx,:);
X1 = reshape(X',1,[]);

X2 = X.^2;
w = sum(X2,2);
[M,I] = max(w);
day = round(I);
<<<<<<< Updated upstream
daymax = max(X(I,:));
daymin = min(X(I,:));
=======

q = quantile(X(day,:),.05);  
daymin = min(X(day,:));
>>>>>>> Stashed changes

figure();
plot(X1)
title(['Month ',num2str(month),' with highest difference in estimates']);
xlabel(['The biggest swan is day ',num2str(day),', Avg = ',num2str(mean(X1))]);
set(gca,'FontSize',15)

figure();
plot(X(day,:))
title(['The biggest (black) swan is day ',num2str(day)]);
<<<<<<< Updated upstream
xlabel(['Lowest return = ',num2str(daymin)]);
=======
xlabel(['Average return = ',num2str(mean(X(day,:))),', Lowest return = ',num2str(daymin),', Threshold return = ',num2str(q)]);
>>>>>>> Stashed changes
set(gca,'FontSize',15)

% month_test to check fitness
plotFitness(X)

% smooth will consistently under-estimate black swans
% avoid market downturns

%% check month with smallest difference

% find the month and extract data
z = Kelly - Smooth;
[M,I] = min(z)
month = I;
idx = month*20-linspace(19,0,20);

% check daily returns
X = data (idx,:);
X1 = reshape(X',1,[]);

X2 = X.^2;
w = sum(X2,2);
[M,I] = max(w);
day = round(I);
daymax = max(X(I,:));
daymin = min(X(I,:));

figure();
<<<<<<< Updated upstream
plot(X1)    %   Smooth cannot capture black swan :(
title(['Month ',num2str(month),' may contain black swans']);
xlabel(['The biggest swan is day ',num2str(day)]);
=======
plot(X1)    
title(['Month ',num2str(month),' with lowest difference in estimates']);
xlabel(['The highest variance is day ',num2str(day)]);
>>>>>>> Stashed changes
set(gca,'FontSize',15)

figure();
plot(X(I,:))
<<<<<<< Updated upstream
title(['The biggest (black) swan is day ',num2str(day)]);
xlabel(['Lowest return = ',num2str(daymin)]);
set(gca,'FontSize',15)
=======
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
>>>>>>> Stashed changes
