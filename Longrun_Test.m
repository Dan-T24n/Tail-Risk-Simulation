
clear; clc;

%% global params 

n = 500;    % #firms
t = 2000;   % #time horizon (20*months)
df_m = 3;   % degree of freedom of marginal t-dist
nu = 5;   % degree of freedom of t-copula
rho = 0.2;  % pairwise correlation

% to check theoretical case, set rho = 0 and nu > 20 (approach Gaussian)
% exponent estimate should approach df_m as sample approach t(df_m)

%% generate data

data = genData(n,t,rho,nu,df_m);


%% compare TR measures for long-run

[Kelly, Smooth]=ComputeTR(data);

z = [Kelly; Smooth];
month = linspace(1,100,100);

figure();
line(month,z);
title(['Long-run Tail Risk measures']);
xlabel(['Mean Kelly = ',num2str(mean(Kelly)),', Mean Smooth = ',num2str(mean(Smooth))]);
set(gca,'FontSize',15)

exponent_Kelly = 1/mean(Kelly);
exponent_Smooth = 1/mean(Smooth);

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
daymax = max(X(I,:));
daymin = min(X(I,:));

figure();
plot(X1)
title(['Month ',num2str(month),' may contain black swans']);
xlabel(['The biggest swan is day ',num2str(day)]);
set(gca,'FontSize',15)

figure();
plot(X(I,:))
title(['The biggest (black) swan is day ',num2str(day)]);
xlabel(['Lowest return = ',num2str(daymin)]);
set(gca,'FontSize',15)

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
plot(X1)    %   Smooth cannot capture black swan :(
title(['Month ',num2str(month),' may contain black swans']);
xlabel(['The biggest swan is day ',num2str(day)]);
set(gca,'FontSize',15)

figure();
plot(X(I,:))
title(['The biggest (black) swan is day ',num2str(day)]);
xlabel(['Lowest return = ',num2str(daymin)]);
set(gca,'FontSize',15)
