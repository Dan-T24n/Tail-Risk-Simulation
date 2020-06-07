
clear; clc;

%%  generate normal distribution

rng(5,'twister'); %set random seed
m = 100000;      %number of data points
b = round(sqrt(m));     %number of bins for histogram

x = normrnd(0,1,m,1);  

q = quantile(x,.05);    % 5% quantile
histogram(x,b,'Normalization','probability');
hold on
plot([q q],[0 2/b],'--r')
hold off

%% use exceedances to check tail exponent 

y = q - x(x<q);
figure;
histogram(y,b,'Normalization','probability'); 
%histogram(log(y),100,'Normalization','probability')

%estimate params of standard Generalized Pareto (2 params only, no location)
[paramEsts,paramCI] = gpfit(y);
kHat      = paramEsts(1)   % "Tail index" parameter - or shape
sigmaHat  = paramEsts(2)    % scale parameter
[nll,acov] = gplike(paramEsts, y);
stdErr = sqrt(diag(acov))
kCI  = paramCI(:,1)


%fit PDF function GP distribution
bins = 0:.01:100;
figure;
h = bar(bins,histc(y,bins)/(length(y)*.01),'histc');
h.FaceColor = [.9 .9 .9];
ygrid = linspace(0,1.1*max(y),100);
line(ygrid,gppdf(ygrid,kHat,sigmaHat),'Color','red');
xlim([0,quantile(y,.99)]);
xlabel('Exceedance');
ylabel('Probability Density');

%fit CDF function of GP
[F,yi] = ecdf(y);
figure;
plot(yi,gpcdf(yi,kHat,sigmaHat),'-');
hold on;
stairs(yi,F,'r');
hold off;
legend('Fitted Generalized Pareto CDF','Empirical CDF')


%compare with HILL
ysort = sort(y,'descend');     
T = 20;             
hill1 = 1/mean(log(ysort(1:T)/-q))

ysort = sort(y,'descend');     
T = 20;             
hill2 = 1/mean(log(ysort(1:T)/ysort(T+1)))

ysort = sort(x,'descend');     
T = 20;             
hill3 = 1/mean(log(ysort(1:T)/ysort(T+1))) %change T and test in rolling windows

%% extract prob and curve fitting - in-house estimation

%extract frequency
hh = histogram(y,100,'Normalization','probability');
counts = hh.BinCounts;
freq = counts'/sum(counts);
k = length(counts)/length(y);
pos = (1:k:length(y));
y_bin = y(pos);
