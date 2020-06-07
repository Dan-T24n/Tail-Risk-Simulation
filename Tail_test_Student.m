
clear; clc;

%%  generate t-distribution

%rng(2,'twister');   %set random seed
m = 10000;      %number of data points
b = round(sqrt(m));     %number of bins for histogram

df = 2;     %degree of freedom
x = trnd(df,m,1);  

q = quantile(x,.05);    % 5% quantile
histogram(x,b,'Normalization','probability')
hold on
plot([q q],[0 4/b],'--r')
hold off


%% use exceedances to check tail exponent 

y = q - x(x<q);
histogram(y,b,'Normalization','probability'); 
%histogram(log(y),100,'Normalization','probability')

%estimate params of standard Generalized Pareto (2 params only, no location)
[paramEsts,paramCI] = gpfit(y);
kHat      = paramEsts(1)   % "Tail index" parameter - or shape
sigmaHat  = paramEsts(2);    % scale parameter
[nll,acov] = gplike(paramEsts, y);
stdErr = sqrt(diag(acov));
kCI  = paramCI(:,1)

%fit CDF function of GP
[F,yi] = ecdf(y);
plot(yi,gpcdf(yi,kHat,sigmaHat),'-');
hold on;
stairs(yi,F,'r');
hold off;
legend('Fitted Generalized Pareto CDF','Empirical CDF')


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


%% compare with HILL

%cannot use Hill on exceedances, completely off
ysort = sort(y,'descend');     
T = 20;             
hill1 = 1/mean(log(ysort(1:T)/-q));

ysort = sort(y,'descend');     
T = 20;             
hill2 = 1/mean(log(ysort(1:T)/ysort(T+1)));

ysort = sort(x,'descend');     
T = 20;             
hill3 = 1/mean(log(ysort(1:T)/ysort(T+1))) %change T and test in rolling windows


% find the relationship f(Hill) = Tail index(GP fit)

% draw a line of g(R-u)=(R-u/u)^exponent


%% build Kelly measure
y = x(x<q);
z=sort(log(y/q));

kelly = mean(log(y/q));

Tail_kelly = 1/kelly %Kelly = Hill3 applied on all values instead T-order


