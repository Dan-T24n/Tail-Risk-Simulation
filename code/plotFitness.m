function [ ] = plotFitness(X)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% compute params
Kelly = CSTR(X);
Smooth = SmoothCSTR(X);
[GP_k,GP_sigma] = GP_Pool(X);

%% Fit CDF function
x = reshape(X,1,[]);    %reshape data into 1 vector
q = quantile(x,.05);    %compute q left tail threshold at 5%
y = x(x<q)/q;

[F,yi] = ecdf(y); %empirical cumulative function

cdf_Kelly = CDF_Tail(yi,Kelly); %reconstruct cdf from TR formula 
cdf_Smooth = CDF_Tail(yi,Smooth);

z = q - x(x<q);
[F1,zi] = ecdf(z);
cdf_GP = gpcdf(zi,GP_k,GP_sigma);

%% Plot fitted curves

figure();
stairs(yi,F,'r'); %yi is just sorted values of y
hold on;
plot(yi,cdf_Kelly,'b-');
plot(yi,cdf_Smooth,'k-');
plot(yi,cdf_GP,'c-');   %plotting vs. yi to have same axis (cdf does not changes)
hold off;
title(['Fitted cumulative functions']);
legend('Empirical CDF','Fitted Kelly CDF','Fitted Smooth CDF','Fitted GP CDF','Location','best');
set(gca,'FontSize',15)

end

