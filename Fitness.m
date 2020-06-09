function [pMSE,pR2] = Fitness(cdf,ecdf)
%Probabilistic Mean Square Error
%   cdf is fitted cdf from hypothesis
%   ecdf is the empirical (true) cdf from data

e = (cdf - ecdf).^2;
pMSE = mean(e)^(1/2);
pR2 = 1 - pMSE;

end

