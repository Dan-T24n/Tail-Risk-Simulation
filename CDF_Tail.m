function [cdf] = CDF_Tail(y,kHAT)

%Rebuild the CDF function of the tail exceedances 
%   according to the power law hypothesis of Kelly
%   y is the sorted* vector of relative exceedances (R/q) with R<q
%   kHat is the estimate of the tail index
   
expo = -1/kHAT;

prob = (y).^expo;

cdf = sort(prob);

end

