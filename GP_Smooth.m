function [k,sigma,kVec] = GP_Smooth(X)
%Estimate the shape parameter of the tail distribution
%   estimate daily tail riskthen average across the month
%   X is 20*n dimension monthly data

days = size(X,1);
temp1 = zeros(days,1);
temp2 = zeros(days,1);


for i = 1 : days
    
    Xtemp = X(i,:);
    q = quantile(Xtemp,.05);    %compute q left tail threshold at 5%
    y = q - Xtemp(Xtemp<q);      %extract the tail

    paramEsts = gpfit(y);
    
    temp1(i) = paramEsts(1);
    temp2(i) = paramEsts(2);
end

kVec = temp1;
k = mean(temp1);
sigma = mean(temp2);

end

