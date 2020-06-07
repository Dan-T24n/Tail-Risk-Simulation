function [kHat, kVec] = SmoothCSTR(X)
%SmoothCSTR Summary 
%   Smooth cross-sectional Tail Risk
%   compute TR daily then averaging across the month
%   X is 20*n dimension data


days = size(X,1);
temp = zeros(days,1);


for i = 1 : days
    Xtemp = X(i,:);

    q = quantile(Xtemp,.05);    %compute q left tail threshold at 5%

    Y = Xtemp(Xtemp<q);      %extract the tail

    temp(i) = mean(log(Y/q));   %daily tail 
end

kVec = temp;

kHat = mean(kVec);

end

