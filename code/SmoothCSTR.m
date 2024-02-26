function [kHat,flag, kVec] = SmoothCSTR(X)
%SmoothCSTR Summary 
%   Smooth cross-sectional Tail Risk
%   compute TR daily then averaging across the month
%   X is 20*n dimension data


days = size(X,1);
temp = zeros(days,1);
kVec = zeros(days,1);
flag = 0;

%compute vectpor of daily TR (Hill)
for i = 1 : days
    Xtemp = X(i,:);

    q = quantile(Xtemp,.05);    %compute q left tail threshold at 5%

    Y = Xtemp(Xtemp<q);      %extract the tail
    
    temp(i) = mean(log(Y/q));   %daily tail 
    
    if isreal(temp(i))      %check if TR is real number (complex if q>0)
    kVec(i) = temp(i);
    end

end

    if ~isreal(temp)
    flag = 1;
    end

kHat = mean(kVec);

end

