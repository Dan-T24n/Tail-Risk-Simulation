function [kHat,qMonth,flag] = TR_MonthCm(X)

%CSTR Summary 
%   Month cumulative cross-sectional Tail Risk
%   X is 20*n dimension data

kHat = 0;
flag = 0;

Xc = CmReturn(X);
Xc = Xc/std(Xc);    %standardize, otherwise uncomparable volatility

X = Xc;    %reshape data into 1 vector

q = quantile(X,.05);    %compute q left tail threshold at 5%

y = X(X<q);      %extract the tail

temp = mean(log(y/q));

    if isreal(temp)      %check if TR is real number (complex if q>0)
    kHat = temp;
    end
    
    if ~isreal(temp)
    flag = 1;
    end

qMonth=q;

end

