function [kHat] = CSTR(X)
%CSTR Summary 
%   Kelly cross-sectional Tail Risk
%   X is 20*n dimension data


X = reshape(X,1,[]);    %reshape data into 1 vector

q = quantile(X,.05);    %compute q left tail threshold at 5%

Y = X(X<q);      %extract the tail

kHat = mean(log(Y/q));

end

