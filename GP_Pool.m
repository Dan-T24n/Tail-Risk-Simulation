function [k,sigma] = GP_Pool(X)

%Estimate the shape parameter of the tail distribution 
%   X is 20*n dimension data


x = reshape(X,1,[]);
q = quantile(x,.05);
y = q - x(x<q);

paramEsts = gpfit(y);
k      = paramEsts(1);   % "Tail index" parameter - or shape
sigma  = paramEsts(2);    % scale parameter

end

