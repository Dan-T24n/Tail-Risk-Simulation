function Xc = CmReturn(X)
%Compute cumulative return in 1 month 
%   X is 20*n dimension data


days = size(X,1);
n = size(X,2);
cumul = ones(1,n);


for i = 1 : days
    
    cumul = cumul .* (1 + X(i,:)./100); 
     
end

Xc = (cumul-1)*100;


end

