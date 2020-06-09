function [P,H] = KS_matrix(Z)
%Kolmogorov-Smirnov test for the distribution of pMSE
%   Z is 3x100 matrix, each row is a vector of pMSE for 1 measure
%   null = same distribution


k = size(Z,1);

P = zeros(k);
H = zeros(k);


for i = 1 : k

    for j = i : k
        
    [H(i,j),P(i,j)] = kstest2(Z(i,:),Z(j,:));
               
    end
    
end



end