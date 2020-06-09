function [P,H] = Wilcoxon_matrix(Z)
%Rank sum test for the location of the distribution of pMSE
%   Z is 3x100 matrix, each row is a vector of pMSE for 1 measure
%   null = same location (mean/median) of distribution


k = size(Z,1);

P = zeros(k);
H = zeros(k);


for i = 1 : k

    for j = i : k
        
    [P(i,j),H(i,j)] = ranksum(Z(i,:),Z(j,:));
               
    end
    
end



end

