function [M] = block_mat(block,size,r)
%Create block matrix for correlations
%   b = nb of blocks
%   s = block size
%   r = vector of rhos, has sizes=(1,b)

%% params

dim = block * size;

M = zeros(dim,dim); 


%% make block matrix

%upper triangle blocks
for b = 1:block
    rho = r(b);
    
    for i = 1:dim
 
       for j = 1:dim
            
            if (j > i)&&(j<=b*size)&&(j>(b-1)*size)&&(i<=b*size)&&(i>(b-1)*size)
                
               M(i,j) = rho;
            else
                if (j == i)
                    M(i,j) = 1;
                   
                end
            end
        end
    end  
end

% reflection over diag

for i = 1:dim
    for j = 1:dim
        if j<i
            M(i,j)=M(j,i);
        end
    end
end


end

