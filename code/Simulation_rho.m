
clear; clc;

%% global params

n = 500;    % #firms
t = 2000;   % #time horizon (20*months)
%df_m = 3;   % degree of freedom of marginal t-dist
%nu = 3;     % degree of freedom of t-copula
%rho = 0.2;  % pairwise correlation


%% run simulations for rho, same df_m and nu

iter = 20; %nb of iterations for each params vector

D = zeros(11,11); %store in matrix

Rho = linspace(0,0.5,11);

for i = 1 : 11
    
    rho = Rho(i);

    for j = 1:11
        
        k = 3 + (j-1)*10;
        nu = k;
        df_m = k;
   
        sum = 0;
        
            for h = 1:iter
            
                delta = Simulate(n,t,rho,nu,df_m);
            
                sum = sum + delta;    
             
            end
            
        D(i,j) = sum/iter;
            
    end

end

% export matrix of deltas

root = pwd;
name = sprintf('compare rho, DoF=%.0f',k);
ext = '.csv';

filename = strcat(name,ext);

path = fullfile(root,filename);

csvwrite(path,D)