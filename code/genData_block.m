function [X] = genData_block(n,b,t,rvec,nu,df_m)

% Generate n*t dimension data, with t-coplua and t-dist marginals
%   ==params==
% n         % #firms
% b         % #blocks
% t         % #time horizon (20*months)
% rho       % pairwise correlation 
% df_m      % degree of freedom of marginal t-dist
% nu        % degree of freedom of t-copula


Rho_mat = block_mat(b,n/b,rvec);    %correlation matrix
T = mvtrnd(Rho_mat,nu,t);     %t-copula, marginal are t-dist(nu)
U = tcdf(T,nu);     %use cdf to get the grade
X = tinv(U,df_m);   %inverse tranform the grades to t-dist(df_m)

end

