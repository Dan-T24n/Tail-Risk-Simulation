clear; clc;

%%  generate bivariate t-dist with Gaussian copula

%global params
dim = 2;
m = 1000;
df_m = 5;   %degree of freedom of marginal t-dist
rho = 0.5;  %pairwise correlation


Corr_mat = (1-rho)*eye(dim) + rho*ones(dim);    %correlation matrix
mu = zeros(1,dim);     

Z = mvnrnd(mu, Corr_mat, m);    %maginals of Z are normal 
U = normcdf(Z,0,1);     %use the grade (feed marginals into their own cdf) to make uniform marginals
subplot(2,2,1);
plot(U(:,1),U(:,2),'.');
title(['Gaussian copula,rho = ',num2str(rho)]);
xlabel('U1');
ylabel('U2');
set(gca,'FontSize',15)


X = [tinv(U(:,1),df_m) tinv(U(:,2),df_m)];      %inverse transform of marginal cdf^(-1) to desired cdf
subplot(2,2,2);
plot(X(:,1),X(:,2),'.');
title(['Gaussian copula,rho = ',num2str(rho)]);
xlabel(['X1 ~ t(',num2str(df_m),')']);
ylabel(['X2 ~ t(',num2str(df_m),')']);
q = quantile(X(:,1),0.99)
xlim([-q,q]);
ylim([-q,q]);
set(gca,'FontSize',15)

%%  generate bivariate t-dist with t-copula

nu = 20;     %degree of freedom of the copula

T = mvtrnd(Corr_mat,nu, m);    %maginals of T are t-dist
U = tcdf(T,nu);     %use same specs cdf to get the grade
subplot(2,2,3);
plot(U(:,1),U(:,2),'.');
title(['t-copula,rho = ',num2str(rho)]);
xlabel('U1');
ylabel('U2');
set(gca,'FontSize',15)


X = [tinv(U(:,1),df_m) tinv(U(:,2),df_m)];     %inverse transform of marginal cdf^(-1) to desired cdf
subplot(2,2,4);
plot(X(:,1),X(:,2),'.');
title(['t-copula,rho = ',num2str(rho)]);
xlabel(['X1 ~ t(',num2str(df_m),')']);
ylabel(['X2 ~ t(',num2str(df_m),')']);
q = quantile(X(:,1),0.99)
xlim([-q,q]);
ylim([-q,q]);
set(gca,'FontSize',15)
