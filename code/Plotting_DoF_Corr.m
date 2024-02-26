
clear; clc;


%% import file correlation

%set folder for import/export
target = '/Users/macbookpro/Desktop/Google Drive/Smooth Tail Risk RA/coding/figure';
source = '/Users/macbookpro/Desktop/Google Drive/Smooth Tail Risk RA/coding/data';
root = pwd;

%set DoF
Dof=30;

%set rho
rho=0.5;

%retrieve path for data file

name = sprintf('Kelly, rho=%.1f, DoF=%.0f',rho,Dof);
ext = '.csv';
filename = strcat(name,ext);
path = fullfile(source,filename);
Delta1 = importfile(path);
clear filename path store name ext;

step = Dof/10;
k=(step:step:Dof);
nu = k ;
df_m = k;

% import file Smooth_correlation
%retrieve path for data file

name = sprintf('Smooth, rho=%.1f, DoF=%.0f',rho,Dof);
ext = '.csv';
filename = strcat(name,ext);
path = fullfile(source,filename);
Delta2 = importfile(path);
clear filename path store name ext;

step = Dof/10;
k=(step:step:Dof);
nu = k ;
df_m = k;



%% Color map

cd(target);

%change color map
c=jet;
% k=0.01;
% 
% for i=1:8 %or i = vector indexes  
%      c(i,3)= c(i,3) - k/(i^1.01);
% end

% Set the limits of the colorbar
bottom = min(min(min(Delta1)),min(min(Delta2)));
top  = max(max(max(Delta1)),max(max(Delta2)));

% plot Kelly
figure1=figure();
colormap(c);

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create colormap
pcolor(df_m,nu,Delta1)
title(['Kelly correlation, $\rho= $ ' num2str(rho)],'Interpreter','latex')
shading interp;
%add manual colorbar
caxis manual
caxis([bottom top]);
% Create xlabel
xlabel('Marginals DoF $(df_m)$','Interpreter','latex','FontSize',16);
xlim([step Dof])
ylim([step Dof])
% Create ylabel
ylabel('Copula DoF $(\nu)$','Interpreter','latex','FontSize',16);
grid(axes1,'on');
% Create colorbar
colorbar('peer',axes1);
hold off
set(gca,'FontSize',16)

name = sprintf('Kelly-corr, rho=%.1f, DoF=%.0f',rho,Dof);
ext = '.png';
filename = strcat(name,ext);
print(filename,'-dpng','-r300');


% plot Smooth
figure1=figure();
colormap(c);
% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
% Create colormap
pcolor(df_m,nu,Delta2)
title(['Smooth correlation, $\rho= $ ' num2str(rho)],'Interpreter','latex')
shading interp;
%add manual colorbar
caxis manual
caxis([bottom top]);
% Create xlabel
xlabel('Marginals DoF $(df_m)$','Interpreter','latex','FontSize',16);
% Create ylabel
ylabel('Copula DoF $(\nu)$','Interpreter','latex','FontSize',16);
xlim([step Dof])
ylim([step Dof])
grid(axes1,'on');
% Create colorbar
colorbar('peer',axes1);
hold off
set(gca,'FontSize',16)

%export
name = sprintf('Smooth-corr, rho=%.1f, DoF=%.0f',rho,Dof);
ext = '.png';
filename = strcat(name,ext);
print(filename,'-dpng','-r300');

cd(root);

%% Relative color map

%standardize values

DD = Delta2 - Delta1;

cd(target);

% plot Smooth
figure1=figure();
colormap(c);
% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
% Create colormap
pcolor(df_m,nu,DD)
title(['Difference in correlation, $\rho= $ ' num2str(rho)],'Interpreter','latex')
shading interp;
%add manual colorbar
caxis manual
caxis([bottom top]);
% Create xlabel
xlabel('Marginals DoF $(df_m)$','Interpreter','latex','FontSize',16);
% Create ylabel
ylabel('Copula DoF $(\nu)$','Interpreter','latex','FontSize',16);
xlim([step Dof])
ylim([step Dof])
grid(axes1,'on');
% Create colorbar
colorbar('peer',axes1);
hold off
set(gca,'FontSize',16)


%export
name = sprintf('Delta-corr, rho=%.1f, DoF=%.0f',rho,Dof);
ext = '.png';
filename = strcat(name,ext);
print(filename,'-dpng','-r300');

cd(root);