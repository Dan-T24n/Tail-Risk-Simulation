
clear; clc;

%% import file


%set DoF
k=103;

name = sprintf('compare rho, DoF=%.0f',k);
%path = filepath

ext = '.csv';
filename = strcat(name,ext);

Delta = importfile_rho(filename);
clear filename;

j = linspace(1,11,11);
Dof = 3 + (j-1)*10;
rho = linspace(0,0.5,11);

%% 3D plot

figure2=figure();
mesh(Dof,rho,Delta);
title('Differences in pMSE','Interpreter','latex');
xlabel('Degree of freedom','Interpreter','latex');
ylabel('Correlation $(\rho)$','Interpreter','latex');
view(40,14);
set(gca,'FontSize',16)


%% Color map

%drop 1st column (too high values) to fit colors
Delta2 = Delta(:,2:end);
j = linspace(2,11,10);
Dof = 3 + (j-1)*10;

%change color map
c=jet;
k=0.2;

for i=1:6 %or i = vector indexes  
     c(i,3)= c(i,3) - k/(i^1.2);
end


figure1=figure();
colormap(c);
% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
% Create colormap
pcolor(Dof,rho,Delta2)
title('Differences in pMSE ','Interpreter','latex')
shading interp;
% Create xlabel
xlabel('DoF','Interpreter','latex','FontSize',16);
xlim([13,103]);
% Create ylabel
ylabel('Correlation $(\rho)$','Interpreter','latex','FontSize',16);
grid(axes1,'on');
% Create colorbar
colorbar('peer',axes1);
hold off
set(gca,'FontSize',16)
