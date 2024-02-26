
clear; clc;

%% global params

n = 1000;    % #firms
t = 2000;   % #time horizon (20*months)
%df_m = 3;   % degree of freedom of marginal t-dist
%nu = 3;     % degree of freedom of t-copula
%rho = 0.2;  % pairwise correlation


%% run simulations for df_m and nu

%set path
root = pwd;
folder = '/Users/macbookpro/Desktop/Google Drive/Smooth Tail Risk RA/coding/data';

iter = 10;

k = 3; %base DoF 
step = 3; %increment for DoF
size = 10; %size of test range

D = zeros(10,10);

Rho = [0.5];


for rho = Rho

    for nu = k:step:size*step
    
        for df_m = k:step:size*step
        
            sum = 0;
        
            for i = 1:iter
            
                delta = Simulate(n,t,rho,nu,df_m);
            
                sum = sum + delta;    
             
            end
        
        D(nu/step,df_m/step) = sum/iter;
        
        end
   
    end

% export matrix of deltas


name = sprintf('Delta, rho=%.1f, DoF=%.0f',rho,size*step);
ext = '.csv';

filename = strcat(name,ext);
store = sprintf('DoF=%.0f',size*step);
path = fullfile(folder,store,filename);

csvwrite(path,D)

cd(root)

end


%% simulations cumul returns, for df_m and nu

%set path
root = pwd;
folder = '/Users/macbookpro/Desktop/Google Drive/Smooth Tail Risk RA/coding/data';

iter = 10;

k = 3; %base DoF 
step = 3; %increment for DoF
size = 10; %size of test range

Rho = [0.1 0.2 0.5 0.7];

%preload holders
D1 = zeros(size,size);
D2 = zeros(size,size);


for rho = Rho

    for nu = k:step:size*step
    
        for df_m = k:step:size*step
        
            sum_Kelly = 0;
            sum_Smooth = 0;
        
            for i = 1:iter
            
               [Kelly_Month,Smooth_Month] = Cumul_Simulate(n,t,rho,nu,df_m);
            
                sum_Kelly = sum_Kelly + Kelly_Month;
                sum_Smooth = sum_Smooth + Smooth_Month;
             
            end
        
        D1(nu/step,df_m/step) = sum_Kelly/iter;
        D2(nu/step,df_m/step) = sum_Smooth/iter;
        
        end
   
    end

% export matrix of deltas

name1 = sprintf('Kelly, rho=%.1f, DoF=%.0f',rho,size*step);
name2 = sprintf('Smooth, rho=%.1f, DoF=%.0f',rho,size*step);
ext = '.csv';

filename1 = strcat(name1,ext);
filename2 = strcat(name2,ext);

path = fullfile(folder,filename1);
csvwrite(path,D1)

path = fullfile(folder,filename2);
csvwrite(path,D2)

cd(root); 

end

%about 30mins/rho per n=500

%% Learn array data structure
