function [Kelly, Smooth] = ComputeTR(data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

t = size(data,1);
M = t/20;

Kelly = zeros(1,M);
Smooth = zeros(1,M);

for i = 1 : M
    
idx = i*20-linspace(19,0,20);
X = data (idx,:);

Kelly(i) = CSTR(X);
Smooth(i) = SmoothCSTR(X);

end

