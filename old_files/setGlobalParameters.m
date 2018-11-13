function [kkk1,kkk2,lamb,kapp1,kapp2,degP,degN]=setGlobalParameters()
% All coeficients used in the model, they are defined as global ones in 
%   order to use through functions in this project

kkk1 = 1; kkk2 = 0.5; % diffusion coefficient
lamb = 1e3; % penalty coefficient
kapp1 = 0.5; kapp2 = 0.5; % kap1+kap2 = 1

% quadrature's information
degP = 2; % 2 Gaussian quadrature points (for polinomial function)
degN = 3; % Gaussian quadrature points (for complicated function)
% deg-npoints : 1-1, 2-3, 3-4, 4-6, 5-7, 6-12
end