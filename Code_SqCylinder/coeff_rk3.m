% Pragya Patel
% 17807477
% Coefficients for RK3 method

function [Cm,dtrk] = coeff_rk3(m,dt)
% This function generates a coefficient matrix
% for a uniform 3D grid
%
% Inputs
%   RK3 step (m)
%   m can take the values 1, 2 and 3
% Output
%   Coefficients for this step (Cm)
% Notation
%   C sup m sub n
%   superscript is the rk step
%   subscript is the nth coefficient
%   each step (m) has 3 coefficients (n=1,2,3)
%   Cm(i) = ith coefficient for mth step

if m == 1
    Cm = [0, 1/3, 0.3333];
elseif m == 2
    Cm = [-5/9, 15/16, 0.4167];
elseif m == 3
    Cm = [-153/128, 8/15, 0.2500];
else
    disp('Error: m can only take values = 1,2 or 3')
end

% time step for mth rk step
dtrk = Cm(3)*dt;
end