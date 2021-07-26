% Pragya Patel
% 17807477

function res = L22norm(x,xnew)
% This function calculates the L2-norm
%
% Inputs
%   matrix initial (phi_old), matrix final (phi_new)
%   size 
% Output
%   L2-norm (res)

sz = size(x);
Nx=sz(1)-2; Nz=sz(2)-2;

diff = xnew - x;
sum = 0;
for k = 2:Nz+1
        for i = 2:Nx+1
            sum = sum + (diff(i,k))^2;
        end
end
res = sqrt(sum/(Nx*Nz));
end