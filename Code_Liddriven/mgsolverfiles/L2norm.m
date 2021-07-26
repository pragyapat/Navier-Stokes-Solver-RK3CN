% Pragya Patel
% 17807477
% Multigrid Solver: L2 norm

function res = L2norm(x,xnew)
% This function calculates the L2-norm
%
% Inputs
%   matrix initial (phi_old), matrix final (phi_new)
%   size 
% Output
%   L2-norm (res)

sz = size(x);
Nx=sz(1)-2; Ny=sz(2)-2; Nz=sz(3)-2;

diff = xnew - x;
sum = 0;
for k = 2:Nz+1
    for j = 2:Ny+1
        for i = 2:Nx+1
            sum = sum + (diff(i,j,k))^2;
        end
    end
end
res = sqrt(sum/(Nx*Ny*Nz));
end