% Pragya Patel
% 17807477
% Gauss Seidel Smoother for errors (uses updatebcr.m)

function x = GSr(x_old,C,rhs,steps)
% This function gives an approximate, partially converged
% solution for a matrix equation of the residual/error
%
% Inputs
%   coefficient matrix (C), initial value (x_old)
%   desired number of gauss seidel iterations (steps)
%   and the rhs (rhs)
% Output
%   partially converged value (x_new)

pc = C.c;
pce = C.ce;
pcn = C.cn;
pca = C.ca;
pcw = C.cw;
pcs = C.cs;
pcb = C.cb;
s = C.s;
Nx=s(1); Ny=s(2); Nz=s(3); % dy=s(5); dz=s(6);

for t = 1:steps
    % Boundary
    x = updatebcr(x_old);
    for k = 2:Nz+1
        for j = 2:Ny+1
            for i = 2:Nx+1
                % Interior
                expn = rhs(i,j,k) ...
                    - pce(i,j,k)*x_old(i+1,j,k) - pcw(i,j,k)*x(i-1,j,k) ...
                    - pcn(i,j,k)*x_old(i,j+1,k) - pcs(i,j,k)*x(i,j-1,k) ...
                    - pca(i,j,k)*x_old(i,j,k+1) - pcb(i,j,k)*x(i,j,k-1);
                x(i,j,k) = expn/pc(i,j,k);
            end
        end
    end
    % Update
    x_old = x;
    % Next step
end
x = updatebcr(x_old);
end