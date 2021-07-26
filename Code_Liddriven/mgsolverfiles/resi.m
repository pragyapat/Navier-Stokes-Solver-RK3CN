% Pragya Patel
% 17807477
% Multigrid Solver: Residual Calculation

function r = resi(x,C,rhs)
% This function calculates the residual of the
% form r = b - Ax
%
% Inputs
%   relevant coefficient matrix A (C), the variable (x)
%   rhs of the exact equation Ax=b (rhs)
% Output
%   residual (r)

pc = C.c;
pce = C.ce;
pcn = C.cn;
pca = C.ca;
pcw = C.cw;
pcs = C.cs;
pcb = C.cb;
s = C.s; % internal size
Nx = s(1); Ny = s(2); Nz = s(3);

% Initialize
r = zeros(Nx+2,Ny+2,Nz+2);

% Residual
for i = 2:Nx+1
    for j = 2:Ny+1
        for k = 2:Nz+1
            r(i,j,k) = rhs(i,j,k) - pc(i,j,k)*x(i,j,k) ...
                - pce(i,j,k)*x(i+1,j,k) - pcw(i,j,k)*x(i-1,j,k) ...
                - pcn(i,j,k)*x(i,j+1,k) - pcs(i,j,k)*x(i,j-1,k) ...
                - pca(i,j,k)*x(i,j,k+1) - pcb(i,j,k)*x(i,j,k-1);
        end
    end
end

% Note: Residual on boundaries is zero
% Update boundary nodes using updatebcr
r = updatebcr(r);
end