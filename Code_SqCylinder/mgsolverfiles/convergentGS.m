% Pragya Patel
% 17807477
% Convergent Gauss Seidel Solution (uses updatebcp.m)

function [xnew,steps] = convergentGS(x,C,rhs)
% This function gives an approximate, partially converged
% solution for a matrix equation of the form Ax=B
%
% Inputs
%   coefficient matrix (C), initial value (x_old)
%   desired number of gauss seidel iterations (steps)
%   and the rhs (rhs)
% Output
%   Converged value (x)
tol = 1.0e-4;
tmax = 10000;

pc = C.c;
pce = C.ce;
pcn = C.cn;
pca = C.ca;
pcw = C.cw;
pcs = C.cs;
pcb = C.cb;
s = C.s;
Nx=s(1); Ny=s(2); Nz=s(3);

for t = 1:tmax
    % Boundary
    xnew = updatebcp(x,C.s);
    for k = 2:Nz+1
        for j = 2:Ny+1
            for i = 2:Nx+1
                % Interior
                expn = rhs(i,j,k) ...
                    - pce(i,j,k)*x(i+1,j,k) - pcw(i,j,k)*xnew(i-1,j,k) ...
                    - pcn(i,j,k)*x(i,j+1,k) - pcs(i,j,k)*xnew(i,j-1,k) ...
                    - pca(i,j,k)*x(i,j,k+1) - pcb(i,j,k)*xnew(i,j,k-1);
                xnew(i,j,k) = expn/pc(i,j,k);
            end
        end
    end
    
%     sor = 0.75;
%     xnew = (1-sor)*x + sor*xnew;
    
    % Check convergence
    rnew = L2norm(x,xnew);
%     disp([num2str(t) ' ' num2str(rnew)]);
    if rnew < tol
        steps = t;
        disp(['Converging at ', num2str(t) ' with res = ', num2str(rnew)])
        break
    elseif rnew > 1000
        disp(['GS diverging at ', num2str(t) ' with res = ', num2str(rnew)])
        break
    elseif t > 1 && rnew == r
        disp('F')
        break
    else
        % Update
        x = xnew;
        r = rnew;
    end
    % Next step
end
if t == tmax
    disp(['Insufficient GS tmax ', num2str(t) ' with res = ', num2str(r)])
end
xnew = updatebcp(x,C.s);
end