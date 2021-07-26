% Pragya Patel
% 17807477
% Convergent Gauss Seidel Solution (uses updatebc.m)
% For uniform grids

function x = GSsolver(d,rhs,tol,steps)
% Inputs
%   grid dimensions (d = [dx,dy,dz]), PDE RHS (rhs),
%   desired tolerance (tol) and maximum steps (steps)
% Output
%   Converged value (x)
% Parameter l = SOR Relaxation

l = 0.7;
s = size(rhs);
nxp2=s(1); dx = d(1); cx = 1/(dx)^2;
nyp2=s(2); dy = d(2); cy = 1/(dy)^2;
nzp2=s(3); dz = d(3); cz = 1/(dz)^2;
c = -2*(cx+cy+cz);
x = zeros(nxp2,nyp2,nzp2);
x0 = x;

for t = 1:steps
    % Boundary
    x = updatebc(x,'p');
    % Interior
    for k = 2:nzp2-1
        for j = 2:nyp2-1
            for i = 2:nxp2-1
                expn = rhs(i,j,k) - (...
                    cx*(x0(i+1,j,k)+x(i-1,j,k)) + ...
                    cy*(x0(i,j+1,k)+x(i,j-1,k)) + ...
                    cz*(x0(i,j,k+1)+x(i,j,k-1)));
                x(i,j,k) = l*(expn/c) + (1-l)*x0(i,j,k);
            end
        end
    end
    
    % Check convergence
    rnew = L2norm(x0,x);
%     disp([num2str(t) ' ' num2str(rnew)]);
    if rnew < tol
        disp(['Converging at ', num2str(t) ' with res = ', num2str(rnew)])
        break
    elseif rnew > 1000
        disp(['GS diverging at ', num2str(t) ' with res = ', num2str(rnew)])
        break
    elseif t > 1 && rnew == r
        disp('GS residual constant')
        break
    else
        % Update
        x0 = x;
        r = rnew;
    end
    % Next step
end
if t == steps
    disp(['Insufficient GS tmax ', num2str(t) ' with res = ', num2str(r)])
end
x = updatebc(x,'p');
end