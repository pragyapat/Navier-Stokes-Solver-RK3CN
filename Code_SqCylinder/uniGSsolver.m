% Pragya Patel
% 17807477
% Convergent Gauss Seidel Solution (uses updatebcp.m)

function x = uniGSsolver(dx,rhs,tol,steps)
% Inputs
%   coefficient matrix (C), PDE RHS (rhs)
% Output
%   Converged value (x)
s = size(rhs);
nxp2=s(1); nyp2=s(2); nzp2=s(3);

x = zeros(nxp2,nyp2,nzp2); x0 = x;
for t = 1:steps
    % Boundary
    x = updatebc(x,'p');
    for k = 2:nzp2-1
        for j = 2:nyp2-1
            for i = 2:nxp2-1
                % Interior
                x(i,j,k) = (1/6)*(-dx^2*rhs(i,j,k)+x(i+1,j,k)+x(i-1,j,k) ...
                           +x(i,j+1,k)+x(i,j-1,k)+x(i,j,k+1)+x(i,j,k-1));
            end
        end
    end
    
    % Check convergence
    rnew = L2norm(x0,x);
    disp([num2str(t) ' ' num2str(rnew)]);
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

% Pragya Patel
% 17807477
% Convergent Gauss Seidel Solution (uses updatebcp.m)

% function x = GSsolver(dx,rhs,tol,steps)
% % Inputs
% %   coefficient matrix (C), PDE RHS (rhs)
% % Output
% %   Converged value (x)
% s = size(rhs);
% nxp2=s(1); nyp2=s(2); nzp2=s(3);
% 
% x0 = zeros(nxp2,nyp2,nzp2);
% for t = 1:steps
%     % Boundary
%     x = updatebc(x0,'p');
%     for k = 2:nzp2-1
%         for j = 2:nyp2-1
%             for i = 2:nxp2-1
%                 % Interior
%                 x(i,j,k) = (1/6)*(-dx^2*rhs(i,j,k)+x0(i+1,j,k)+x(i-1,j,k) ...
%                            +x0(i,j+1,k)+x(i,j-1,k)+x0(i,j,k+1)+x(i,j,k-1));
%             end
%         end
%     end
%     
%     % Check convergence
%     rnew = L2norm(x0,x);
%     disp([num2str(t) ' ' num2str(rnew)]);
%     if rnew < tol
%         steps = t;
%         disp(['Converging at ', num2str(t) ' with res = ', num2str(rnew)])
%         break
%     elseif rnew > 1000
%         disp(['GS diverging at ', num2str(t) ' with res = ', num2str(rnew)])
%         break
%     elseif t > 1 && rnew == r
%         disp('F')
%         break
%     else
%         % Update
%         x0 = x;
%         r = rnew;
%     end
%     % Next step
% end
% if t == steps
%     disp(['Insufficient GS tmax ', num2str(t) ' with res = ', num2str(r)])
% end
% x = updatebc(x0,'p');
% end