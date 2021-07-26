% Pragya Patel (17807477)
% Flow past a Square Cylinder
% RHS FOR PRESSURE CORRECTION

function utemp = urhs(u,v,w) %,p,T,xc,ye,ze)
% This function calculates the rhs of the poisson relation
% at for every interior i,j,k at every rk-substep
%
% Inputs
%   u,v,w at the previous RK3 step
% Output
%   u_temp for the next RK3 step

global Re
% Initialize
sz = size(u);
nxp2 = sz(1); nyp2 = sz(2); nzp2 = sz(3);
utemp = zeros(nxp2,nyp2,nzp2);
dx = 1/(nxp2-2);
dy = 1/(nyp2-2);
dz = 1/(nzp2-2);

% Non-linear terms
for i = 2:nxp2-1
    for j = 2:nyp2-1
        for k = 2:nzp2-1
            % d(uu)/dx
            uav1 = 0.5*(u(i,j,k)+u(i+1,j,k));
            uav2 = 0.5*(u(i,j,k)+u(i-1,j,k));
            utemp(i,j,k) = utemp(i,j,k) + (uav1^2 - uav2^2)/dx;
            % d(uv)/dy
            uav1 = 0.5*(u(i,j,k)  + u(i,j+1,k));
            uav2 = 0.5*(u(i,j,k)  + u(i,j-1,k));
            vav1 = 0.5*(v(i,j,k)  + v(i+1,j,k));
            vav2 = 0.5*(v(i,j-1,k)+ v(i+1,j-1,k));
            utemp(i,j,k) = utemp(i,j,k) + (uav1*vav1 - uav2*vav2)/dy;
            % d(uw)/dz
            uav1 = 0.5*(u(i,j,k)  + u(i,j,k+1));
            uav2 = 0.5*(u(i,j,k)  + u(i,j,k-1));
            wav1 = 0.5*(w(i,j,k)  + w(i+1,j,k));
            wav2 = 0.5*(w(i,j,k-1)+ w(i+1,j,k-1));
            utemp(i,j,k) = utemp(i,j,k) + (uav1*wav1 - uav2*wav2)/dz;    
        end
    end
end

% Sub from Linear terms
for i = 2:nxp2-1
    for j = 2:nyp2-1
        for k = 2:nzp2-1    
            d2u = (u(i+1,j,k)-2*u(i,j,k)+u(i-1,j,k))/(dx)^2 ...
                + (u(i,j+1,k)-2*u(i,j,k)+u(i,j-1,k))/(dy)^2;
            utemp(i,j,k) = -utemp(i,j,k) + d2u/Re;
        end
    end
end
end

% For non-uniform grid
% Add these lines in each loop
%             dx = xc(i+1)-xc(i);
%             dy = ye(j)-ye(j-1);
%             dz = ze(k)-ze(k-1);
% Further,
% Linear Terms:
% d2u/dx2 = (delu/delx (at i+1) - (at i))/(xe(i+1)-xe(i))
% where
%       delu/delx (at i+1) = (u(i+1) - u(i))/(xe(i+1)-xe(i))
%       delu/delx (at i) = (u(i) - u(i-1))/(xe(i)-xe(i-1))