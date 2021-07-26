% Pragya Patel (17807477)
% Flow past a Square Cylinder
% RK3CN: EXPLICIT TERMS

function vtemp = vrhs(d,u,v,w)
% This function calculates the explicit rhs
% to be used in Thomas Alg. for vstar
%
% Inputs
%   grid dimensions (d = [dx,dy,dz]), 
%   and u,v,w at the previous RK3 step
% Output
%   explicit terms for rhs (into qv1->vstar->rhsp->v)
global Re
% Initialize
sz = size(v);
nxp2 = sz(1); dx = d(1);
nyp2 = sz(2); dy = d(2);
nzp2 = sz(3); dz = d(3);
vtemp = zeros(nxp2,nyp2,nzp2);

% Non-linear terms
for i = 2:nxp2-1
    for j = 2:nyp2-1
        for k = 2:nzp2-1
            % d(uv)/dx
            uav1 = 0.5*(u(i,j,k)  + u(i,j+1,k));
            uav2 = 0.5*(u(i-1,j,k)+ u(i-1,j+1,k));
            vav1 = 0.5*(v(i,j,k)  + v(i+1,j,k));
            vav2 = 0.5*(v(i,j,k)  + v(i-1,j,k));
            vtemp(i,j,k) = vtemp(i,j,k) + (uav1*vav1 - uav2*vav2)/dx;
            % d(vv)/dy
            vav1 = 0.5*(v(i,j,k)+v(i,j+1,k));
            vav2 = 0.5*(v(i,j,k)+v(i,j-1,k));
            vtemp(i,j,k) = vtemp(i,j,k) + (vav1^2 - vav2^2)/dy;
            % d(vw)/dz
            vav1 = 0.5*(v(i,j,k)  + v(i,j,k+1));
            vav2 = 0.5*(v(i,j,k)  + v(i,j,k-1));
            wav1 = 0.5*(w(i,j,k)  + w(i,j+1,k));
            wav2 = 0.5*(w(i,j,k-1)+ w(i,j+1,k-1));
            vtemp(i,j,k) = vtemp(i,j,k) + (vav1*wav1 - vav2*wav2)/dz;    
        end
    end
end

% Sub from Linear terms
for i = 2:nxp2-1
    for j = 2:nyp2-1
        for k = 2:nzp2-1    
            d2v = (v(i+1,j,k)-2*v(i,j,k)+v(i-1,j,k))/(dx)^2 ...
                + (v(i,j+1,k)-2*v(i,j,k)+v(i,j-1,k))/(dy)^2;
            vtemp(i,j,k) = -vtemp(i,j,k) + d2v/Re;
        end
    end
end
end