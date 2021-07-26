% Pragya Patel (17807477)
% Flow past a Square Cylinder
% RK3CN: EXPLICIT TERMS

function wtemp = wrhs(d,u,v,w)
% This function calculates the explicit rhs
% to be used in Thomas Alg. for wstar
%
% Inputs
%   grid dimensions (d = [dx,dy,dz]), 
%   and u,v,w at the previous RK3 step
% Output
%   explicit terms for rhs (into wu1->wstar->rhsp->w)
global Re
% Initialize
sz = size(w);
nxp2 = sz(1); dx = d(1);
nyp2 = sz(2); dy = d(2);
nzp2 = sz(3); dz = d(3);
wtemp = zeros(nxp2,nyp2,nzp2);

% Non-linear terms
for i = 2:nxp2-1
    for j = 2:nyp2-1
        for k = 2:nzp2-1
            % d(uw)/dx
            uav1 = 0.5*(u(i,j,k)  + u(i,j,k+1));
            uav2 = 0.5*(u(i-1,j,k)+ u(i-1,j,k+1));
            wav1 = 0.5*(w(i,j,k)  + w(i+1,j,k));
            wav2 = 0.5*(w(i,j,k)  + w(i-1,j,k));
            wtemp(i,j,k) = wtemp(i,j,k) + (uav1*wav1 - uav2*wav2)/dx;
            % d(vw)/dy
            vav1 = 0.5*(v(i,j,k)  + v(i,j,k+1));
            vav2 = 0.5*(v(i,j-1,k)+ v(i,j-1,k+1));
            wav1 = 0.5*(w(i,j,k)  + w(i,j+1,k));
            wav2 = 0.5*(w(i,j,k)  + w(i,j-1,k));
            wtemp(i,j,k) = wtemp(i,j,k) + (vav1*wav1 - vav2*wav2)/dy;
            % d(ww)/dz
            wav1 = 0.5*(w(i,j,k)+w(i,j,k+1));
            wav2 = 0.5*(w(i,j,k)+w(i,j,k-1));
            wtemp(i,j,k) = wtemp(i,j,k) + (wav1^2 - wav2^2)/dz;
        end
    end
end

% Sub from Linear terms
for i = 2:nxp2-1
    for j = 2:nyp2-1
        for k = 2:nzp2-1    
            d2w = (w(i+1,j,k)-2*w(i,j,k)+w(i-1,j,k))/(dx)^2 ...
                + (w(i,j+1,k)-2*w(i,j,k)+w(i,j-1,k))/(dy)^2;
            wtemp(i,j,k) = -wtemp(i,j,k) + d2w/Re;
        end
    end
end
end