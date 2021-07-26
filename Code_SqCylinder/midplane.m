% Pragya Patel
% 17807477
% Multigrid Solver: Midplane plotter

function midplane(data,s)
% This function calculates the planar phi
% at j = j(midplane) and returns a contour plot
%
% Inputs
%   matrix to be plotted (data), string for plot title (s)
% Output
%   plot data(i,k) at j = Ny/2

L = 1;
sz = size(data);
nxp2 = sz(1); nyp2 = sz(2); nzp2 = sz(3);
dx = L/(nxp2-2); dz = L/(nzp2-2);

planar = zeros(nxp2-2,nzp2-2);
x = zeros(nxp2-2,1); z = zeros(nzp2-2,1);
for i = 2:nxp2-1
    for k = 2:nzp2-1
        planar(k-1,i-1) = data(i,(nyp2-2)/2,k);
        z(k-1) = dz*(k-1);
    end
    x(i-1) = dx*(i-1); 
end
lev = -1:0.004:1;
figure
hold on
contourf(x,z,planar,lev)
xlabel('x')
ylabel('z')
colorbar
title([s ' at j=j_{mid}'])
hold off
drawnow
end