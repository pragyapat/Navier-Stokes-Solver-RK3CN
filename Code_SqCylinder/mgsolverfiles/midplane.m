% Pragya Patel
% 17807477
% Multigrid Solver: Midplane plotter

function midplane(phi,Nx,Ny,Nz,s)
% This function calculates the planar phi
% at j = j(midplane) and returns a contour plot
%
% Inputs
%   matrix to be plotted (phi), interior size (Nx,Ny,Nz)
% Output
%   plot phi(i,k) at j = Ny/2

phiplanar = zeros(Nx,Nz);
for i = 2:Nx+1
    for k = 2:Nz+1
        phiplanar(k-1,i-1) = phi(i,Ny/2,k);
    end
end

figure 
hold on
contourf(phiplanar)
colorbar
title(['\phi(x,z) for ' s ' at j=j_{mid}'])
hold off

end