% Pragya Patel
% 17807477
% Grid Maker: For uniform 3D grids

function C = coeff_uni(Nx,Ny,Nz,Lx,Ly,Lz)
% This function generates a coefficient matrix
% for a uniform 3D grid
%
% Inputs
%   Size (Nx Ny Nz Lx Ly Lz)
% Output
%   Coefficient matrices for uniform mesh (C)

dx = Lx/(Nx-1); dy = Ly/(Ny-1); dz = Lz/(Nz-1);
s = [Nx,Ny,Nz,dx,dy,dz,Lx,Ly,Lz];

% Coefficient Matrices
c = ones(Nx+2,Ny+2,Nz+2)*(-2)*(1/(dx*dx) + 1/(dy*dy) + 1/(dz*dz));
ce = ones(Nx+2,Ny+2,Nz+2)/(dx*dx);
cn = ones(Nx+2,Ny+2,Nz+2)/(dy*dy);
ca = ones(Nx+2,Ny+2,Nz+2)/(dz*dz);
cw = ce; cs = cn; cb = ca;

f1 = 'c'; v1 = c;   % coefficient center
f2 = 'ce'; v2 = ce; % coefficient east
f3 = 'cn'; v3 = cn; % coefficient north
f4 = 'ca'; v4 = ca; % coefficient above
f5 = 'cw'; v5 = cw; % coefficient west
f6 = 'cs'; v6 = cs; % coefficient south
f7 = 'cb'; v7 = cb; % coefficient below
f8 = 's'; v8 = s; % size and dimensions

% OUTPUT
C = struct(f1,v1,f2,v2,f3,v3,f4,v4,f5,v5,... 
    f6,v6,f7,v7,f8,v8);
end