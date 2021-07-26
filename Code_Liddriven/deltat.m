% Pragya Patel
% 17807477
% Multigrid Solver: Update BC for the problem

function dt = deltat(u,v,w,dim)
% This function updates the boundary conditions
% as specified in the problem
%
% Inputs
%   velocity matrices (u,v,w) and dimensions (dx,dy,dz)
% Output
%   dt according to the CFL criteria

global CFL
dx = dim(1); dy = dim(2); dz = dim(3);
dtx = 100; dty = 100; dtz = 100;
umax = max(max(max(u)));
vmax = max(max(max(v)));
wmax = max(max(max(w)));

if umax > 0
    dtx = CFL*dx/umax;
elseif vmax > 0
    dty = CFL*dy/vmax;
elseif wmax > 0
    dtz = CFL*dz/wmax;
end

% minimum value
dt = min([dtx,dty,dtz]);
end