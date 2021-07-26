% Pragya Patel
% 17807477
% Time Step Calculator

function dt = deltat(u,v,w,h)
% This function calculates the time step
% according to CFL criteria
%
% Inputs
%   velocity matrices (u,v,w) and dimensions (dx,dy,dz)
% Output
%   dt according to the CFL criteria

global CFL
dx = h(1); dy = h(2); dz = h(3);
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