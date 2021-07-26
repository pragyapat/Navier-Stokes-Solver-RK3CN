% Pragya Patel (17807477)
% Flow past a Square Cylinder
% RK3CN: RK3 part

function [urk1,vrk1,wrk1,prk1,qu1,qv1,qw1] ...
    = rk3(sub,d,sqn,urk0,vrk0,wrk0,qu0,qv0,qw0)
% This function calculates the variables after one RK substep
% Inputs
%   substep number (sub), grid dimensions (d = [dx,dy,dz]) and
%   matrices of the variables obtained in the previous substep (var0)
% Output
%   variables for the next rk-step (var1)

% Initialize
sz = size(urk0);
nxp2 = sz(1); dx = d(1);
nyp2 = sz(2); dy = d(2);
nzp2 = sz(3); dz = d(3);
O = zeros(nxp2,nyp2,nzp2);
urk1 = O; vrk1 = O; wrk1 = O;
ustar = O; vstar = O; wstar = O;
ustar = updatebc(sqn,ustar,'u');
vstar = updatebc(sqn,vstar,'v');
wstar = updatebc(sqn,wstar,'w');

% Constants
global dt Re
[Crk,dtrk] = coeff_rk3(sub,dt);
a = -dtrk/(Re*dz^2);
b = 1+2*dtrk/(Re*dz^2);
c = -dtrk/(Re*dz^2);

% Explicit terms
qu1 = dt*urhs(d,urk0,vrk0,wrk0) + Crk(1)*qu0;
qv1 = dt*vrhs(d,urk0,vrk0,wrk0) + Crk(1)*qv0;
qw1 = dt*wrhs(d,urk0,vrk0,wrk0) + Crk(1)*qw0;
utemp = urk0 + Crk(2)*qu1;
vtemp = vrk0 + Crk(2)*qv1;
wtemp = wrk0 + Crk(2)*qw1;

% Intermediate velocity fields via Thomas Algorithm
% USTAR
for i = 2:nxp2-2
    for j = 2:nyp2-1
        ustar(i,j,:) = velstar(a,b,c,utemp(i,j,:),'u');
    end
end
% VSTAR
for i = 2:nxp2-1
    for j = 2:nyp2-2
        vstar(i,j,:) = velstar(a,b,c,vtemp(i,j,:),'v');
    end
end
% WSTAR
for i = 2:nxp2-1
    for j = 2:nyp2-1
        wstar(i,j,1:nzp2-1) = velstar(a,b,c,wtemp(i,j,:),'w');
    end
end
ustar = updatebc(sqn,ustar,'u');
vstar = updatebc(sqn,vstar,'v');
wstar = updatebc(sqn,wstar,'w');
% [ustar,vstar,wstar] = fixit(sqn,ustar,vstar,wstar);

% Set up and solve for pressure correction
utemp = O; vtemp = O; wtemp = O; % reusing temp variables
for i = 2:nxp2-1
    for j = 2:nyp2-1
        for k = 2:nzp2-1
            utemp(i,j,k) = (ustar(i,j,k)-ustar(i-1,j,k))/dx;
            vtemp(i,j,k) = (vstar(i,j,k)-vstar(i,j-1,k))/dy;
            wtemp(i,j,k) = (wstar(i,j,k)-wstar(i,j,k-1))/dz;
        end
    end
end
rhsp = (utemp+vtemp+wtemp)/dtrk;    % RHS for pressure-poisson
prk1 = GSsolver(d,sqn,rhsp,1e-5,1e4);  % Solve Ax=b

% Update the intermediate velocities using correction terms
for i = 2:nxp2-1
    for j = 2:nyp2-1
        for k = 2:nzp2-1
            urk1(i,j,k) = ustar(i,j,k) ...
                - dtrk*(prk1(i+1,j,k)-prk1(i,j,k))/dx; % dxc
            vrk1(i,j,k) = vstar(i,j,k) ...
                - dtrk*(prk1(i,j+1,k)-prk1(i,j,k))/dy; % dyc
            wrk1(i,j,k) = wstar(i,j,k) ...
                - dtrk*(prk1(i,j,k+1)-prk1(i,j,k))/dz; % dzc
        end
    end
end
urk1 = updatebc(sqn,urk1,'u');
vrk1 = updatebc(sqn,vrk1,'v');
wrk1 = updatebc(sqn,wrk1,'w');
% [urk1,vrk1,wrk1] = fixit(sqn,urk1,vrk1,wrk1);
end