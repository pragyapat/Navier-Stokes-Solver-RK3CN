% ME634 End-Semester
% Pragya Patel (17807477)
% Flow past a Square Cylinder
% using RK3-CN Algorithm
clear global
%% Setup
% Parameters
global tol CFL dt U Re
tol = 1e-4;
CFL = 1.2;
U   = 1;
Re  = 20;
steps = 10000;

% Domain
Nx = 178;   nxp2 = Nx+2;
Ny = 32;    nyp2 = Ny+2;
Nz = 80;    nzp2 = Nz+2;

L  = 24; Lx = L; dx = Lx/Nx;
A  = 6;  Ly = A; dy = Ly/Ny;
H  = 10; Lz = H; dz = Lz/Nz;
d  = [dx,dy,dz];

% Square Configuration
B  = 1; nu = U*B/Re;    % Given
nxsq = round(B/dx);     % number of cells within the square in x
nzsq = round(B/dz);     % number of cells within the square in z
Bx = nxsq*dx;           % adjusted dimensions of the square
Bz = nzsq*dz;           % adjusted dimensions of the square

La = 6;                 % Given
LLE = La-(Bx/2);        % Leading edge adjusted location
LBE = (Lz-Bz)/2;        % Bottom edge adjusted location
cnx = LLE/dx;           % Bottom-left corner nx
cnz = LBE/dz;           % Bottom-left corner nz

% Square Nodes
x1 = cnx;      z1 = cnz;
x2 = cnx+nxsq; z2 = cnz+nzsq;
sqn = [x1;  % left
    x2;     % right
    z1;     % bottom
    z2];    % top
disp('Square Nodes')
disp(sqn)

%% Initialize
O  = zeros(nxp2,nyp2,nzp2);
u  = O; v = O; w = O; p = O; t = 0;
qu0 = O; qv0 = O; qw0 = O;
u = updatebc(sqn,u,'u');
v = updatebc(sqn,v,'v');
w = updatebc(sqn,w,'w');

%% RK3 CN Loop
for step = 1:steps
    dt = deltat(u,v,w,d); t = t + dt;
    disp(['Main Step: ' num2str(step) ' Time (s) = ' num2str(t)]);
    
    % RK Substep 1 to 3
    [urk1,vrk1,wrk1,prk1,qu1,qv1,qw1] = rk3(1,d,sqn,u,v,w,qu0,qv0,qw0);
    [urk2,vrk2,wrk2,prk2,qu2,qv2,qw2] = rk3(2,d,sqn,urk1,vrk1,wrk1,qu1,qv1,qw1);
    [urk3,vrk3,wrk3,prk3,qu3,qv3,qw3] = rk3(3,d,sqn,urk2,vrk2,wrk2,qu2,qv2,qw2);
    
    % Check convergence
    resu = L2norm(u,urk3);
    res = resu;
    disp(['Main Residual = ' num2str(res)])
    if res < tol
        disp(['Converged. Final Residual = ' num2str(res)])
        break
    elseif res > 2000
        disp(['Not converging. Main Residual = ' num2str(res)])
        break
    else
        % Update
        u = urk3;
        v = vrk3;
        w = wrk3;
        p = prk3;
    end
end

% DUMP
% % Square Nodes
% nblx = cnx;      nblz = cnz;
% nbrx = cnx+nxsq; nbrz = cnz;
% ntlx = cnx;      ntlz = cnz+nzsq;
% ntrx = cnx+nxsq; ntrz = cnz+nzsq;
% sqn = [nblx,nblz;   % Bottom left
%     nbrx,nbrz;      % Bottom right
%     ntlx,ntlz;      % Top left
%     ntrx,ntrz];    % Top right
% disp('Square Nodes')
% disp(sqn);
% % disp(['(' num2str(ntlx) ',' num2str(ntlz) ')--' '(' num2str(ntrx) ',' num2str(ntrz) ')'])
% % disp(['(' num2str(nblx) ',' num2str(nblz) ')--' '(' num2str(nbrx) ',' num2str(nbrz) ')'])