addpath(genpath('mgsolverfiles'))

global Re U
U = 1; L = 1; Lx = L; Ly = L; Lz = L;
N = 32; Re = 100; Nx = N; Ny = N; Nz = N; % Case 1
dx = Lx/Nx; dy = Ly/Ny; dz = Lz/Nz;
% N = 64; Re = 400; Nx = N; Ny = N; Nz = N; % Case 2
% N = 128; Re = 1000; Nx = N; Ny = N; Nz = N; % Case 3
nu = U*L/Re;
nxp2 = Nx+2; nyp2 = Ny+2; nzp2 = Nz+2;
O = zeros(nxp2,nyp2,nzp2);

%% Initialize
u = O; v = O; w = O; p = O; t = 0;
qu0 = O; qv0 = O; qw0 = O;
u = updatebc(u,'u');
v = updatebc(v,'v');
w = updatebc(w,'w');

global tol CFL dt
tol = 1e-6;
CFL = 1.2;
steps = 10000;

%% RK3 CN Loop
for step = 1:steps
    dt = deltat(u,v,w,[dx,dy,dz]);
    t = t + dt; disp(num2str(t));
    
    % RK Substep 1 to 3
    [urk1,vrk1,wrk1,prk1,qu1,qv1,qw1] = rk(1,u,v,w,qu0,qv0,qw0);
    [urk2,vrk2,wrk2,prk2,qu2,qv2,qw2] = rk(2,urk1,vrk1,wrk1,qu1,qv1,qw1);
    [urk3,vrk3,wrk3,prk3,qu3,qv3,qw3] = rk(3,urk2,vrk2,wrk2,qu2,qv2,qw2);
    midplane(urk3,['step ' num2str(step) ' t = ' num2str(t)])
    
    % Check convergence
    resu = L2norm(u,urk3);
    res = resu;
    disp(['Main res = ' num2str(res)])
    if res < tol
        disp(['Converged. Main res = ' num2str(res)])
        break
    elseif res > 5000
        disp(['Not converging. Main res = ' num2str(res)])
        break
    else
        % Update
        u = urk3;
        v = vrk3;
        w = wrk3;
        p = prk3;
    end
end