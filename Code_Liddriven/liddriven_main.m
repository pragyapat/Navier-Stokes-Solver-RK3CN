% ME634 Assignment 4
% 3D Lid Driven Cavity (Laminar)
% Cuboid, Uniform Grid: Periodic in y, lid at zmax
% Non-dimensional solution
% (1/Re) = (nu/UL)

addpath(genpath('mgsolverfiles'))
global Re U L
U = 1; L = 1; Lx = L; Ly = L; Lz = L;
N = 32; Re = 100; Nx = N; Ny = N; Nz = N; % Case 1
% N = 64; Re = 400; Nx = N; Ny = N; Nz = N; % Case 2
% N = 128; Re = 1000; Nx = N; Ny = N; Nz = N; % Case 3
nu = U*L/Re;
nxp2 = Nx+2; nyp2 = Ny+2; nzp2 = Nz+2;
O = zeros(nxp2,nyp2,nzp2);

% Initialize
u = O; v = O; w = O; p = O; t = 0;
% [u, v, w] = updatebc(u, v, w);
u = updatebc(u,'u');
v = updatebc(v,'v');
w = updatebc(w,'w');

%% Grid Coefficient Matrix
global tol CFL dt
C1 = coeff_uni(Nx  ,Ny  ,Nz  ,Lx,Ly,Lz);
C2 = coeff_uni(Nx/2,Ny/2,Nz/2,Lx,Ly,Lz);
C3 = coeff_uni(Nx/4,Ny/4,Nz/4,Lx,Ly,Lz);
C4 = coeff_uni(Nx/8,Ny/8,Nz/8,Lx,Ly,Lz);
tol = 1e-4;
CFL = 1.2;
steps = 100;

%% RK3 CN Loop
for step = 1:steps
    dt = deltat(u,v,w,(C1.s));
    t = t + dt;
    
    % RK Substep 1 to 3
    [urk1,vrk1,wrk1,prk1,qu1,qv1,qw1] ...
        = rk31(u,v,w,p,C1,C2,C3,C4);
    [urk2,vrk2,wrk2,prk2,qu2,qv2,qw2] ...
        = rk32(urk1,vrk1,wrk1,prk1,qu1,qv1,qw1,C1,C2,C3,C4);
    [unew,vnew,wnew,pnew] ...
        = rk33(urk2,vrk2,wrk2,prk2,qu2,qv2,qw2,C1,C2,C3,C4);
    
    % Check convergence
    resu = L2norm(u,unew);
    res = resu;
    disp(['Main res = ' num2str(res)])
    
    if res < tol
        disp(['Converged. Main res = ' num2str(res)])
        break
    elseif res > 100
        disp(['Not converging. Main res = ' num2str(res)])
        break
    else
        % Update
        u = unew;
        v = vnew;
        w = wnew;
        p = pnew;
    end
    % Next step
end

%% Plots
% midplane(u,'U for Re = 100 at ')
% vort = curl(u,v,w);
% midplane(vort,'Vorticity for Re = 100 at ')

% figure
% hold on
% xlabel('u')
% ylabel('z')
% umid(:) = u(nxp2/2,nxp2/2,2:nxp2-1);
% plot(umid, z)
% hold off