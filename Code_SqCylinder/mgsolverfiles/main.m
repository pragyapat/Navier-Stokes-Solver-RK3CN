% Pragya Patel
% 17807477
% Multigrid Solver: Main
% Functions used:
%   coeff_uni.m, pderhs.m, GSp.m, V4.m, L2norm.m
%   and other functions used within these

%% Initialize
global tol 
tol = 1.0e-6;
Lx = 1; Ly = 1; Lz = 1;

% Grid size
N = 32; 
Nx1 = N; Ny1 = N; Nz1 = N;

%
nxp21 = Nx1+2; nyp21 = Ny1+2; nzp21 = Nz1+2;
Nx2 = Nx1/2; Ny2 = Ny1/2; Nz2 = Nz1/2;
Nx3 = Nx2/2; Ny3 = Ny2/2; Nz3 = Nz2/2;
Nx4 = Nx3/2; Ny4 = Ny3/2; Nz4 = Nz3/2;
phi1 = zeros(nxp21,nyp21,nzp21);

%% Coefficient matrices
% global C1 C2 C3 C4
C1 = coeff_uni(Nx1,Ny1,Nz1,Lx,Ly,Lz);
C2 = coeff_uni(Nx2,Ny2,Nz2,Lx,Ly,Lz);
C3 = coeff_uni(Nx3,Ny3,Nz3,Lx,Ly,Lz);
C4 = coeff_uni(Nx4,Ny4,Nz4,Lx,Ly,Lz);
g = pderhs([Nx1,Ny1,Nz1,Lx,Ly,Lz]);

p1 = mgsolver(g,C1,C2,C3,C4);
p2 = mgsolver2(g,C1,C2,C3,C4);

% %% Execution
% % Level 1
% phi1 = GSp(phi1,C1,g,10);
% phi0 = phi1;
% 
% % V-cycle begins
% tmax = 300;
% for t = 1:tmax
% %     phi1 = V3(C1,C2,C3,g,phi0,t);
%     phi1 = V4(C1,C2,C3,C4,g,phi0,t);
%     rms = L2norm(phi0,phi1);
%     if rms < tol
%         disp(['Number of iterations ', num2str(t)])
%         disp(['res = ', num2str(rms)])
%         break
%     else
%         phi0 = phi1;
%     end
% end

% midplane(p1,'Level 1')
% colorbar