% Pragya Patel
% 17807477
% Multigrid Solver: V-cycle upto the 4th level
% Functions used:
%   resi.m, updatebcr.m, restrict.m,
%   updatebcp.m, GSp.m and prolong.m

function p1 = V4(C1,C2,C3,C4,rhs1,p0,iter)
% This function represents the multigrid accelerator
% upto 4-levels of the V-cycle
%
% Inputs
%   finest coefficient matrix (C1), level 2,3,4 coefficient matrices (C2,
%   C2,C3,C4), rhs (g)
% Output
%   p1 (final)

s = C1.s; Nx1=s(1); Ny1=s(2); Nz1=s(3);

% V-cycle begins
p1 = GSp(p0,C1,rhs1,10);

% DOWNSWEEP
% 1 to 2
res = resi(p1,C1,rhs1);
res = updatebcr(res);
rhs2 = restrict(res);
rhs2 = updatebcr(rhs2);
% Level 2
phi2 = zeros(Nx1/2+2,Ny1/2+2,Nz1/2+2);
phi2 = GSr(phi2,C2,rhs2,10);
% 2 to 3
res = resi(phi2,C2,rhs2);
res = updatebcr(res);
rhs3 = restrict(res);
rhs3 = updatebcr(rhs3);
% Level 3
phi3 = zeros(Nx1/4+2,Ny1/4+2,Nz1/4+2);
phi3 = GSr(phi3,C3,rhs3,4);
% 3 to 4
res = resi(phi3,C3,rhs2);
res = updatebcr(res);
rhs4 = restrict(res);
rhs4 = updatebcr(rhs4);
% Level 4
phi4 = zeros(Nx1/8+2,Ny1/8+2,Nz1/8+2);
phi4 = GSr(phi4,C4,rhs4,4);

% UPSWEEP
% 4 to 3
e = prolong(phi4);
e = updatebcr(e);
phi3 = phi3 + e;
phi3 = GSr(phi3,C3,rhs3,4);
% 3 to 2
e = prolong(phi3);
e = updatebcr(e);
phi2 = phi2 + e;
phi2 = GSr(phi2,C2,rhs2,8);
% 2 to 1
e = prolong(phi2);
e = updatebcr(e);
p1 = p1 + e;
p1 = GSp(p1,C1,rhs1,10);
end