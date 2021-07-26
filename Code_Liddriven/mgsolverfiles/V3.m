% Pragya Patel
% 17807477
% Multigrid Solver: V-cycle upto the 3rd level
% Functions used:
%   resi.m, updatebcr.m, restrict.m,
%   updatebcp.m, GSp.m and prolong.m

function p1 = V3(C1,C2,C3,rhs1,p0,iter)
% This function represents the multigrid accelerator
% upto 3-levels of the V-cycle (Total GS within = 10+10+4+8+10 = 42)
%
% Inputs
%   finest coefficient matrix (C1), level 2 coefficient matrix (C2),
%   level 3 coefficient matrix (C3), rhs
% Output
%   p1(final)

s = C1.s; Nx1=s(1); Ny1=s(2); Nz1=s(3);

% V-cycle begins
p1 = GSp(p0,C1,rhs1,10);

% DOWNSWEEP
% 1 to 2
res1 = resi(p1,C1,rhs1);
res1 = updatebcr(res1);
rhs2 = restrict(res1);
rhs2 = updatebcr(rhs2);

% Level 2
phi2 = zeros(Nx1/2+2,Ny1/2+2,Nz1/2+2);
phi2 = GSr(phi2,C2,rhs2,10);

% 2 to 3
res2 = resi(phi2,C2,rhs2);
res2 = updatebcr(res2);
rhs3 = restrict(res2);
rhs3 = updatebcr(rhs3);

% Level 3
phi3 = zeros(Nx1/4+2,Ny1/4+2,Nz1/4+2);
phi3 = GSr(phi3,C3,rhs3,4);

% UPSWEEP
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