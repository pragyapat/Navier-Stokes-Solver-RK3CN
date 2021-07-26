% Pragya Patel
% 17807477
% Multigrid Solver: Restriction

function bc = restrict(rf)
% This function restricts the residual of a finer level
% during the downsweep, to a coarser grid
%
% Inputs
%   residual of the finer level (rf)
% Output
%   restricted matrix for coarser level (bc)

% Size of the finer level
sf = size(rf);
Nxf=sf(1); Nyf=sf(2); Nzf=sf(3);

% Size of the coarser level
Nxc = Nxf/2+1; Nyc = Nyf/2+1; Nzc = Nzf/2+1;
sc = [Nxc,Nyc,Nzc];
bc = zeros(sc);

% Interior
for i = 2:Nxc-1
    for j = 2:Nyc-1
        for k = 2:Nzc-1
            ifi = 2*i-2;
            jfi = 2*j-2;
            kfi = 2*k-2;
            bc(i,j,k) = 0.125*(rf(ifi,jfi,kfi) + rf(ifi+1,jfi,kfi) ...
                + rf(ifi+1,jfi+1,kfi) + rf(ifi+1,jfi+1,kfi+1) ...
                + rf(ifi,jfi+1,kfi+1) + rf(ifi,jfi+1,kfi) ...
                + rf(ifi,jfi,kfi+1) + rf(ifi+1,jfi,kfi+1));
        end
    end
end

% Update boundaries using updatebcr
bc = updatebcr(bc);
end