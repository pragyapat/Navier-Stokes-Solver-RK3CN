% Pragya Patel
% 17807477
% Multigrid Solver: Prolongation

function xf = prolong(ec)
% This function prolongs the error of a coarser level
% to a finer grid during the upsweep
%
% Inputs
%   error of coarser grid (ec)
% Output
%   initial guess for the finer GS (xf)

% Size of the coarser level
sc = size(ec);
Nxc=sc(1); Nyc=sc(2); Nzc=sc(3);

% Size of the finer level
Nxf = 2*(Nxc-1); Nyf = 2*(Nyc-1); Nzf = 2*(Nzc-1);
sf = [Nxf,Nyf,Nzf];
xf = zeros(sf);

for k = 2:Nzf-1
    for j = 2:Nyf-1
        for i = 2:Nxf-1
            icc = ceil(i/2) + (1 - rem(i,2)); % closest coarse
            jcc = ceil(j/2) + (1 - rem(j,2));
            kcc = ceil(k/2) + (1 - rem(k,2));
            icf = icc - 1 + 2*rem(i,2);
            jcf = jcc - 1 + 2*rem(j,2);
            kcf = kcc - 1 + 2*rem(k,2);
            xf(i,j,k) = 0.015625*(27*ec(icc,jcc,kcc) ...
                + 9*ec(icf,jcc,kcc) + 9*ec(icc,jcf,kcc) + 9*ec(icc,jcc,kcf) ...
                + 3*ec(icf,jcf,kcc) + 3*ec(icf,jcc,kcf) + 3*ec(icc,jcf,kcf) ...
                + ec(icf,jcf,kcf));
        end
    end
end

% Update boundaries using updatebcr
end