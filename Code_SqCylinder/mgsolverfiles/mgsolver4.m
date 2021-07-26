% Pragya Patel
% 17807477
% MG Solver - 4 levels in V cycle

function p1 = mgsolver4(rhsp,C1,C2,C3,C4)
% This function solves the pressure poisson eqn
%
% Inputs
%   RHS for the equation (rhsp
% Output
%   Updated matrix (p1)

% Initialize
global tol
p1 = zeros(size(rhsp));

% Level 1
p1 = GSp(p1,C1,rhsp,10);
p0 = p1;

% V-cycles begin
tmax = 400;
for t = 1:tmax
    p1 = V4(C1,C2,C3,C4,rhsp,p0,t);
    r = L2norm(p0,p1);
    if r < tol
        disp(['MG4 converging at res = ' num2str(r)])
        break
    elseif r > 20
        disp(['MG4 diverging at res = ' num2str(r)])
        break
    else
        p0 = p1;
    end
end
if t == tmax
    disp(['Insufficient MG4 tmax ', num2str(t) ' with res = ', num2str(r)])
end
end