% Pragya Patel
% 17807477
% Multigrid Solver: Update BC for the problem

function x = updatebcp(x,~)
% This function updates the boundary conditions
% as specified in the problem
%
% Inputs
%   x_old, C
% Output
%   x_new
x = updatebc(x,'p');
end