% Pragya Patel
% 17807477
% Solving for intermediate velocity u*,v*,w*

function x = velstar(a,b,c,rhscolumn,s)
% This function is a linear solver and utilizes
% Thomas Algorithm for Tridiagonal Matrix Eqn
%
% Inputs
%   TD matrix coefficients (a,b,c) and the column RHS (rhscolumn)
%   s string for specifying u, v or w star (for BCs)
% Output
%   Solution vector (X)

% COEFFICIENT MATRICES
% Tri-diagonal matrix
% a-lower, b-diagonal, c-upper, G-gamma
sz = size(rhscolumn);
nzp2 = sz(3);
A = zeros(nzp2,1);
B = zeros(nzp2,1);
C = zeros(nzp2,1);
D = zeros(nzp2,1);

% Coeffs for the interior nodes
for k = 2:nzp2-1
    A(k) = a; B(k) = b; C(k) = c;
    D(k) = rhscolumn(k);
end

% Coeffs according to BCs
% 1 = Bottom wall, nzp2 = Top lid
if s == 'u' % for ustar
    A(1) = 0; A(nzp2) = 1;
    B(1) = 1; B(nzp2) = 1;
    C(1) = 1; C(nzp2) = 0;
    D(1) = 0; D(nzp2) = 2;
elseif s == 'v' % for vstar
    A(1) = 0; A(nzp2) = 1;
    B(1) = 1; B(nzp2) = 1;
    C(1) = 1; C(nzp2) = 0;
    D(1) = 0; D(nzp2) = 0;
elseif s == 'w' % for wstar
    A(1) = 0; A(nzp2-1) = 0; A(nzp2) = [];
    B(1) = 1; B(nzp2-1) = 1; B(nzp2) = [];
    C(1) = 1; C(nzp2-1) = 0; C(nzp2) = [];
    D(1) = 0; D(nzp2-1) = 0; D(nzp2) = [];
end

x = thomas(A,B,C,D); % for one column
end