% Pragya Patel
% 17807477
% Thomas Algorithm

function X = thomas(A,B,C,D)
% This function is a linear solver and utilizes
% Thomas Algorithm for Tridiagonal Matrix Eqn
%
% Inputs
%   All the required matrices (A,B,C-LHS) and (D-RHS)
% Output
%   Solution vector (X)

% Tri-diagonal matrix
% a-lower, b-diagonal, c-upper, G-gamma
N = max(size(D));
X = zeros(N,1);
% Downsweep
for k = 2:N
    R = A(k)/B(k-1);
    B(k) = B(k)-R*C(k-1);
    D(k) = D(k)-R*D(k-1);
end
X(N) = D(N)/B(N);
% Upsweep
for i = 1:N-1
    k = N-i;
    X(k) = (D(k)-C(k)*X(k+1))/B(k);
end
end