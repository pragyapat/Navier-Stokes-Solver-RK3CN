% N = 100;
% a = -1;
% b = 2;
% c = -1;
% D = rand(N,1);
% 
% P = diag(b*ones(1,N)) + diag(c*ones(1,N-1),1) + diag(a*ones(1,N-1),-1);
% tic,x = P\D;toc
% 
% A = a*ones(N,1);
% B = b*ones(N,1);
% C = c*ones(N,1);
% tic,xmy = thomas(A,B,C,D);toc

% % Init
%     G(1) = C(1)/B(1);
%     R(1) = D(1)/B(1);
% % Down
%         G(k) = C(k)/(B(k)-A(k)*G(k-1));
%         R(k) = (D(k)-R(k-1)*A(k))/(B(k)-G(k-1)*A(k));
% % Up
%         X(k) = R(k)-G(k)*X(k+1);

P = zeros(5,5,5);
for i = 1:5
    for j = 1:5
        for k = 1:5
            P(i,j,k) = i*100+j*10+k;
        end
    end
end

Pbc = updatebc(P,'p');
Pbc3 = updatebc(Pbc,'3');
