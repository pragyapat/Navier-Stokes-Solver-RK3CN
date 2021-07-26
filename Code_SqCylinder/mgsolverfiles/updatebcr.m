% Pragya Patel
% 17807477
% Multigrid Solver: Update boundary conditions for the residual

function res = updatebcr(res)
% This function updates the boundary conditions
% of the residuals/errors
%
% Inputs
%   residual (res)
% Output
%   updated variable (res)
res = updatebc(res,'r');
end

% s = size(res);
% nxp2 = s(1); nyp2 = s(2); nzp2 = s(3);
% % res_old = res;
% for i  = 2:nxp2-1
%     for j = 2:nyp2-1
%         res(i,j,1) = 0;%-res(i,j,2);
%         res(i,j,nzp2) =0;% -res(i,j,nzp2-1);
%     end
% end
% 
% for i = 2:nxp2-1
%     for k = 2:nzp2-1
%         res(i,1,k) = 0;%-res(i,2,k);
%         res(i,nyp2,k) =0;% -res(i,nyp2-1,k);
%     end
% end
% 
% for j = 2:nyp2-1
%     for k = 2:nzp2-1
%         res(1,j,k) =0;% -res(2,j,k);
%         res(nxp2,j,k) =0;% -res(nxp2-1,j,k);
%     end
% end

% for k = 1:Nzp
%     for j = 1:Nyp
%         for i = 1:Nxp
%             % Boundary
%             if i == 1
%                 res(i,j,k) = -res_old(i+1,j,k);
%             elseif i == Nxp
%                 res(i,j,k) = -res_old(i-1,j,k);
%             elseif j == 1
%                 res(i,j,k) = -res_old(i,j+1,k);
%             elseif j == Nyp
%                 res(i,j,k) = -res_old(i,j-1,k);
%             elseif k == 1
%                 res(i,j,k) = -res_old(i,j,k+1);
%             elseif k == Nzp
%                 res(i,j,k) = -res_old(i,j,k-1);
%             end
%         end
%     end
% end

% DUMP
%
% boundary
% res(1,:,:)    = -res_old(2,:,:);
% res(Nxp,:,:)  = -res_old(Nxp-1,:,:);
% res(:,1,:)    =  res_old(:,Nyp-1,:);
% res(:,Nyp,:)  =  res_old(:,2,:);
% res(:,:,1)    = -res_old(:,:,2);
% res(:,:,Nzp)  = -res_old(:,:,Nzp-1);