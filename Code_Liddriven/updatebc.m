% Pragya Patel
% 17807477
% Lid Driven: Update Boundary Conditions

function x = updatebc(x,id)
% This function calculates the rhs of the poisson relation
% at for every interior i,j,k at time step 'n.'
%
% Given
%   Lid velocity (U) = 1, periodic in y-direction
% Inputs
%   Quantity to be updated (x_old) and it's identity (id)
%   id is a string which can take the values 'u','v','w', and 'p'
% Output
%   Updated matrix (x)

global U % U = 1;


% Update
% U velocity
if id == 'u'
    x(1,:,:)        = 0;                % Left wall
    x(end-1,:,:)    = 0;                % Right wall
    x(:,:,1)        = -x(:,:,2);        % Bottom wall
    x(:,:,end)      = 2*U-x(:,:,end-1); % Top lid
    x(:,1,:)        = x(:,end-1,:);     % y
    x(:,end,:)      = x(:,2,:);         % y
    % V velocity
elseif id == 'v'
    x(1,:,:)        = -x(2,:,:);        % Left wall
    x(end,:,:)      = -x(end-1,:,:);    % Right wall
    x(:,:,1)        = -x(:,:,2);        % Bottom wall
    x(:,:,end)      = -x(:,:,end-1);    % Top lid
    x(:,1,:)        = x(:,end-2,:);     % y
    x(:,end-1,:)    = x(:,2,:);         % y
    % W velocity
elseif id == 'w'
    x(1,:,:)        = -x(2,:,:);        % Left wall
    x(end,:,:)      = -x(end-1,:,:);    % Right wall
    x(:,:,1)        = 0;                % Bottom wall
    x(:,:,end-1)    = 0;                % Top lid
    x(:,1,:)        = x(:,end-1,:);     % y
    x(:,end,:)      = x(:,2,:);         % y
    % Pressure Gradient
elseif id == 'p'
    x(1,:,:)        = x(2,:,:);         % Left wall
    x(end,:,:)      = x(end-1,:,:);     % Right wall
    x(:,:,1)        = x(:,:,2);         % Bottom wall
    x(:,:,end)      = x(:,:,end-1);     % Top lid
    x(:,1,:)        = x(:,end-1,:);     % y
    x(:,end,:)      = x(:,2,:);         % y
    x(1,1,1)        = 0;                % Corner
elseif id == 'r'
    x(1,:,:)        = 0;                % Left wall
    x(end,:,:)      = 0;                % Right wall
    x(:,:,1)        = 0;                % Bottom wall
    x(:,:,end)      = 0;                % Top lid
    x(:,1,:)        = 0;                % y
    x(:,end,:)      = 0;                % y
elseif id == 'e'
    x(1,:,:)        = -x(1,:,:);        % Left wall
    x(end,:,:)      = -x(end-1,:,:);    % Right wall
    x(:,:,1)        = -x(:,:,2);        % Bottom wall
    x(:,:,end)      = -x(:,:,end-1);    % Top lid
    x(:,1,:)        = -x(:,2,:);        % y
    x(:,end,:)      = -x(:,end-1,:);    % y
elseif id == '2'
    x(1,:)          = x(1,:);           % Left wall
    x(end,:)        = x(end-1,:);       % Right wall
    x(:,1)          = x(:,2);           % Bottom wall
    x(:,end)        = x(:,end-1);       % Top lid
elseif id == '3'
    x(1,1,1)        = 0;                % Corner
    x(1,:,:)        = x(2,:,:);         % Left wall
    x(end,:,:)      = x(end-1,:,:);     % Right wall
    x(:,:,1)        = x(:,:,2);         % Bottom wall
    x(:,:,end)      = x(:,:,end-1);     % Top lid
    x(:,1,:)        = x(:,end-1,:);     % y
    x(:,end,:)      = x(:,2,:);         % y
    x(1,1,1)        = 0;                % Corner
    x(2,:,:)        = x(1,:,:);         % Left wall
    x(end-1,:,:)    = x(end,:,:);       % Right wall
    x(:,:,2)        = x(:,:,1);         % Bottom wall
    x(:,:,end-1)    = x(:,:,end);       % Top lid
    x(:,1,:)        = x(:,end-1,:);     % y
    x(:,end,:)      = x(:,2,:);         % y
    x(1,1,1)        = 0;                % Corner
    x(1,1,1)        = 0;                % Corner
    x(2,:,:)        = x(1,:,:);         % Left wall
    x(end-1,:,:)    = x(end,:,:);       % Right wall
    x(:,:,2)        = x(:,:,1);         % Bottom wall
    x(:,:,end-1)    = x(:,:,end);       % Top lid
    x(:,1,:)        = x(:,end-1,:);     % y
    x(:,end,:)      = x(:,2,:);         % y
    x(1,1,1)        = 0;                % Corner
end
end

% sz = size(x);
% nxp2 = sz(1); nyp2 = sz(2);
% 
% % Update
% % U velocity
% if id == 'u'
%     nzp2 = sz(3);
%     x(1,:,:)        = 0;                    % Left wall
%     x(nxp2-1,:,:)   = 0;                    % Right wall
%     x(:,:,1)        = -x(:,:,2);            % Bottom wall
%     x(:,:,nzp2)     = 2*U-x(:,:,nzp2-1);    % Top lid
%     x(:,1,:)        = x(:,nyp2-1,:);        % y
%     x(:,nyp2,:)     = x(:,2,:);             % y
%     % V velocity
% elseif id == 'v'
%     nzp2 = sz(3);
%     x(1,:,:)        = -x(2,:,:);            % Left wall
%     x(nxp2,:,:)     = -x(nxp2-1,:,:);       % Right wall
%     x(:,:,1)        = -x(:,:,2);            % Bottom wall
%     x(:,:,nzp2)     = -x(:,:,nzp2-1);       % Top lid
%     x(:,1,:)        = x(:,nyp2-2,:);        % y
%     x(:,nyp2-1,:)   = x(:,2,:);             % y
%     % W velocity
% elseif id == 'w'
%     nzp2 = sz(3);
%     x(1,:,:)        = -x(2,:,:);            % Left wall
%     x(nxp2,:,:)     = -x(nxp2-1,:,:);       % Right wall
%     x(:,:,1)        = 0;                    % Bottom wall
%     x(:,:,nzp2-1)   = 0;                    % Top lid
%     x(:,1,:)        = x(:,nyp2-1,:);        % y
%     x(:,nyp2,:)     = x(:,2,:);             % y
%     % Pressure Gradient
% elseif id == 'p'
%     nzp2 = sz(3);
%     x(1,1,1)        = 0;                    % Corner
%     x(1,:,:)        = x(2,:,:);             % Left wall
%     x(nxp2,:,:)     = x(nxp2-1,:,:);        % Right wall
%     x(:,:,1)        = x(:,:,2);             % Bottom wall
%     x(:,:,nzp2)     = x(:,:,nzp2-1);        % Top lid
%     x(:,1,:)        = x(:,nyp2-1,:);        % y
%     x(:,nyp2,:)     = x(:,2,:);             % y
%     x(1,1,1)        = 0;                    % Corner
% elseif id == 'r'
%     nzp2 = sz(3);
%     x(1,:,:)        = 0;                    % Left wall
%     x(nxp2,:,:)     = 0;                    % Right wall
%     x(:,:,1)        = 0;                    % Bottom wall
%     x(:,:,nzp2)     = 0;                    % Top lid
%     x(:,1,:)        = 0;                    % y
%     x(:,nyp2,:)     = 0;                    % y
% elseif id == 'e'
%     nzp2 = sz(3);
%     x(1,:,:)        = -x(1,:,:);            % Left wall
%     x(nxp2,:,:)     = -x(nxp2-1,:,:);       % Right wall
%     x(:,:,1)        = -x(:,:,2);            % Bottom wall
%     x(:,:,nzp2)     = -x(:,:,nzp2-1);       % Top lid
%     x(:,1,:)        = -x(:,2,:);            % y
%     x(:,nyp2,:)     = -x(:,nyp2-1,:);       % y
% elseif id == '2'
%     nzp2 = sz(2);
%     x(1,:)          = x(1,:);               % Left wall
%     x(nxp2,:)       = x(nxp2-1,:);          % Right wall
%     x(:,1)          = x(:,2);               % Bottom wall
%     x(:,nzp2)       = x(:,nzp2-1);          % Top lid
% end

% x1 = x0;
% 
% % Update
% % U velocity
% if id == 'u'
%     nzp2 = sz(3);
%     x1(1,:,:)        = 0;                    % Left wall
%     x1(nxp2-1,:,:)   = 0;                    % Right wall
%     x1(:,:,1)        = -x0(:,:,2);            % Bottom wall
%     x1(:,:,nzp2)     = 2*U-x0(:,:,nzp2-1);    % Top lid
%     x1(:,1,:)        = x0(:,nyp2-1,:);        % y
%     x1(:,nyp2,:)     = x0(:,2,:);             % y
%     % V velocity
% elseif id == 'v'
%     nzp2 = sz(3);
%     x1(1,:,:)        = -x0(2,:,:);            % Left wall
%     x1(nxp2,:,:)     = -x0(nxp2-1,:,:);       % Right wall
%     x1(:,:,1)        = -x0(:,:,2);            % Bottom wall
%     x1(:,:,nzp2)     = -x0(:,:,nzp2-1);       % Top lid
%     x1(:,1,:)        = x0(:,nyp2-2,:);        % y
%     x1(:,nyp2-1,:)   = x0(:,2,:);             % y
%     % W velocity
% elseif id == 'w'
%     nzp2 = sz(3);
%     x1(1,:,:)        = -x0(2,:,:);            % Left wall
%     x1(nxp2,:,:)     = -x0(nxp2-1,:,:);       % Right wall
%     x1(:,:,1)        = 0;                    % Bottom wall
%     x1(:,:,nzp2-1)   = 0;                    % Top lid
%     x1(:,1,:)        = x0(:,nyp2-1,:);        % y
%     x1(:,nyp2,:)     = x0(:,2,:);             % y
%     % Pressure Gradient
% elseif id == 'p'
%     nzp2 = sz(3);
%     x1(1,1,1)        = 0;                    % Corner
%     x1(1,:,:)        = x0(2,:,:);             % Left wall
%     x1(nxp2,:,:)     = x0(nxp2-1,:,:);        % Right wall
%     x1(:,:,1)        = x0(:,:,2);             % Bottom wall
%     x1(:,:,nzp2)     = x0(:,:,nzp2-1);        % Top lid
%     x1(:,1,:)        = x0(:,nyp2-1,:);        % y
%     x1(:,nyp2,:)     = x0(:,2,:);             % y
%     x1(1,1,1)        = 0;                    % Corner
% elseif id == 'r'
%     nzp2 = sz(3);
%     x1(1,:,:)        = 0;                    % Left wall
%     x1(nxp2,:,:)     = 0;                    % Right wall
%     x1(:,:,1)        = 0;                    % Bottom wall
%     x1(:,:,nzp2)     = 0;                    % Top lid
%     x1(:,1,:)        = 0;                    % y
%     x1(:,nyp2,:)     = 0;                    % y
% elseif id == 'e'
%     nzp2 = sz(3);
%     x1(1,:,:)        = -x0(1,:,:);            % Left wall
%     x1(nxp2,:,:)     = -x0(nxp2-1,:,:);       % Right wall
%     x1(:,:,1)        = -x0(:,:,2);            % Bottom wall
%     x1(:,:,nzp2)     = -x0(:,:,nzp2-1);       % Top lid
%     x1(:,1,:)        = -x0(:,2,:);            % y
%     x1(:,nyp2,:)     = -x0(:,nyp2-1,:);       % y
% elseif id == '2'
%     nzp2 = sz(2);
%     x1(1,:)          = x0(1,:);               % Left wall
%     x1(nxp2,:)       = x0(nxp2-1,:);          % Right wall
%     x1(:,1)          = x0(:,2);               % Bottom wall
%     x1(:,nzp2)       = x0(:,nzp2-1);          % Top lid
% end


% % % Pragya Patel
% % % 17807477
% % % Lid Driven: Update Boundary Conditions
% %
% % function x = updatebc(x_old,id)
% % % This function calculates the rhs of the poisson relation
% % % at for every interior i,j,k at time step 'n.'
% % %
% % % Given
% % %   Lid velocity (U) = 1, periodic in y-direction
% % % Inputs
% % %   Quantity to be updated (x_old) and it's identity (id)
% % %   id is a string which can take the values 'u','v','w', and 'p'
% % % Output
% % %   Updated matrix (x)
% %
% % U = 1;
% %
% % % Initialize
% % sz = size(x_old);
% % nxp2 = sz(1); nyp2 = sz(2); nzp2 = sz(3);
% % x = x_old;
% %
% % % Update
% % for k = 2:nzp2-1
% %     for j = 2:nyp2-1
% %         for i = 2:nxp2-1
% %             % U velocity
% %             if id == 'u'
% %                 x(1,j,k)        = 0;
% %                 x(nxp2-1,j,k)   = 0;
% %                 x(i,1,k)        = x_old(i,nyp2-1,k);
% %                 x(i,nyp2,k)     = x_old(i,2,k);
% %                 x(i,j,1)        = -x_old(i,j,2);
% %                 x(i,j,nzp2)     = 2*U-x_old(i,j,nzp2-1);
% %                 % V velocity
% %             elseif id == 'v'
% %                 x(1,j,k)        = -x_old(2,j,k);
% %                 x(nxp2,j,k)     = -x_old(nxp2-1,j,k);
% %                 x(i,1,k)        = x_old(i,nyp2-2,k);
% %                 x(i,nyp2-1,k)   = x_old(i,2,k);
% %                 x(i,j,1)        = -x_old(i,j,2);
% %                 x(i,j,nzp2)     = -x_old(i,j,nzp2-1);
% %                 % W velocity
% %             elseif id == 'w'
% %                 x(1,j,k)        = -x_old(2,j,k);
% %                 x(nxp2,j,k)     = -x_old(nxp2-1,j,k);
% %                 x(i,1,k)        = x_old(i,nyp2-1,k);
% %                 x(i,nyp2,k)     = x_old(i,2,k);
% %                 x(i,j,2)        = 0;
% %                 x(i,j,nzp2-1)   = 0;
% %                 % Pressure Gradient
% %             elseif id == 'p'
% %                 x(1,j,k)        = x_old(2,j,k);
% %                 x(nxp2,j,k)     = x_old(nxp2-1,j,k);
% %                 x(i,1,k)        = x_old(i,nyp2-1,k);
% %                 x(i,nyp2,k)     = x_old(i,2,k);
% %                 x(i,j,1)        = x_old(i,j,2);
% %                 x(i,j,nzp2)     = x_old(i,j,nzp2-1);
% %             end
% %         end
% %     end
% % end
% % end
% 
% 
% % function x = updatebc(x,id)
% % sz = size(x);
% % nxp2 = sz(1); nyp2 = sz(2); nzp2 = sz(3);
% % for j = 2:nyp2-1
% %     for k = 2:nzp2-1
% %         if id == 'u'
% %             x(1,j,k) = 0;
% %             x(nxp2-1,j,k) = 0;
% %         elseif id == 'v'
% %             x(1,j,k) = -x(2,j,k);
% %             x(nxp2,j,k) = -x(nxp2-1,j,k);
% %         elseif id == 'w'
% %             x(1,j,k) = -x(2,j,k);
% %             x(nxp2,j,k) = -x(nxp2-1,j,k);
% %         end
% %     end
% % end
% % for i = 2:nxp2-1
% %     for j = 2:nyp2-1
% %         if id == 'u'
% %             x(i,j,1) = -x(i,j,2);
% %             x(i,j,nzp2) = 2 - x(i,j,nzp2-1);
% %         elseif id == 'v'
% %             x(i,j,1) = -x(i,j,2);
% %             x(i,j,nzp2) = -x(i,j,nzp2-1);
% %         elseif id == 'w'
% %             x(i,j,1) = 0;
% %             x(i,j,nzp2-1) = 0;
% %         end
% %     end
% % end
% % for i = 2:nxp2-1
% %     for k = 2:nzp2-1
% %         if id == 'u'
% %             x(i,1,k) = x(i,nyp2-1,k);
% %             x(i,nyp2,k) = x(i,2,k);
% %         elseif id == 'v'
% %             x(i,1,k) = x(i,nyp2-2,k);
% %             x(i,nyp2-1,k) = x(i,2,k);
% %         elseif id == 'w'
% %             x(i,1,k) = x(i,nyp2-1,k);
% %             x(i,nyp2,k) = x(i,2,k);
% %         end
% %     end
% % end
% %
% % if id == 'p'
% %     for j = 2:nyp2-1
% %         for k = 2:nzp2-1
% %             x(1,j,k) = x(2,j,k);
% %             x(nxp2,j,k) = x(nxp2-1,j,k);
% %         end
% %     end
% %
% %     for i = 2:nxp2 -1
% %         for j = 2:nyp2 -1
% %             x(i,j,1) = x(i,j,2);
% %             x(i,j,nzp2) = x(i,j,nzp2-1);
% %         end
% %     end
% %
% %     for i = 2:nxp2-1
% %         for k = 2:nzp2-1
% %             x(i,1,k) = x(i,nyp2-1,k);
% %             x(i,nyp2,k) = x(i,2,k);
% %         end
% %     end
% % end
% % end
% 
% % DUMP
% % for k = 1:nzp2
% %     for j = 1:nyp2
% %         for i = 1:nxp2
% %             % U velocity
% %             if id == 'u'
% %                 if i == 1
% %                     x(i,j,k) = 0;
% %                 elseif i == nxp2-1
% %                     x(i,j,k) = 0;
% %                 elseif j == 1
% %                     x(i,j,k) = x_old(i,nyp2-1,k);
% %                 elseif j == nyp2
% %                     x(i,j,k) = x_old(i,2,k);
% %                 elseif k == 1
% %                     x(i,j,k) = -x_old(i,j,k+1);
% %                 elseif k == nzp2
% %                     x(i,j,k) = 2*U - x_old(i,j,k-1);
% %                 end
% %                 % V velocity
% %             elseif id == 'v'
% %                 if i == 1
% %                     x(i,j,k) = -x_old(i+1,j,k);
% %                 elseif i == nxp2
% %                     x(i,j,k) = -x_old(i-1,j,k);
% %                 elseif j == 1
% %                     x(i,j,k) = x_old(i,nyp2-2,k);
% %                 elseif j == nyp2-1
% %                     x(i,j,k) = x_old(i,2,k);
% %                 elseif k == 1
% %                     x(i,j,k) = -x_old(i,j,k+1);
% %                 elseif k == nzp2
% %                     x(i,j,k) = -x_old(i,j,k-1);
% %                 end
% %                 % W velocity
% %             elseif id == 'w'
% %                 if i == 1
% %                     x(i,j,k) = -x_old(i+1,j,k);
% %                 elseif i == nxp2
% %                     x(i,j,k) = -x_old(i-1,j,k);
% %                 elseif j == 1
% %                     x(i,j,k) = x_old(i,nyp2-1,k);
% %                 elseif j == nyp2
% %                     x(i,j,k) = x_old(i,2,k);
% %                 elseif k == 2
% %                     x(i,j,k) = 0;
% %                 elseif k == nzp2-1
% %                     x(i,j,k) = 0;
% %                 end
% %                 % Pressure Gradient Boundary Conditions
% %             elseif id == 'p'
% %                 if i == 1
% %                     x(i,j,k) = x_old(i+1,j,k);
% %                 elseif i == nxp2
% %                     x(i,j,k) = x_old(i-1,j,k);
% %                 elseif j == 1
% %                     x(i,j,k) = x_old(i,nyp2-1,k);
% %                 elseif j == nyp2
% %                     x(i,j,k) = x_old(i,2,k);
% %                 elseif k == 1
% %                     x(i,j,k) = x_old(i,j,k+1);
% %                 elseif k == nzp2
% %                     x(i,j,k) = x_old(i,j,k-1);
% %                 end
% %             end
% %         end
% %     end
% % end
% 
% % function [u, v, w] = updatebc(u, v, w)
% % sz = size(u);
% % nxp2 = sz(1); nyp2 = sz(2); nzp2 = sz(3);
% % for j = 2:nyp2-1
% %     for k = 2:nzp2-1
% %         u(1,j,k) = 0;
% %         u(nxp2-1,j,k) = 0;
% %         v(1,j,k) = -v(2,j,k);
% %         v(nxp2,j,k) = -v(nxp2-1,j,k);
% %         w(1,j,k) = -w(2,j,k);
% %         w(nxp2,j,k) = -w(nxp2-1,j,k);
% %     end
% % end
% %
% % for i = 2:nxp2-1
% %     for j = 2:nyp2-1
% %         u(i,j,1) = -u(i,j,2);
% %         u(i,j,nzp2) = 2 - u(i,j,nzp2-1);
% %         v(i,j,1) = -v(i,j,2);
% %         v(i,j,nzp2) = -v(i,j,nzp2-1);
% %         w(i,j,1) = 0;
% %         w(i,j,nzp2-1) = 0;
% %     end
% % end
% %
% % for i = 2:nxp2-1
% %     for k = 2:nzp2-1
% %         u(i,1,k) = u(i,nyp2-1,k);
% %         u(i,nyp2,k) = u(i,2,k);
% %         v(i,1,k) = v(i,nyp2-2,k);
% %         v(i,nyp2-1,k) = v(i,2,k);
% %         w(i,1,k) = w(i,nyp2-1,k);
% %         w(i,nyp2,k) = w(i,2,k);
% %     end
% % end
% % end