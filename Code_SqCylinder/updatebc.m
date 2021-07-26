% Pragya Patel (17807477)
% Flow past a Square Cylinder
% BOUNDARY CONDITIONS

function x = updatebc(sqn,x,id)
% Update boundary conditions for this flow
%
% Given
%   Inlet velocity (U) = 1, periodic in y-direction
% Inputs
%   Quantity to be updated (x) and its identity (id)
%   id is a string which can take the values 'u','v','w', and 'p'
% Output
%   Updated matrix (x)
% global U

x1 = sqn(1); x2 = sqn(2);
z1 = sqn(3); z2 = sqn(4);

% Update
% U velocity
if id == 'u'
    % Square
    x(x1,:,z1:z2)           = -x(x1-1,:,z1:z2);     % Left face
    x(x2-1,:,z1+1:z2-1)     = -x(x2,:,z1+1:z2-1);   % Right face
    x(x1+1:x2-2,:,z1+1:z2-1)= 0;                    % Interior
    x(x1+1:x2-2,:,z1)       = 0;                    % Bottom face
    x(x1+1:x2-2,:,z2)       = 0;                    % Top face
    % Domain
    x(1,:,:)        = 1; % 2*U-x(2,:,:);                 % Inflow
    x(end,:,:)      = x(end-1,:,:);                 % Outflow
    x(:,:,1)        = x(:,:,2);                     % Bottom
    x(:,:,end)      = x(:,:,end-1);                 % Top
    x(:,1,:)        = x(:,end-1,:);                 % y
    x(:,end,:)      = x(:,2,:);                     % y
    
    % V velocity
elseif id == 'v'
    % Square
    x(x1:x2,:,z1:z2)= 0;                            % All
    % Domain
    x(1,:,:)        = 0;                            % Inflow
    x(end,:,:)      = x(end-1,:,:);                 % Outflow
    x(:,:,1)        = x(:,:,2);                     % Bottom
    x(:,:,end)      = x(:,:,end-1);                 % Top
    x(:,1,:)        = x(:,end-2,:);                 % y
    x(:,end-1,:)    = x(:,2,:);                     % y
    
    % W velocity
elseif id == 'w'
    % Square
    x(x1,:,z1:z2-2)         = 0;                    % Left face
    x(x2,:,z1:z2-2)         = 0;                    % Right face
    x(x1+1:x2-1,:,z1:z2-2)  = 0;                    % Interior
    x(x1:x2-2,:,z1)         = -x(x1:x2-2,:,z1-1);   % Bottom face
    x(x1:x2-2,:,z2-1)       = -x(x1:x2-2,:,z2);     % Top face
    % Domain
    x(1,:,:)        = 0;                            % Inflow
    x(end,:,:)      = x(end-1,:,:);                 % Outflow
    x(:,:,1)        = x(:,:,2);                     % Bottom
    x(:,:,end-1)    = x(:,:,end-2);                 % Top
    x(:,1,:)        = x(:,end-1,:);                 % y
    x(:,end,:)      = x(:,2,:);                     % y
    
    % Pressure
elseif id == 'p'
    % Square
    x(x1,:,z1:z2)           = x(x1-1,:,z1:z2);      % Square left
    x(x1+1,:,z1:z2)         = x(x1-1,:,z1:z2);      % Square left
    x(x1+1:x2-1,:,z1)       = x(x1+1:x2-1,:,z1-1);  % Square bottom
    x(x1+1:x2-1,:,z2+1)     = x(x1+1:x2-1,:,z1-1);  % Square bottom
    x(x1+2:x2-2,:,z1+1:z2-1)= 0;                    % Square interior
    x(x2-1,:,z1:z2)         = x(x2+1,:,z1:z2);      % Square right
    x(x2,:,z1:z2)           = x(x2+1,:,z1:z2);      % Square right
    x(x1+1:x2-1,:,z2-1)     = x(x1+1:x2-1,:,z2+1);  % Square top
    x(x1+1:x2-1,:,z2)       = x(x1+1:x2-1,:,z2+1);  % Square top
    % Domain
    x(1,:,:)        = x(2,:,:);                     % Inflow
    x(end,:,:)      = 0;                            % Outflow
    x(:,:,1)        = 0;                            % Bottom
    x(:,:,end)      = 0;                            % Top
    x(:,1,:)        = x(:,end-1,:);                 % y
    x(:,end,:)      = x(:,2,:);                     % y
    x(1,1,1)        = 0;                            % Corner
end


% nblx = sqn(1,1); nblz = sqn(1,2);
% nbrx = sqn(2,1); nbrz = sqn(2,2);
% ntlx = sqn(3,1); ntlz = sqn(3,2);
% ntrx = sqn(4,1); ntrz = sqn(4,2);

% u
%     x(40,:,36:44)   = -x(39,:,36:44);   % Square left face
%     x(41:46,:,36)   = 0;                % Square bottom face
%     x(41:46,:,36:44)= 0;                % Square interior
%     x(47,:,36:44)   = -x(48,:,36:44);   % Square right face
%     x(41:46,:,44)   = 0;                % Square top face
% v
%     x(40,:,36:44)   = 0;                % Square left
%     x(40:47,:,36)   = 0;                % Square bottom
%     x(40:48,:,36:44)= 0;                % Square interior
%     x(48,:,36:44)   = 0;                % Square right
%     x(40:47,:,44)   = 0;                % Square top
% w
%     x(40,:,36:43)   = 0;                % Square left
%     x(40:48,:,36)   = -x(40:48,:,35);   % Square bottom
%     x(40:48,:,37:42)= 0;                % Square interior
%     x(48,:,36:43)   = 0;                % Square right
%     x(40:48,:,43)   = -x(40:48,:,44);   % Square top
% p
%     x(40,:,36:44)   = x(39,:,36:44);    % Square left
%     x(41,:,36:44)   = x(39,:,36:44);    % Square left
%     x(41:47,:,36)   = x(41:47,:,35);    % Square bottom
%     x(41:47,:,37)   = x(41:47,:,35);    % Square bottom
%     x(42:46,:,37:43)= 0;                % Square interior
%     x(47,:,36:44)   = x(49,:,36:44);    % Square right
%     x(48,:,36:44)   = x(49,:,36:44);    % Square right
%     x(41:47,:,43)   = x(41:47,:,45);    % Square top
%     x(41:47,:,44)   = x(41:47,:,45);    % Square top
% if MGSOLVE
% elseif id == 'r'
%     x(1,:,:)        = 0;                            % Inflow
%     x(end,:,:)      = 0;                            % Outflow
%     x(:,:,1)        = 0;                            % Bottom
%     x(:,:,end)      = 0;                            % Top
%     x(:,1,:)        = 0;                            % y
%     x(:,end,:)      = 0;                            % y
% elseif id == 'e'
%     x(1,:,:)        = -x(1,:,:);                    % Inflow
%     x(end,:,:)      = -x(end-1,:,:);                % Outflow
%     x(:,:,1)        = -x(:,:,2);                    % Bottom
%     x(:,:,end)      = -x(:,:,end-1);                % Top
%     x(:,1,:)        = -x(:,2,:);                    % y
%     x(:,end,:)      = -x(:,end-1,:);                % y