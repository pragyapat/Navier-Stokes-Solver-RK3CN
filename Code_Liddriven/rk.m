function [urk1,vrk1,wrk1,prk1,qu1,qv1,qw1] = rk(sub,urk0,vrk0,wrk0,qu0,qv0,qw0,C1,C2,C3,C4)
sz = size(urk0);
nxp2 = sz(1); nyp2 = sz(2); nzp2 = sz(3);
dx = 1/(nxp2-2);
dy = 1/(nyp2-2);
dz = 1/(nzp2-2); d = [dx,dy,dz];
urk1 = zeros(nxp2,nyp2,nzp2); vrk1 = urk1; wrk1 = urk1;
ustar = urk1; vstar = urk1; wstar = urk1; prk1 = urk1;
utemp2 = urk1; vtemp2 = urk1; wtemp2 = urk1;
ustar = updatebc(ustar,'u');
vstar = updatebc(vstar,'v');
wstar = updatebc(wstar,'w');
prk1 = updatebc(prk1,'p');

global dt Re
[Crk,dtrk] = coeff_rk3(sub,dt);
a = -dtrk/(Re*dz^2);
b = 1+2*dtrk/(Re*dz^2);
c = -dtrk/(Re*dz^2);

qu1 = dt*urhs(urk0,vrk0,wrk0) + Crk(1)*qu0;
qv1 = dt*vrhs(urk0,vrk0,wrk0) + Crk(1)*qv0;
qw1 = dt*wrhs(urk0,vrk0,wrk0) + Crk(1)*qw0;

utemp1 = urk0 + Crk(2)*qu1;
vtemp1 = vrk0 + Crk(2)*qv1;
wtemp1 = wrk0 + Crk(2)*qw1;

% USTAR
for i = 2:nxp2-2
    for j = 2:nyp2-1
        ustar(i,j,:) = velstar(a,b,c,utemp1(i,j,:),'u');
    end
end
% VSTAR
for i = 2:nxp2-1
    for j = 2:nyp2-2
        vstar(i,j,:) = velstar(a,b,c,vtemp1(i,j,:),'v');
    end
end
% WSTAR
for i = 2:nxp2-1
    for j = 2:nyp2-1
        wstar(i,j,1:nzp2-1) = velstar(a,b,c,wtemp1(i,j,:),'w');
    end
end
ustar = updatebc(ustar,'u');
vstar = updatebc(vstar,'v');
wstar = updatebc(wstar,'w');
% [ustar,vstar,wstar] = fixit(ustar,vstar,wstar);

% for i = 2:nxp2-1
%     for j = 2:nyp2-1
%         for k = 2:nzp2-1
%             utemp2(i,j,k) = (ustar(i+1,j,k)-ustar(i,j,k))/dx;
%             vtemp2(i,j,k) = (vstar(i,j+1,k)-vstar(i,j,k))/dy;
%             wtemp2(i,j,k) = (wstar(i,j,k+1)-wstar(i,j,k))/dz;
%         end
%     end
% end
for i = 2:nxp2-1
    for j = 2:nyp2-1
        for k = 2:nzp2-1
            utemp2(i,j,k) = (ustar(i,j,k)-ustar(i-1,j,k))/dx;
            vtemp2(i,j,k) = (vstar(i,j,k)-vstar(i,j-1,k))/dy;
            wtemp2(i,j,k) = (wstar(i,j,k)-wstar(i,j,k-1))/dz;
        end
    end
end
rhsp = (utemp2 + vtemp2 + wtemp2)/dtrk;
prk1 = GSsolver(d,rhsp,1e-6,6e3);

for i = 2:nxp2-1
    for j = 2:nyp2-1
        for k = 2:nzp2-1
            urk1(i,j,k) = ustar(i,j,k) ...
                - dtrk*(prk1(i+1,j,k)-prk1(i,j,k))/dx; % dxc
            vrk1(i,j,k) = vstar(i,j,k) ...
                - dtrk*(prk1(i,j+1,k)-prk1(i,j,k))/dy; % dyc
            wrk1(i,j,k) = wstar(i,j,k) ...
                - dtrk*(prk1(i,j,k+1)-prk1(i,j,k))/dz; % dzc
        end
    end
end
urk1 = updatebc(urk1,'u');
vrk1 = updatebc(vrk1,'v');
wrk1 = updatebc(wrk1,'w');
% [urk1,vrk1,wrk1] = fixit(urk1,vrk1,wrk1);

end