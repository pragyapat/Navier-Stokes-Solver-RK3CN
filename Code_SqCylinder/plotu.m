uplot1(:,:) = urk1(:,17,:);
uplot2(:,:) = urk2(:,17,:);
uplot3(:,:) = urk3(:,17,:);
figure (1)
hold on
title(['u_{mid} at 1^{st} RK-substep (t=' num2str(t) 's)'])
contourf(uplot1')
axis equal
colorbar
xlabel('x')
ylabel('z')
hold off

figure (2)
hold on
title(['u_{mid} at 2^{nd} RK-substep (t=' num2str(t) 's)'])
contourf(uplot2')
axis equal
colorbar
xlabel('x')
ylabel('z')
hold off

figure (3)
hold on
title(['u_{mid} at 3^{rd} RK-substep (t=' num2str(t) 's)'])
contourf(uplot3')
axis equal
colorbar
xlabel('x')
ylabel('z')
hold off