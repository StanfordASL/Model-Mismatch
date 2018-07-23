%% Load file
clear; close all; clc;

% cmap = cbrewer('seq','YlGnBu',64);
load('car5D_soln_2.mat','V_sol','rho_sol','E_sol');
v_hat = 1;

%% Sylvia 3D meshgrid values (isosurface value: rho_sol)

euc_limits = sqrt(diag(E_sol\eye(size(E_sol,1))));

n_grid = 80; %grid points for (x,y,theta)
n_grid_other = 30; %grid points for other states (unioned)

x = linspace(-euc_limits(1),euc_limits(1),n_grid);
y = linspace(-euc_limits(2),euc_limits(2),n_grid);
theta = linspace(-euc_limits(3),euc_limits(3),n_grid);

%(x_r,y_r,theta_r) meshgrid
[X,Y,Theta] = meshgrid(x,y,theta);

omega = linspace(-euc_limits(6),euc_limits(6),n_grid_other);

%now alternate between (theta_r,vx) and (theta_r,vy)
x_list = reshape(X,1,n_grid^3);
y_list = reshape(Y,1,n_grid^3);
theta_list = reshape(Theta,1,n_grid^3);

V = zeros(n_grid,n_grid,n_grid,2*n_grid_other^2);

%(theta_r,vx)
vx = linspace(-euc_limits(4),euc_limits(4),n_grid_other);
itr = 1;
for i = 1:length(vx)
    vx_list = vx(i)*ones(1,n_grid^3);
    vy_list = (vx_list+v_hat).*tan(theta_list);
    for j = 1:length(omega)
        omega_list = omega(j)*ones(1,n_grid^3);
        Z_list = V_sol([x_list;y_list;theta_list;vx_list;vy_list;omega_list]);
        V(:,:,:,itr) = reshape(Z_list,n_grid,n_grid,n_grid);
        itr = itr + 1;
    end
end

%(theta_r,vy)
vy = linspace(-euc_limits(5),euc_limits(5),n_grid_other);
for i = 1:length(vy)
    vy_list = vy(i)*ones(1,n_grid^3);
    vx_list = (vy_list./tan(theta_list))-v_hat;
    for j = 1:length(omega)
        omega_list = omega(j)*ones(1,n_grid^3);
        Z_list = V_sol([x_list;y_list;theta_list;vx_list;vy_list;omega_list]);
        V(:,:,:,itr) = reshape(Z_list,n_grid,n_grid,n_grid);
        itr = itr + 1;
    end
end

%Final solution: (X,Y,Theta,V); iso_surface val: rho_sol
%kth slice of V: V(:,:,:,k)
%% 3D Plot
figure()
for k = 1:size(V,4)
    hold on
    visSetIm3D(X,Y,Theta,V(:,:,:,k),rho_sol,'r',1);
end

%%
        

% %% Grid 1: (x_r,y_r)
% 
% euc_limits = sqrt(diag(E_sol\eye(size(E_sol,1))));
% n_grid = 200;
% x = linspace(-euc_limits(1),euc_limits(1),n_grid);
% y = linspace(-euc_limits(2),euc_limits(2),n_grid);
% 
% n_grid_other = 10;
% theta = linspace(-euc_limits(3),euc_limits(3),n_grid_other);
% vx = linspace(-euc_limits(4),euc_limits(4),n_grid_other);
% vy = linspace(-euc_limits(5),euc_limits(5),n_grid_other);                                     
% omega = linspace(-euc_limits(6),euc_limits(6),n_grid_other);
% % theta = 0; vx = 0; vy = 0; omega = 0;
% 
% [X,Y] = meshgrid(x,y);
% Z = zeros(n_grid,n_grid,3*n_grid_other^3);
% x_list = reshape(X,1,n_grid^2);
% y_list = reshape(Y,1,n_grid^2);
% 
% %alternate between fixing (th,vx), (th,vy), (vx,vy)
% itr = 1;
% for i = 1:length(theta)
%     theta_list = theta(i)*ones(1,n_grid*n_grid);
%     for j = 1:length(vx)
%         vx_list = vx(j)*ones(1,n_grid*n_grid);
%         vy_list = (vx(j)+v_hat)*tan(theta(i))*ones(1,n_grid*n_grid);
%         for k = 1:length(omega)    
%             omega_list = omega(k)*ones(1,n_grid*n_grid);
%             Z_list = V_sol([x_list;y_list;theta_list;vx_list;vy_list;omega_list]);
%             Z(:,:,itr) = reshape(Z_list,n_grid,n_grid);
%             itr = itr + 1;
%         end
%     end
% end
% 
% for i = 1:length(theta)
%     theta_list = theta(i)*ones(1,n_grid*n_grid);
%     for j = 1:length(vy)
%         vy_list = vy(j)*ones(1,n_grid*n_grid);
%         vx_list = ((vy(j)/tan(theta(i)))-v_hat)*ones(1,n_grid*n_grid);
%         for k = 1:length(omega)    
%             omega_list = omega(k)*ones(1,n_grid*n_grid);
%             Z_list = V_sol([x_list;y_list;theta_list;vx_list;vy_list;omega_list]);
%             Z(:,:,itr) = reshape(Z_list,n_grid,n_grid);
%             itr = itr + 1;
%         end
%     end
% end
% 
% for i = 1:length(vx)
%     vx_list = vx(i)*ones(1,n_grid*n_grid);
%     for j = 1:length(vy)
%         vy_list = vy(j)*ones(1,n_grid*n_grid);
%         theta_list = atan2(vy(j),(vx(i)+v_hat))*ones(1,n_grid*n_grid);
%         for k = 1:length(omega)    
%             omega_list = omega(k)*ones(1,n_grid*n_grid);
%             Z_list = V_sol([x_list;y_list;theta_list;vx_list;vy_list;omega_list]);
%             Z(:,:,itr) = reshape(Z_list,n_grid,n_grid);
%             itr = itr + 1;
%         end
%     end
% end
% 
% Z(Z>rho_sol) = nan;
% 
% figure()
% proj_Ellipse([1:2],E_sol,1,[0;0],30,'g')
% xlabel('$x_r$','interpreter','latex'); ylabel('$y_r$','interpreter','latex');
% set(findall(gcf,'type','text'),'FontSize',38);set(gca,'FontSize',38)
% set(gcf,'Color','w');
% hold on
% for k = 1:size(Z,3)
%     contourf(X,Y,Z(:,:,k),'linestyle','none'); hold on
% end
% colormap(cmap);
% grid on
% axis equal
% colorbar
% 
% clear Z
% 
% %% Grid 2: (v_par, v_perp)
% 
% vx = linspace(-euc_limits(4),euc_limits(4),n_grid);
% vy = linspace(-euc_limits(5),euc_limits(5),n_grid);
% 
% x = linspace(-euc_limits(1),euc_limits(1),n_grid_other);
% y = linspace(-euc_limits(2),euc_limits(2),n_grid_other);
% omega = linspace(-euc_limits(6),euc_limits(6),n_grid_other);
% 
% [VX,VY] = meshgrid(vx,vy);
% Z = zeros(n_grid,n_grid,n_grid_other^3);
% vx_list = reshape(VX,1,n_grid^2);
% vy_list = reshape(VY,1,n_grid^2);
% theta_list = atan2(vy_list,(vx_list+v_hat));
% 
% itr = 1;
% for i = 1:length(x)
%     x_list = x(i)*ones(1,n_grid*n_grid);
%     for j = 1:length(y)
%         y_list = y(j)*ones(1,n_grid*n_grid);
%         for k = 1:length(omega)    
%             omega_list = omega(k)*ones(1,n_grid*n_grid);
%             Z_list = V_sol([x_list;y_list;theta_list;vx_list;vy_list;omega_list]);
%             Z(:,:,itr) = reshape(Z_list,n_grid,n_grid);
%             itr = itr + 1;
%         end
%     end
% end
% Z(Z>rho_sol) = nan;
% 
% figure()
% proj_Ellipse([4:5],E_sol,1,[0;0],30,'g')
% xlabel('$v_{\parallel}$','interpreter','latex'); ylabel('$v_{\perp}$','interpreter','latex');
% set(findall(gcf,'type','text'),'FontSize',38);set(gca,'FontSize',38)
% set(gcf,'Color','w');
% hold on
% for k = 1:size(Z,3)
%     contourf(VX,VY,Z(:,:,k),'linestyle','none'); hold on
% end
% colormap(cmap);
% grid on
% axis equal
% 
% %% Grid 3: (th,omega)
% 
% theta = linspace(-euc_limits(3),euc_limits(3),n_grid);
% omega = linspace(-euc_limits(6),euc_limits(6),n_grid);
% 
% vx = linspace(-euc_limits(4),euc_limits(4),n_grid_other);
% vy = linspace(-euc_limits(5),euc_limits(5),n_grid_other);
% x = linspace(-euc_limits(1),euc_limits(1),n_grid_other);
% y = linspace(-euc_limits(2),euc_limits(2),n_grid_other);
% 
% 
% [TH,OM] = meshgrid(theta,omega);
% Z = zeros(n_grid,n_grid,n_grid_other^3);
% theta_list = reshape(TH,1,n_grid^2);
% omega_list = reshape(OM,1,n_grid^2);
% 
% %alternate between fixing (th,vx), (th,vy)
% itr = 1;
% for i = 1:length(x)
%     x_list = x(i)*ones(1,n_grid*n_grid);
%     for j = 1:length(vx)
%         vx_list = vx(j)*ones(1,n_grid*n_grid);
%         vy_list = (vx_list+v_hat).*tan(theta_list);
%         for k = 1:length(y)    
%             y_list = y(k)*ones(1,n_grid*n_grid);
%             Z_list = V_sol([x_list;y_list;theta_list;vx_list;vy_list;omega_list]);
%             Z(:,:,itr) = reshape(Z_list,n_grid,n_grid);
%             itr = itr + 1;
%         end
%     end
% end
% 
% for i = 1:length(x)
%     x_list = x(i)*ones(1,n_grid*n_grid);
%     for j = 1:length(vy)
%         vy_list = vy(j)*ones(1,n_grid*n_grid);
%         vx_list = ((vy_list./tan(theta_list))-v_hat);
%         for k = 1:length(y)    
%             y_list = y(k)*ones(1,n_grid*n_grid);
%             Z_list = V_sol([x_list;y_list;theta_list;vx_list;vy_list;omega_list]);
%             Z(:,:,itr) = reshape(Z_list,n_grid,n_grid);
%             itr = itr + 1;
%         end
%     end
% end
% 
% Z(Z>rho_sol) = nan;
% 
% figure()
% proj_Ellipse([3;6],E_sol,1,[0;0],30,'g')
% xlabel('$\theta_r$','interpreter','latex'); ylabel('$\omega$','interpreter','latex');
% set(findall(gcf,'type','text'),'FontSize',38);set(gca,'FontSize',38)
% set(gcf,'Color','w');
% hold on
% for k = 1:size(Z,3)
%     contourf(TH,OM,Z(:,:,k),'linestyle','none'); hold on
% end
% colormap(cmap);
% colorbar;
% grid on
% axis equal
