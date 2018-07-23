function [h] = visualize5D3D(g, data, teb_level, view3D, viewXY, viewVV, viewWT)
% dimensions = [xr, yr, tr, v, w]
data = -data;


if nargin <3
    teb_level = min(data(:)) + .2;
end

if nargin <4
    view3D = 1;
end

if nargin <5
    viewXY = 1;
end

if nargin <6
    viewVV = 1;
end

if nargin <7
    viewWT = 1;
end


%% Pretty 3D figure
% This isn't required by Sumeet but I think it helps with the
% interpretation/intuition more.

if view3D
    
    % dimensions = [xr, yr, tr, v, w]
    
    [g3D{1}, data3D{1}] = proj(g, data, [0 0 0 1 1], 'min');
    [g3D{2}, data3D{2}] = proj(g, data, [0 0 1 0 1], 'min');
    [g3D{3}, data3D{3}] = proj(g, data, [0 0 1 1 0], 'min');
    
    figure(1)
    clf
    hold on
    
    zaxis = {'$\theta$';'$v$';'$\omega$'};
    
    for ii = 1:3
        g3D_temp = g3D{ii};
        data3D_temp = data3D{ii};
        
        subplot(1, 3, ii)
        h{1} = visSetIm(g3D_temp, data3D_temp, 'red', teb_level);
                
        axis([g3D_temp.min(1) g3D_temp.max(1) g3D_temp.min(2) ...
            g3D_temp.max(2) g3D_temp.min(3) g3D_temp.max(3)])
        xlabel('$e_x''$','interpreter','latex')
        ylabel('$e_y''$','interpreter','latex');
        
        
        zlabel(zaxis{ii},'interpreter','latex');
        set(gcf,'Color','w');
        title('3D Projection (Union)','interpreter','latex')
        axis equal
        
        lighting phong
        camlight left
        camlight right
    end
end

%% Projection over (x,y): union over the rest

if viewXY
    
    cmap = cbrewer('seq','YlGnBu',64);
    colormap(cmap);
    
    figure(2)
    clf
    
    % project to relevant dimensions (x,y)
    % dimensions = [xr, yr, tr, v, w]
    [g2D, data2D_temp] = proj(g, data, [0 0 1 1 1],'min');
    
    % remove values above teb level
    data2D_temp(data2D_temp>teb_level) = NaN;
    
    % plot
    h{2} = contourf(g2D.xs{1},g2D.xs{2},...
        data2D_temp,'linestyle','none'); hold on
    
    xlabel('$e_x''$','interpreter','latex'); ylabel('$e_y''$','interpreter','latex');
    set(findall(gcf,'type','text'),'FontSize',38);set(gca,'FontSize',38)
    set(gcf,'Color','w');
    grid on
    axis equal
    axis([g2D.min(1) g2D.max(1) g2D.min(2) g2D.max(2)])
    colorbar;
    
end
%% Projection over (v_parallel, v_perp): sweep over the rest
% v_parallel = vcos(theta)
% v_perp = vsin(theta)

% NOTE: not "sweeping" over this because it sucks to do.  Instead just
% projecting down to the 2 relevant dimensions (theta and v) and then doing
% the transform

if viewVV
    
    theta = linspace(g.min(3), g.max(3), 100);
    vel = linspace(g.min(4), g.max(4), 100);
    
    % Define new reference frame
    v_par = vel.*cos(theta);
    v_perp = vel.*sin(theta);
    
    % Get the grid for the new reference frame
    [V_Par, V_Perp] = meshgrid(v_par, v_perp);
    
    
    cmap = cbrewer('seq','YlGnBu',64);
    colormap(cmap);
    
    figure(3)
    clf
    
    % Project down to theta, v
    % dimensions = [xr, yr, tr, v, w]
    [g2D, data2D_temp] = proj(g, data, [1 1 0 0 1], 'min');
    
    % Set any values above teb_level to NaN so they don't visualize
    data2D_temp(data2D_temp>teb_level) = NaN;
    
    % tranform the data to the new reference fram
    data_transform = zeros(length(v_par), length(v_perp));
    for ll = 1:length(v_par)
        for mm = 1:length(v_perp)
            
            % want to transform to v_par, v_perp,
            % where v_par = vcos(theta), v_perp = vsin(theta)
            % v_perp/v_par = (vsin(theta))/(vcos(theta)) = tan (theta)
            % theta = arctan (v_perp/v_par)
            % v_par = vcos(arctan(v_perp/v_par)) --> math -->
            % v = v_par*sqrt((v_perp^2/v_par^2)+1)
            
            theta_eval = atan(v_perp(mm)/v_par(ll));
            vel_eval = ...
                v_par(ll) *sqrt((v_perp(mm)^2/v_par(ll)^2)+1);
            
            
            data_transform(ll, mm) = eval_u(g2D, data2D_temp,...
                [theta_eval, vel_eval]);
        end
    end
    
    % plot
    h{3} = contourf(V_Par,V_Perp,data_transform,'linestyle','none'); 
    hold on
    
    xlabel('$v_{\parallel}$','interpreter','latex');
    ylabel('$v_{\bot}''$','interpreter','latex');
    set(findall(gcf,'type','text'),'FontSize',38);set(gca,'FontSize',38)
    set(gcf,'Color','w');
    grid on
    axis equal
    axis([min(v_par) max(v_par) min(v_perp) max(v_perp)])
    colorbar;
end

%% Projection over (omega, theta_r): sweep over the rest

if viewWT
    cmap = cbrewer('seq','YlGnBu',64);
    colormap(cmap);
    
    figure(4)
    clf
    
    % project to relevant dimensions (theta, omega)
    % dimensions = [xr, yr, tr, v, w]
    [g2D, data2D_temp] = proj(g, data, [1 1 0 1 0],'min');
    
    % remove values above teb level
    data2D_temp(data2D_temp>teb_level) = NaN;
    
    % plot
    h{4} = contourf(g2D.xs{1},g2D.xs{2},...
        data2D_temp,'linestyle','none'); hold on
    
    xlabel('$\theta_r$','interpreter','latex'); 
    ylabel('$\omega$','interpreter','latex');
    set(findall(gcf,'type','text'),'FontSize',38);set(gca,'FontSize',38)
    set(gcf,'Color','w');
    grid on
    axis equal
    axis([g2D.min(1) g2D.max(1) g2D.min(2) g2D.max(2)])
    colorbar;
end

