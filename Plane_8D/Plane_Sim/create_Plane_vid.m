close all; 

fig = figure();
plot3(X_ref(:,1),X_ref(:,2),X_ref(:,3),'b-','linewidth',2); hold on
plot3(X(:,1),X(:,2),X(:,3),'k-','linewidth',2.0);
goal.plot('color','blue','alpha',0.3); 

%% FPV view:
fpv_view = 0;
set(gca,'CameraViewAngleMode','manual')
camproj('perspective')
camva(45)

%%
t_step = round((1/10)/dT);
t_step_b = 0.5/dT;
T_b = 2; %2 second lookahead for bound
n_b = 1+(T_b/0.5);

i = 1;
pos = X(i,1:3)';
att = X(i,[4;6;7]);

%% create objects

hold on
if (~fpv_view)
    h_b = [];
    for j = 1:n_b
        j_idx = min(i + (j-1)*t_step_b, length(t_ref));
        h_b(j) = proj_Ellipse([1:3],E_sol,1,X_ref(j_idx,1:3)',30,'b',X_ref(j_idx,4));
    end
end
[h_p,rot] = plot_plane(plane,pos,att,0,fig);

plot3dObstacles(tree_obstacles,'g',0.6); hold on
plot3dObstacles(tower_obstacles,[0.5,0.5,0.5],0.8);
plot3dObstacles(obstacles_infl,'r',0.05,1);
plot3dObstacles([0 0 0; World_dim],'k',0);
patch([0,World_dim(1),World_dim(1),0],...
      [0, 0, World_dim(2), World_dim(2)],...
      [0, 0,  0, 0],'FaceColor',[0,0.5,0],'FaceAlpha',0.5);
  
hold off

view(3);

%set view
if (fpv_view)
    pos_tail = [-2;0;1.0];
    pos_view = [3;0;-0.5];
    campos(pos'+(rot*pos_tail)')
    camtarget(pos'+(rot*(pos_tail+pos_view))')
else
    campos(pos'+[-7.5,-9.9,7.18]);
    camtarget(pos');
end
drawnow;

%% record movie

record_vid = 1;
if (record_vid)
    writerObj = VideoWriter('Plane_sim_2.mp4');
    writerObj.FrameRate = 1/(t_step*dT);
    writerObj.Quality = 100;
    open(writerObj);
    set(gcf, 'renderer', 'zbuffer')
end

%% Go
for i = 1:t_step:length(t_ref)
% for i = (20/dT)+1:10*t_step:length(t_ref)
    pos = X(i,1:3)';
    att = X(i,[4;6;7]);
    
    delete(h_p); 
    hold on
    %update plane and bound
    if (~fpv_view)
        delete(h_b);
        h_b = [];
        for j = 1:n_b
            j_idx = min(i + (j-1)*t_step_b, length(t_ref));
            h_b(j) = proj_Ellipse([1:3],E_sol,1,X_ref(j_idx,1:3)',30,'b',X_ref(j_idx,4));
        end
    end
    [h_p,rot] = plot_plane(plane,pos,att,0,fig);
    
    hold off

    %update view
    if (fpv_view)
        campos(pos'+(rot*pos_tail)')
        camtarget(pos'+(rot*(pos_tail+pos_view))')
    else
        campos(pos'+[-7.5,-9.9,7.18]);
        camtarget(pos');
    end
    drawnow;
    
    %record
    if (record_vid)
        thisFrame = getframe(gcf);
        %Write this frame out to a new video file.
        writeVideo(writerObj, thisFrame);
    end
end 

if (record_vid)
    close(writerObj);
end
    