function [h_p,R] = plot_plane(plane,r,rot,show_box,varargin)

yaw = rot(1);
pitch = rot(2);
roll = rot(3);

R_z = [cos(yaw), -sin(yaw), 0;
       sin(yaw), cos(yaw), 0;
       0,0,1];
   
pitch = -pitch; %since not using NED

R_y = [cos(pitch), 0, sin(pitch);
       0, 1, 0;
      -sin(pitch), 0, cos(pitch)];
  
R_x = [1,0,0;
       0, cos(roll), -sin(roll);
       0, sin(roll), cos(roll)];
  
R = R_z*R_y*R_x;

%% distort and translate

body_V = ( R*plane.body.V' + repmat(r,1,plane.body.nV) )';
wing_V = ( R*plane.wing.V' + repmat(r,1,plane.wing.nV) )';
stab_V = ( R*plane.stab.V' + repmat(r,1,plane.stab.nV) )';
tail_V = ( R*plane.tail.V' + repmat(r,1,plane.tail.nV) )';

%% plot

if nargin > 4
    figure(varargin{1})
else
    figure()
end
h_p(1) = trisurf(plane.body.hull,body_V(:,1),body_V(:,2),body_V(:,3),'FaceColor','r','EdgeColor','k','Linewidth',1); hold on
h_p(2) = trisurf(plane.wing.hull,wing_V(:,1),wing_V(:,2),wing_V(:,3),'FaceColor','r','EdgeColor','k','Linewidth',1);
h_p(3) = trisurf(plane.stab.hull,stab_V(:,1),stab_V(:,2),stab_V(:,3),'FaceColor','r','EdgeColor','k','Linewidth',1);
h_p(4) = trisurf(plane.tail.hull,tail_V(:,1),tail_V(:,2),tail_V(:,3),'FaceColor','r','EdgeColor','k','Linewidth',1);

if (show_box)
    box_V = ( R*plane.box.V' + repmat(r,1,plane.box.nV) )';
    trimesh(plane.box.hull,box_V(:,1),box_V(:,2),box_V(:,3),'FaceColor','none','EdgeColor','k');
end
axis equal;

end