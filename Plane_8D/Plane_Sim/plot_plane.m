function [h_p,R] = plot_plane(plane,r,rot,show_box,varargin)

psi = rot(1);
gamma = rot(2);
phi = rot(3);

R_z = [cos(psi), -sin(psi), 0;
       sin(psi), cos(psi), 0;
       0,0,1];
   
gamma = -gamma;

R_y = [cos(gamma), 0, sin(gamma);
       0, 1, 0;
      -sin(gamma), 0, cos(gamma)];
  
R_x = [1,0,0;
       0, cos(phi), -sin(phi);
       0, sin(phi), cos(phi)];
  
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
h_p(1) = trisurf(plane.body.hull,body_V(:,1),body_V(:,2),body_V(:,3),'FaceColor','r','EdgeColor','k','Linewidth',2); hold on
h_p(2) = trisurf(plane.wing.hull,wing_V(:,1),wing_V(:,2),wing_V(:,3),'FaceColor','r','EdgeColor','k','Linewidth',2);
h_p(3) = trisurf(plane.stab.hull,stab_V(:,1),stab_V(:,2),stab_V(:,3),'FaceColor','r','EdgeColor','k','Linewidth',2);
h_p(4) = trisurf(plane.tail.hull,tail_V(:,1),tail_V(:,2),tail_V(:,3),'FaceColor','r','EdgeColor','k','Linewidth',2);

if (show_box)
    box_V = ( R*plane.box.V' + repmat(r,1,plane.box.nV) )';
    trimesh(plane.box.hull,box_V(:,1),box_V(:,2),box_V(:,3),'FaceColor','none','EdgeColor','k');
end
axis equal;

end