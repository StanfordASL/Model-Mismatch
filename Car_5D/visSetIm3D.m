function h = visSetIm3D(X,Y,Z,V,level,color, applyLight)
% h = visSetIm3D(g, data, color, level, applyLight)
% Visualizes a 3D reachable set
 
 
% [ mesh_xs, mesh_data ] = gridnd2mesh(g, data);
 
h = patch(isosurface(X,Y,Z,V,level));
isonormals(X,Y,Z,V,h);
h.FaceColor = color;
h.EdgeColor = 'none';
 
if applyLight
  lighting phong
  camlight left
  camlight right
end
 
view(3)
end