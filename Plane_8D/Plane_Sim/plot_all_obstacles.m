plot3dObstacles(tree_obstacles,'g',0.6); hold on
plot3dObstacles(tower_obstacles,'r',0.6);
plot3dObstacles(obstacles_infl,'r',0.05);
plot3dObstacles([0 0 0; World_dim],'k',0);
patch([0,World_dim(1),World_dim(1),0],...
      [0, 0, World_dim(2), World_dim(2)],...
      [0, 0,  0, 0],'FaceColor',[0,0.5,0],'FaceAlpha',0.5);
  
axis equal
axis tight