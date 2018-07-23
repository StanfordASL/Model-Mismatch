
%% Generate some towers and trees

numTrees = 4;
tree_pos_xy = [65,5;
               50,50;
               50,80;
               80,50];
           
numTowers = 2;
tower_width = 8; tower_height = 12;
tower_pos_xy = [60,40;
                90,75];
tower_obstacles = [];

for i = 1:numTowers
    tower = createBoxObs([tower_pos_xy(i,1),tower_pos_xy(i,2),0],[tower_width,tower_width,tower_height]);
    tower_obstacles = [tower_obstacles; tower];
end           
                      
pine_rad = 2;
pine_height = 6;
tree_rad = 3;
tree_height = 6;
numTop = 6;
numDisc = 4;

tree_bounds = [];
tree_obstacles = [];

for i = 1:numTrees
    if rand(1) > 0.5
        tree_obstacles = [tree_obstacles;
                         createPineTreeObs([tree_pos_xy(i,:) 0], pine_rad, pine_height, numTop, numDisc)];
        tree_bounds = [tree_bounds;
                     tree_pos_xy(i,1)-pine_rad,tree_pos_xy(i,2)-pine_rad,0;
                     tree_pos_xy(i,1)+pine_rad,tree_pos_xy(i,2)+pine_rad,pine_height];
    else
        tree_obstacles = [tree_obstacles;
                         createTreeObs([tree_pos_xy(i,:) 0],tree_rad, tree_height, numTop, numDisc)];
        tree_bounds = [tree_bounds;
                     tree_pos_xy(i,1)-tree_rad,tree_pos_xy(i,2)-tree_rad,0;
                     tree_pos_xy(i,1)+tree_rad,tree_pos_xy(i,2)+tree_rad,tree_height];
    end
end  

figure(fig_FMT)
plot3dObstacles(tree_obstacles,'g',0.6); hold on
plot3dObstacles(tree_bounds,'r',0.1);
plot3dObstacles(tower_obstacles,'k',0.4);


               