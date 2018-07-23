
%% Create Tree obstacles

%location of the trees
numTrees = 30;
tree_pos_xy = generateHaltonSamples(2,numTrees);
tree_pos_xy = tree_pos_xy.*repmat(World_dim(1:2),numTrees,1);

pine_rad = 1;
pine_height = 5;
tree_rad = 1.5;
tree_height = 7;
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

%% Create towers

numTowers = 10;
tower_pos_xy = rand(numTowers,2);
tower_width = 10; tower_height = 15;
tower_pos_xy = tower_pos_xy.*repmat(World_dim(1:2)-tower_width,numTowers,1);

%let's add a few more in the path of the goal
numTowers = numTowers + 2;
tower_pos_xy = [tower_pos_xy;
                105,80;
                105,40];

tower_obstacles = [];

for i = 1:numTowers
    tower = createBoxObs([tower_pos_xy(i,1),tower_pos_xy(i,2),0],[tower_width,tower_width,tower_height]);
    tower_obstacles = [tower_obstacles; tower];
end

%% Inflate

%tightest bounding box form
obstacles = [tree_bounds; tower_obstacles];
n_obstacles = size(obstacles,1)/2;

%add plane box inflation                        
obstacles_infl = obstacles     - kron(repmat([vol_infl_xy,0],n_obstacles,1),[1;0]) +...
                                 kron(repmat([vol_infl_xy,vol_infl_z(1)],n_obstacles,1),[0;1]);                       

%get collision obstacles
obstacles_coll = get_coll_obstacles(obstacles,n_obstacles);

%% Plot 

fig = figure();
plot_all_obstacles();
plot_plane(plane,[1;1;5],[0;0;0],0,fig);

