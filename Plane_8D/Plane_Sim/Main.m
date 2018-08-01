clear; close all; clc
addpath(genpath('trees'));
%% Load TEB solution

load('../Plane8D_soln_3.mat');

%% Define plane object, parameters, sensing radius

%inflation for TEB (plane CoM)
euc_limits = sqrt(diag(E_sol\eye(size(E_sol,1))));
teb_infl = euc_limits(1:3)';

%inflation for plane body (volume)
vol_infl_xy = [1,1];
vol_infl_z = [0.3;0.4]; %how much to add on top & bottom of obstacle

plane = define_plane_obj(); %for visualization only
teb_box = define_teb_obj(teb_infl);

global grav air_den area mass cd_0 Kd v_hat dT min_turn_r vz_p_lim;

grav = 9.81;
air_den = 1.225;
area = 0.7;
mass = 1;
cd_0 = .015;
Kd = .025;
alpha_0 = 5*pi/180;

C_con = pi*air_den*area/mass;
v_hat = sqrt(grav/(C_con*alpha_0));

phi_max = pi/4;
gamma_max = pi/6;

%tracking control limits
acc_nom = (air_den*area)*(v_hat^2)*(cd_0+4*(pi^2)*Kd*(alpha_0^2));
acc_lim = 10*acc_nom; %linear acceleration lim
phi_d_lim = (pi/3)*2; %phi_dot lim
al_d_lim = 2*gamma_max; %alpha_dot lim

u_lim = [acc_lim;
         phi_d_lim;
         al_d_lim];

sensing_radius = 20;

%dynamics
B = [zeros(4,3);
     1,0,0;
     0,0,0;
     0,1,0;
     0,0,1];

%% Define planning params

World_dim = [120, 120, 28]; 

plan_radius = 50;

%time resolution of paths
dT = 0.01;

om_p_lim = 0.2*(C_con*alpha_0*v_hat*sin(phi_max));
vz_p_lim = 0.1*v_hat;

%planner: x_dot = v*cos, y_dot = v*sin, z_dot = vz, th_dot = om
%normalized: ds= v*dt, where ds is infinitesimal arc length
%dx_ds = cos, dy_ds = sin, dz_ds = vz/v_hat, dth_ds = (om/v_hat) := om_s

om_s_lim = om_p_lim/v_hat;
min_turn_r = 1/om_s_lim; %for planner

%Generate samples and determine nearest neighbors (if not already done)
if exist('DubinsP_FMT_Nodes.mat','file')~=2
    %gen samples
    N_pos_samples = 500;
    N_yaw_samples = 10;
    N_nodes = N_pos_samples*N_yaw_samples;
    
    pos_samples = generateHaltonSamples(3,N_pos_samples);
    pos_samples = pos_samples.*repmat(World_dim,N_pos_samples,1);
    yaw_samples = (linspace(-pi,pi,N_yaw_samples))';
    nodes = [kron(pos_samples,ones(N_yaw_samples,1)),...
             kron(ones(N_pos_samples,1),yaw_samples)];
    [N_nn,T_nn,C_nn] = compute_NN(nodes,plan_radius);
    P_nn = get_steer_paths(nodes,N_nn);
    save('DubinsP_FMT_Nodes.mat','N_nn','P_nn','T_nn','C_nn','nodes','N_nodes');
    vis_nn = 0;
else
    load('DubinsP_FMT_Nodes.mat');
    vis_nn = 0;
end
%visualize steering connections
if (vis_nn)
   figure()
   scatter3(nodes(:,1),nodes(:,2),nodes(:,3)); hold on
   for i = 1:N_nodes
       NN_i = find(N_nn(i,:));
       for j = 1:length(NN_i)
           path = P_nn{i,NN_i(j)};
           plot3(path(:,1),path(:,2),path(:,3),'k-','linewidth',2);
       end
   end
end

%% Define global set of obstacles and goal

new_obstacles = 0;
if (new_obstacles)
    define_obstacles;
    %Check view - if happy, proceed
%     keyboard;
    save('Obstacles.mat','tree_obstacles','tower_obstacles',...
          'obstacles','obstacles_infl','n_obstacles','obstacles_coll');
else
    load('Obstacles.mat');
end

%% Setup simulation

xp_0 = [1,1,5,0]'; %For Dubins
x_0 = [1,1,5,0,v_hat,0,0,alpha_0]'; %For plane

%% Mock simulation
do_mock = 0;

if (do_mock)
    test_Plane8D;
end

%% Find path

%define goal
goal_box = [[World_dim(1:2),15]-[10,10,10];
            World_dim];
goal_V =  get_coll_obstacles(goal_box,1);
goal = Polyhedron('V',goal_V);
        
[FMT_time, X_ref] =  FMTStar_Dubins_Plane(N_nodes,nodes,N_nn,C_nn,T_nn,P_nn,...
                        plan_radius,teb_box,obstacles_coll,xp_0',goal);

X_ref(:,4) = wrapToPi(X_ref(:,4));

%% Plot path
fig_FMT = figure();
title(sprintf('time: %.4f',FMT_time));
plot3(X_ref(:,1),X_ref(:,2),X_ref(:,3),'r-','linewidth',2); hold on
goal.plot('color','blue','alpha',0.3); 
plot_all_obstacles();
xlabel('X'); ylabel('Y'); zlabel('Z');
axis tight
axis equal

%% Simulate

T = (length(X_ref)-1)*dT;
t_ref = 0:dT:T;

simulate_Plane8D;

                    






