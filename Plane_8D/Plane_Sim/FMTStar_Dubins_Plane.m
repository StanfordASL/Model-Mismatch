function [FMT_time, path_final] = FMTStar_Dubins_Plane(N,V,N_nn,C_nn,T_nn,P_nn,EPS,teb,obs,x0,goal)

% EPS = 0.5;
% n: number of total samples
% V: N x 4 matrix of (x,y,z,theta) tuples
% N_mm: N x N sparse matrix; N_idx(i,:) set of nearest neighbors for node i
% C_nn: N x N sparse matrix; Cost(i,:) connection costs to neighbors for node i
% x0: start node (x,y,z,theta) tuple

%% Augment Nodes
V = [x0;V];
N = N+1;

N_nn = blkdiag(0,N_nn);
C_nn = blkdiag(0,C_nn);
T_nn = blkdiag(0,T_nn);

%% Initialize

%unvisited nodes (all except start)
V_unvisited = ones(N,1); 
V_unvisited(1) = 0;

%open nodes (none except start)
V_open = zeros(N,1); V_open(1) = 1;
C_open = Inf(N,1); C_open(1) = 0;

%graph nodes parent node indexes
V_graph = zeros(N,1);

%expand node (idx)
z = 1;
%get N_nn just for start node
[N_nn(1,:),C_nn(1,:),T_nn(1,:),P_init] = get_NN(1,V,N,EPS);
P_pad = cell(N-1,1);
P_nn = [P_init;
        P_pad,P_nn];

z_in_goal = 0;

%% Search

tic
while ~z_in_goal && sum(V_open) > 0
    
    %fwd reachable nodes of z
    NF_z = find(N_nn(z,:)); %idx of nn(z)
    
    for i = 1:length(NF_z)
        x = NF_z(i);
        
        %get unvisited neighbors of z
        if V_unvisited(x) 
            
            %bwd reachable nodes of x
            BN_x = find(N_nn(:,x)); %idx of nn(x)
            
            %find lowest cost node in V_open to x
            c_min = Inf; y_min = Inf;
            for j = 1:length(BN_x)
                y = BN_x(j);
                
                if V_open(y) 
                    c = C_open(y) + C_nn(y,x);
                    if c < c_min
                        c_min = c;
                        y_min = y; %store min node index
                    end
                end
            end
            if (y_min < Inf) && ~InCollision(P_nn{y_min,x},teb,obs)
                %add x to graph
                V_graph(x) = y_min; 
                %add x to V_open
                V_open(x) = 1; C_open(x) = c_min;
                %remove x from V_unvisited
                V_unvisited(x) = 0;
            end
        end
    end

    %remove z from V_open
    V_open(z) = 0;
    C_open(z) = Inf;
    
    %find next min cost node in V_open
    [~,z] = min(C_open);
    
    %check if in goal
    z_in_goal = goal.contains(V(z,1:3)');
    
end
FMT_time = toc;

if (~z_in_goal)
    %exited with failure
    disp('Failure');
    path_final = 0;
    return;
end
%% Plot tree
% fig_fmt = figure();
% hold on
% for i = 2:N
%     p = V_graph(i);
%     if (p~=0)
%         path_loc = P_nn{p,i};
%         plot3(path_loc(:,1),path_loc(:,2),path_loc(:,3),'k-');
%     end
% end
%% Recover optimal path

% Search backwards from goal to start to find the optimal least cost path
path_idx = []; path_idx(1) = z; %start at end
q = z;
while (V_graph(q)~=0)
   %get parent
   p = V_graph(q);
   path_idx = [path_idx; p];
   q = p;   
end
%flip
path_idx = flipud(path_idx);
path_final = P_nn{path_idx(1),path_idx(2)};
for i = 2:length(path_idx)-1
    path_final = [path_final;
                  P_nn{path_idx(i),path_idx(i+1)}(2:end,:)];
end

end
