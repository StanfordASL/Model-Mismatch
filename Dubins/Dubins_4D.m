clear; close all; clc; 

%% Generate control constraints

acc_lim = 1;
om_lim = 1;

control_A = [1 ,0;
            -1 ,0;
             0 ,1;
             0, -1];
control_b = [acc_lim;
             acc_lim;
             om_lim;
             om_lim];  
         
%% Set up indeterminates

%states: r = (rx,ry,v)
r = msspoly('r',3);

%planner controls:
up = msspoly('up',2);

%% Setup dynamics

dynamics.h = [r(3)-up(1);
               -up(2);
               0];

dynamics.B = [0, r(2);
              0, -r(1);
              1, 0];

dynamics.c = r;
dynamics.Jc = eye(3);
rt = msspoly('rt',3);

%% Setup constraints

up_lim = 0.1;

%planner constraints:
gp = up(1)^2 + up(2)^2 - up_lim^2;
     
%control constraints: A*u - b <= 0
gs = struct('A',control_A,'b',control_b,'N_s',4);	

%% Setup tolerances

all_deg = struct('ly_order',4,'ls_order',4,'lp_order',6,...
                 'le_order',4,...
                 'u_order',2,'V_order',4);

toler = struct('delta_E',0.01,...
               'slack_tol',0.00005,'lambda',5.0,...
               'th_1',0.01,'th_2',0.05,...
               'alpha_u',0.1,'alpha_l',0.5,...
               'bt',0.75,...
               'all_deg',all_deg);
          
%% Initialize V

%Run initialization algorithm or load existing feasible solution:
[V_1,rho_1,E_1,u_sol] = initialize_Dubins(dynamics,r,rt,up,gp,gs,toler);
% load 'Dubins_V.mat';

%% Start loop

%ellipsoid weighting
E_w = diag([1;1;1]);

%convergence params
obj_prev = log(det(E_w*E_1));
converged = 0;
itr = 1;

%Upper-bound on backtrack alpha
alpha_u = toler.alpha_u;

%initialize solution
V_sol = V_1; rho_sol = rho_1; E_sol = E_1; 

while ( ~converged )
    %% Call problem 1:
    
    fprintf('*****************\n');
    fprintf('Iteration: %d \n',itr);
    
    [solved_1,L,u,E_1,gamma_1] = solve_DUBINS_1(dynamics,r,rt,up,gp,gs,V_1,rho_1,E_w,toler);
    
    fprintf('Step 1:(%d:%4f), gamma: %.6f \n ',solved_1,log(det(E_w*E_1)),gamma_1);
    %% Save solution
    
    if (solved_1 == 1) 
        V_sol = V_1; rho_sol = rho_1; E_sol = E_1; u_sol = u;
        Ly = L{1}; Ls = L{2}; L_E = L{3};
    else
        break;
    end
     
    %% Call problem 2:
    
    fprintf('-------\n');
    
    soln_degrade = 0; 
    alpha = alpha_u;
    while (~soln_degrade)
        [solved_2,E_2,V_2,rho_2,gamma_2,u_s] = solve_DUBINS_2(dynamics,r,rt,up,gp,gs,u,Ly,Ls,L_E,E_sol,alpha,E_w,gamma_1,toler);
        
        fprintf('Step 2:(%d:%4f), (gamma,u_s):(%.6f,%.6f) \n',solved_2,log(det(E_w*E_2)),gamma_2,u_s);
        if max(gamma_2,u_s) <= toler.slack_tol
            %save solution
            V_sol = V_2; rho_sol = rho_2; E_sol = E_2; 
        else
            fprintf('start BT\n');
            soln_degrade = 1;
        end
    end
    %backtrack
    while (soln_degrade) && (alpha>toler.th_2*alpha_u)
        alpha = toler.bt*alpha;
        [solved_2,E_2,V_2,rho_2,gamma_2,u_s] = solve_DUBINS_2(dynamics,r,rt,up,gp,gs,u,Ly,Ls,L_E,E_sol,alpha,E_w,gamma_1,toler);
        
        fprintf('Step 2:(%d:%4f), (gamma,u_s):(%.6f,%.6f) \n',solved_2,log(det(E_w*E_2)),gamma_2,u_s);
        if max(gamma_2,u_s) <= toler.slack_tol
            %backtracked sufficiently; save solution
            V_sol = V_2; rho_sol = rho_2; E_sol = E_2; 
            soln_degrade = 0;
        end        
    end  
    
    %% Reset for prob 1
    V_1 = V_sol; rho_1 = rho_sol;
    
    %% Evaluate termination based on cost improvement
    
    obj = log(det(E_w*E_sol));
    dec_frac = ((obj-obj_prev)/abs(obj_prev));
    fprintf('eps: %f \n', dec_frac);
    if (dec_frac <= toler.th_1)
        break;
    else
        obj_prev = obj;
        itr = itr + 1;
    end
    
end

%% Check solution

E_sol = check_DUBINS(dynamics,r,rt,up,V_sol,rho_sol,u_sol,gp,gs,toler);
save('DUBINS_soln_raw.mat','V_sol','rho_sol','E_sol','u_sol');

%% Extract solution

% pause;

V_sol = mss2fnc(V_sol,rt,randn(length(rt),2));
u_sol = mss2fnc(u_sol,r,randn(length(r),2));

save('Dubins_soln.mat','V_sol','rho_sol','E_sol','u_sol');

%% Plot ellipsoids

close all; 
figure()
subplot(2,1,1)
proj_Ellipse([1:2],E_sol,1,[0;0],30,'r')
xlabel('x'); ylabel('y');
subplot(2,1,2)
proj_Ellipse([2:3],E_sol,1,[0;0],30,'r')
xlabel('y'); ylabel('v')

keyboard;
%% Grid computation
n_grid = 300;
euc_limits = sqrt(diag(E_sol\eye(3)));
x = linspace(-euc_limits(1),euc_limits(1),n_grid);
y = linspace(-euc_limits(2),euc_limits(2),n_grid);
n_V = 20;
v = linspace(-euc_limits(3),euc_limits(3),n_V);

[X,Y] = meshgrid(x,y);
Z = zeros(n_grid,n_grid,n_V);
Z_list = zeros(n_V,n_grid^2);
x_list = reshape(X,1,n_grid^2);
y_list = reshape(Y,1,n_grid^2);
for k = 1:length(v)
    v_list = v(k)*ones(1,n_grid*n_grid);
    Z_list(k,:) = V_sol([x_list;y_list;v_list]);
    Z(:,:,k) = reshape(Z_list(k,:),n_grid,n_grid);
end
Z(Z>rho_sol) = nan;

theta = linspace(-pi,pi,100);

%% Full plot

cmap = cbrewer('seq','YlGnBu',64);
%Single (v,theta) slice
figure()
th = 0;
k = 11; %v = 0.013
contourf(cos(th)*X-sin(th)*Y,sin(th)*X+cos(th)*Y,Z(:,:,k),'linestyle','none'); hold on
colormap(cmap);
xlabel('$e_x''$','interpreter','latex'); ylabel('$e_y''$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',40);set(gca,'FontSize',40)
set(gcf,'Color','w');
grid on
axis e
colorbar

%all (v,theta)    
figure()
for j = 1:length(theta)
    proj_Ellipse([1:2],E_sol,1,[0;0],30,'r',theta(j)); hold on;
end
for k = 1:n_V
    for j = 1:length(theta)
        th = theta(j);
        contourf(cos(th)*X-sin(th)*Y,sin(th)*X+cos(th)*Y,Z(:,:,k),'linestyle','none'); hold on
    end   
end
colormap(cmap);
grid on
axis equal
colorbar
xlabel('$e_x''$','interpreter','latex'); ylabel('$e_y''$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',40);set(gca,'FontSize',40)
set(gcf,'Color','w');

