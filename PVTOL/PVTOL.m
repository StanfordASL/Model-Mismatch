clear; close all; clc; 

%% Inertial prop

global grav;
mass = 0.486;
Inertia = 0.00383;
len = 0.25;
grav = 9.81;

%% Generate control constraints
thrust_low = 0.1*mass*grav;
thrust_up = 1*mass*grav;
u1 = @(t1,t2) (t1+t2)/mass;
u2 = @(t1,t2) (t1-t2)*(len/Inertia);
thrust_V = [u1(thrust_low,thrust_low),u2(thrust_low,thrust_low);
            u1(thrust_low,thrust_up),u2(thrust_low,thrust_up);
            u1(thrust_up,thrust_low),u2(thrust_up,thrust_low);
            u1(thrust_up,thrust_up),u2(thrust_up,thrust_up)];
thrust_poly = Polyhedron('V',thrust_V).minHRep();

%A*u <= b: g_s(u) = A*u - b <= 0

%% Set up indeterminates

%states: r = (rx,rz,rvx,rvz,phi,om)
r = msspoly('r',6);

%planner controls:
up = msspoly('up',2);

%% Setup dynamics

% Function approximations

% sin_x = @(x)  0.7264*(x/(pi/4)) - 0.01942*(4*(x/(pi/4))^3 - 3*(x/(pi/4)));
% cos_x = @(x) 0.8516 - 0.1464*(2*(x/(pi/4))^2 -1);

sin_x = @(x) 0.9101*(x/(pi/3)) - 0.04466*(4*(x/(pi/3))^3 - 3*(x/(pi/3)));
cos_x = @(x) 0.7441 -0.2499*(2*(x/(pi/3))^2 -1);

sin_p = sin_x(r(5));
cos_p = cos_x(r(5));

% Dynamics

dynamics.h = [r(3);
              r(4);
             -up(1);
             -grav-up(2);
             r(6);
              0];

dynamics.B = [zeros(2);
             -sin_p, 0;
              cos_p, 0;
              0, 0;
              0, 1];
          
dynamics.c = r;
dynamics.Jc = eye(6);

rt = msspoly('rt',length(dynamics.c));

%% Setup constraints

p_lim = pi/3;
up_lim = 0.05*grav;

%planner constraints:
gp = up(1)^2 + up(2)^2 - up_lim^2;

%state-space constraints:  
g = rt(5)^2-p_lim^2;

%control constraints:
gs = struct('A',thrust_poly.A,'b',thrust_poly.b,'N_s',size(thrust_poly.A,1));

%% Setup tolerances

all_deg = struct('ly_order',4,'ls_order',2,'lp_order',4,...
                 'lg_order',2,'le_order',2,...
                 'u_order',2,'V_order',2);
             
toler = struct('delta_E',0.01,'delta_rho',0.01,...
               'slack_tol',0.00005,'lambda',1.0,...
               'th_1',0.01,'th_2',0.05,...
               'alpha_u',0.1,'alpha_l',0.5,...
               'bt',0.75,'all_deg',all_deg);

%% Initialize V

%Run initialization algorithm or load existing feasible solution:
[V_1,rho_1,E_1,u_sol] = initialize_PVTOL(dynamics,r,rt,up,gp,gs,g,toler);
% load 'PVTOL_V.mat';

%% Start loop

%ellipsoid weighting
E_w = diag([1;1;1;1;1;1]);

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
    
    [solved_1,L,u,E_1,gamma_1] = solve_PVTOL_1(dynamics,r,rt,up,gp,gs,g,V_1,rho_1,E_w,toler);
    
    fprintf('Step 1:(%d:%4f), gamma: %.6f \n',solved_1,log(det(E_w*E_1)),gamma_1);
    %% Save solution
    
    if (solved_1 == 1) 
        V_sol = V_1; rho_sol = rho_1; E_sol = E_1; u_sol = u;
        Ly = L{1}; Lg = L{2}; Ls = L{3}; L_E = L{4};
    else
        break;
    end
     
    %% Call problem 2:
    
    fprintf('-------\n');
    
    soln_degrade = 0; 
    alpha = alpha_u;
    while (~soln_degrade)
        [solved_2,E_2,V_2,rho_2,gamma_2,u_s] = solve_PVTOL_2(dynamics,r,rt,up,gp,gs,g,u_sol,Ly,Lg,Ls,L_E,E_sol,alpha,E_w,toler);
        
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
        [solved_2,E_2,V_2,rho_2,gamma_2,u_s] = solve_PVTOL_2(dynamics,r,rt,up,gp,gs,g,u,Ly,Lg,Ls,L_E,E_sol,alpha,E_w,toler);
        
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

E_sol = check_PVTOL(dynamics,r,rt,up,gp,gs,g,V_sol,rho_sol,u_sol,toler);
save('PVTOL_soln_raw.mat','V_sol','rho_sol','E_sol','u_sol');

%% Extract solution

% pause;

V_sol = mss2fnc(V_sol,rt,randn(length(rt),2));
u_sol = mss2fnc(u_sol,r,randn(length(r),2));

save('PVTOL_soln.mat','V_sol','rho_sol','E_sol','u_sol');

%% Plots

close all; 
figure()
proj_Ellipse([1:2],E_sol,1,[0;0],30,'r')
xlabel('$r_x$','interpreter','latex'); ylabel('$r_z$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',30);set(gca,'FontSize',30)
set(gcf,'Color','w');

figure()
proj_Ellipse([3:4],E_sol,1,[0;0],30,'r')
xlabel('$\dot{r}_x$','interpreter','latex'); ylabel('$\dot{r}_z$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',30);set(gca,'FontSize',30)
set(gcf,'Color','w');

figure()
proj_Ellipse([5:6],E_sol,1,[0;0],30,'r')
xlabel('$\theta$','interpreter','latex'); ylabel('$\dot{\theta}$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',30);set(gca,'FontSize',30)
set(gcf,'Color','w');

keyboard;
%% Grid plot

euc_limits = sqrt(diag(E_sol\eye(6)));
n_grid = 200;
x = linspace(-euc_limits(1),euc_limits(1),n_grid);
y = linspace(-euc_limits(2),euc_limits(2),n_grid);
n_grid_other = 10;
vx = linspace(-euc_limits(3),euc_limits(3),n_grid_other);
vy = linspace(-euc_limits(4),euc_limits(4),n_grid_other);
phi = linspace(-euc_limits(5),euc_limits(5),n_grid_other);
phid = linspace(-euc_limits(6),euc_limits(6),n_grid_other);

[X,Y] = meshgrid(x,y);
Z = zeros(n_grid,n_grid,n_grid_other^4);
Z_list = zeros(n_grid_other,n_grid_other,n_grid_other,n_grid_other,n_grid^2);
x_list = reshape(X,1,n_grid^2);
y_list = reshape(Y,1,n_grid^2);
itr = 1;
for i = 1:length(vx)
    for j = 1:length(vy)
        for k = 1:length(phi)
            for ii = 1:length(phid)
                vx_list = vx(i)*ones(1,n_grid*n_grid);
                vy_list = vy(j)*ones(1,n_grid*n_grid);
                phi_list = phi(k)*ones(1,n_grid*n_grid);
                phid_list = phid(ii)*ones(1,n_grid*n_grid);
                
                Z_list(i,j,k,ii,:) = V_sol([x_list;y_list;vx_list;vy_list;phi_list;phid_list]);
                Z(:,:,itr) = reshape(Z_list(i,j,k,ii,:),n_grid,n_grid);
                itr = itr + 1;
            end
        end
    end
end
Z(Z>rho_sol) = nan;

%% Plot
cmap = cbrewer('seq','YlGnBu',64);
figure()
proj_Ellipse([1:2],E_sol,1,[0;0],30,'g')
xlabel('$x_r$','interpreter','latex'); ylabel('$z_r$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',38);set(gca,'FontSize',38)
set(gcf,'Color','w');
hold on
for k = 1:size(Z,3)
    contourf(X,Y,Z(:,:,k),'linestyle','none'); hold on
end
grid on
axis equal
colormap(cmap);
colorbar