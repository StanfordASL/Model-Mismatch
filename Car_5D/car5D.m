clear; close all; clc; 

%% Generate control constraints

acc_lim = 1;
alp_lim = 3;

control_A = [1 ,0;
            -1 ,0;
             0 ,1;
             0, -1];
control_b = [acc_lim;
             acc_lim;
             alp_lim;
             alp_lim]; 

%A*u <= b: g_s(u) = A*u - b <= 0

%% Set up indeterminates

%states: r = (x_r,y_r,th_r,v,om)
r = msspoly('r',5);

%planner controls:
up = msspoly('up',1);

%% Setup dynamics

% Function approximations

% sin_x = @(x)  0.7264*(x/(pi/4));
sin_x = @(x)  0.50586*(x/(pi/6));
% cos_x = @(x) 0.8516 - 0.1464*(2*(x/(pi/4))^2 -1);
cos_x = @(x) 0.9326 - 0.066986*(2*(x/(pi/6))^2 -1);

sin_th = sin_x(r(3));
cos_th = cos_x(r(3));

% Dynamics
v_hat = 1.0;

dynamics.h = [r(4)*cos_th-v_hat+up(1)*r(2);
              r(4)*sin_th-up(1)*r(1);
              r(5)-up(1);
              0;
              0];

dynamics.B = [zeros(3,2);
              eye(2)];
          
%cost function vector          
dynamics.c = [r(1);
              r(2);
              r(3);
              r(4)*cos_th-v_hat;
              r(4)*sin_th;
              r(5)];
          
dynamics.Jc = [eye(3), zeros(3,2);
               0,0,-r(4)*sin_th,cos_th,0;
               0,0, r(4)*cos_th,sin_th,0;
               zeros(1,4),1];

rt = msspoly('rt',length(dynamics.c));

%% Setup constraints

%planner constraints:
up_lim = 0.1;
gp = up(1)^2 - up_lim^2;

%control constraints:
gs = struct('A',control_A,'b',control_b,'N_s',4);

%state-space constraints: 
th_lim = pi/6;
g = rt(3)^2-th_lim^2;

%% Setup tolerances and degrees

all_deg = struct('ly_order',2,'ls_order',2,'lp_order',4,...
                 'lg_order',4,'le_order',4,...
                 'u_order',2,'V_order',2);

toler = struct('delta_E',0.01,...
               'slack_tol',0.00005,'lambda',1.0,...
               'th_1',0.01,'th_2',0.05,...
               'alpha_u',0.1,'alpha_l',0.5,...
               'bt',0.75,...
               'all_deg',all_deg);

%% Initialize V

%Run initialization algorithm or load existing feasible solution:
[V_1,rho_1,E_1,u_sol] = initialize_car5D(dynamics,r,rt,up,gp,gs,g,toler);
% load 'car5D_V_2.mat';

%% Start loop

%ellipsoid weighting
E_w = diag(ones(length(dynamics.c),1));

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
    
    [solved_1,L,u,E_1,gamma_1] = solve_car5D_1(dynamics,r,rt,up,gp,gs,g,V_1,rho_1,E_w,toler);
    
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
        [solved_2,E_2,V_2,rho_2,gamma_2,u_s] = solve_car5D_2(dynamics,r,rt,up,gp,gs,g,u_sol,Ly,Lg,Ls,L_E,E_sol,alpha,E_w,toler);
        
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
        [solved_2,E_2,V_2,rho_2,gamma_2,u_s] = solve_car5D_2(dynamics,r,rt,up,gp,gs,g,u,Ly,Lg,Ls,L_E,E_sol,alpha,E_w,toler);
        
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

save('car5D_soln_raw_2.mat','V_sol','rho_sol','E_sol','u_sol');
E_sol = check_car5D(dynamics,r,rt,up,gp,gs,g,V_sol,rho_sol,u_sol,toler);

%% Extract solution

% pause;

V_sol = mss2fnc(V_sol,rt,randn(length(rt),2));
u_sol = mss2fnc(u_sol,r,randn(length(r),2));

c_sol = @(r)[r(1);
             r(2);
             r(3);
             r(4)*cos(r(3))-v_hat;
             r(4)*sin(r(3));
             r(5)];
         
save('car5D_soln_2.mat','V_sol','rho_sol','E_sol','u_sol','c_sol');

%% Plot ellipsoids

close all; 
figure()
proj_Ellipse([1:2],E_sol,1,[0;0],30,'r')
xlabel('$x_r$','interpreter','latex'); ylabel('$y_r$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',30);set(gca,'FontSize',30)
set(gcf,'Color','w');

figure()
proj_Ellipse([3;6],E_sol,1,[0;0],30,'r')
xlabel('$\theta_r$','interpreter','latex'); ylabel('$\omega$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',30);set(gca,'FontSize',30)
set(gcf,'Color','w');

figure()
proj_Ellipse([4;5],E_sol,1,[0;0],30,'r')
xlabel('$v_{\parallel}$','interpreter','latex'); ylabel('$v_{\perp}$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',30);set(gca,'FontSize',30)
set(gcf,'Color','w');

