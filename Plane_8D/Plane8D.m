clear; close all; clc; 

%% Constants

global grav air_den area mass cd_0 Kd v_hat;

grav = 9.81;
air_den = 1.225;
area = 0.7;
mass = 1;
cd_0 = .015;
Kd = .025;
alpha_0 = 5*pi/180;

C_con = pi*air_den*area/mass;
v_hat = sqrt(grav/(C_con*alpha_0));
acc_nom = (air_den*area)*(v_hat^2)*(cd_0+4*(pi^2)*Kd*(alpha_0^2));

phi_max = pi/4;
gamma_max = pi/6;

%% Generate control constraints

acc_lim = 10*acc_nom; %linear acceleration lim
phi_d_lim = (pi/3)*2; %phi_dot lim
al_d_lim = 2*gamma_max; %alpha_dot lim

control_A = [1 ,0, 0;
            -1 ,0, 0;
             0 ,1, 0;
             0, -1, 0;
             0, 0, 1;
             0, 0, -1]; 
control_b = [acc_lim;
             acc_lim;
             phi_d_lim;
             phi_d_lim;
             al_d_lim;
             al_d_lim]; 

%A*u <= b: g_s(u) = A*u - b <= 0

%% Set up indeterminates

%states: r = (x_r,y_r,z_r,psi_r, v, gamma, phi, alpha)
r = msspoly('r',8);

%planner controls:
up = msspoly('up',2); %up_1: om_hat, up_2: vz_hat

%% Setup dynamics

% Function approximations

ps_max = pi/6;
g_max = pi/6;
ph_max = pi/4;

sin_ps = 0.50586*(r(4)/(pi/6));
cos_ps = 0.9326 - 0.066986*(2*(r(4)/(pi/6))^2 -1);

sin_g = 0.50586*(r(6)/(pi/6));
cos_g = 0.9326 - 0.066986*(2*(r(6)/(pi/6))^2 -1);
sec_g = 1.07502 + 0.07728*(2*(r(6)/(pi/6))^2 -1);

sin_ph = 0.7264*(r(7)/(pi/4));
cos_ph = 0.8516 - 0.1464*(2*(r(7)/(pi/4))^2 -1);

%for v_hat ~ 6.5
v_min = 3; v_max = 10;
bma_v=0.5*(v_max-v_min);
bpa_v=0.5*(v_min+v_max);
v_inv = 0.182574 -0.1067*(r(5)-bpa_v)/bma_v + ...
        0.031181*(2*((r(5)-bpa_v)/bma_v)^2 -1) +...
       -0.00911*(4*((r(5)-bpa_v)/bma_v)^3 - 3*((r(5)-bpa_v)/bma_v));

% Dynamics
%NOTE: NO MORE DRAG or GRAVITY AS IT IS CANCELLED BY THRUST

dynamics.h = [r(5)*cos_g*cos_ps-v_hat+up(1)*r(2); %5
    r(5)*cos_g*sin_ps-up(1)*r(1); %4
    r(5)*sin_g-up(2); %2
    -C_con*r(8)*r(5)*sin_ph*sec_g-up(1); %5
    -grav*sin_g; %4
    C_con*r(5)*r(8)*cos_ph-grav*cos_g*v_inv; %5
    0;
    0];

dynamics.B = [zeros(4,3);
    1,0,0;
    0,0,0;
    0,1,0;
    0,0,1];

dynamics.c = [r(1:4);
              r(5)-v_hat;
              r(6:7);
              r(8)-alpha_0];
dynamics.Jc = eye(8);
           

rt = msspoly('rt',length(dynamics.c));

%% Setup constraints

%planner constraints:
up1_lim = 0.21*(C_con*alpha_0*v_hat*sin(phi_max));
up2_lim = 0.1*v_hat;
gp = [up(1)^2-up1_lim^2;
      up(2)^2-up2_lim^2];

%control constraints:
gs = struct('A',control_A,'b',control_b,'N_s',6);

%state-space constraints:
g = [rt(4)^2-ps_max^2;
     rt(5)^2-bma_v^2;
     rt(6)^2-g_max^2;
     rt(7)^2-ph_max^2];

%% Setup tolerances and degrees

all_deg = struct('ly_order',4,'ls_order',2,'lp_order',4,...
                 'lg_order',2,'le_order',2,...
                 'u_order',2,'V_order',2);

toler = struct('delta_E',0.01,'delta_rho',0.01,...
               'slack_tol',0.08,'lambda',1.0,...
               'th_1',0.01,'th_2',0.05,...
               'alpha_u',0.1,'alpha_l',0.5,...
               'bt',0.75,...
               'all_deg',all_deg);

%% Initialize V

%Run initialization algorithm or load existing feasible solution:
% [V_1,rho_1,E_1,u_sol] = initialize_Plane8D(dynamics,r,rt,up,gp,gs,g,toler);
load 'Plane8D_V_3.mat';

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
L_sol = {};

while ( ~converged )
    %% Call problem 1:
    
    fprintf('*****************\n');
    fprintf('Iteration: %d \n',itr);
    
    [solved_1,L,u,E_1,gamma_1] = solve_Plane8D_1(dynamics,r,rt,up,gp,gs,g,V_1,rho_1,E_w,toler);
    
    fprintf('Step 1:(%d:%4f), gamma: %.6f \n',solved_1,log(det(E_w*E_1)),gamma_1);
    %% Save solution
    
    if (solved_1 == 1) 
        V_sol = V_1; rho_sol = rho_1; E_sol = E_1; u_sol = u;
        L_sol = L;
        Ly = L{1}; Lg = L{2}; Ls = L{3}; L_E = L{4};
    else
        break;
    end
     
    %% Call problem 2:
    
    fprintf('-------\n');
    
    soln_degrade = 0; 
    alpha = alpha_u;
    while (~soln_degrade)
        [solved_2,E_2,V_2,rho_2,gamma_2,u_s] = solve_Plane8D_2(dynamics,r,rt,up,gp,gs,g,u_sol,Ly,Lg,Ls,L_E,E_sol,alpha,E_w,toler);
        
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
        [solved_2,E_2,V_2,rho_2,gamma_2,u_s] = solve_Plane8D_2(dynamics,r,rt,up,gp,gs,g,u,Ly,Lg,Ls,L_E,E_sol,alpha,E_w,toler);
        
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

E_sol = check_Plane8D(dynamics,r,rt,up,gp,gs,g,V_sol,rho_sol,u_sol,toler);

%% Extract solution

% pause;

dV_dr = subs(diff(V_sol,rt),rt,dynamics.c)*dynamics.Jc;

dVsol_dr = mss2fnc(dV_dr',r,randn(length(r),2));
V_sol = mss2fnc(V_sol,rt,randn(length(rt),2));
u_sol = mss2fnc(u_sol,r,randn(length(r),2));

c_sol = @(r) [r(1:4);
              r(5)-v_hat;
              r(6:7);
              r(8)-alpha_0];
         
save('Plane8D_soln_2.mat','V_sol','dVsol_dr','rho_sol','E_sol','u_sol','c_sol');

%% Plot ellipsoids

close all; 
figure()
proj_Ellipse([1:3],E_sol,1,[0;0;0],30,'r')
xlabel('$x_r$','interpreter','latex'); ylabel('$y_r$','interpreter','latex');
zlabel('$z_r$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',30);set(gca,'FontSize',30)
set(gcf,'Color','w');

figure()
proj_Ellipse([4;6;7],E_sol,1,[0;0;0],30,'r')
xlabel('$\psi_r$','interpreter','latex'); ylabel('$\gamma$','interpreter','latex');
zlabel('$\phi$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',30);set(gca,'FontSize',30)
set(gcf,'Color','w');

figure()
proj_Ellipse([5;8],E_sol,1,[0;0],30,'r')
xlabel('$v-\hat{v}$','interpreter','latex'); ylabel('$\alpha$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',30);set(gca,'FontSize',30)
set(gcf,'Color','w');


