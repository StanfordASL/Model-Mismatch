clear all; close all; clc;

%% Load SOS solution
load '../car5D_soln_2.mat'; %V_sol, u_sol, rho_sol, E_sol, c_sol

%% Load 3D Dubins trajectory

% load 'Dubins_traj.mat';

xp_0 = [0;0;0];
xp_f = [7;7;pi/2];
v_hat = 1;
om_lim = 0.1;
min_turn_r = v_hat/om_lim;
dt = 0.01;

[X_ref, ~, ~] = dubins_curve(xp_0,xp_f,min_turn_r,v_hat*dt,1);
X_ref(:,3) = wrapToPi(X_ref(:,3));
T = (length(X_ref)-1)*dt;
t_traj = 0:dt:T;

%% Simulate

% (x,y,th,v,om)
N = 10; %number of sims
start_rand_p = [-0.02;-0.1] + [0.04;0.2].*rand(2,N);
start_rand_th = zeros(1,N);
start_rand_om = zeros(1,N);
x_0 =   [kron(ones(1,N),X_ref(1,1:2)')+start_rand_p;
         kron(ones(1,N),X_ref(1,3))+start_rand_th;
         v_hat*ones(1,N);
         start_rand_om]; 

T_steps = length(t_traj);

%Store actual state
x = zeros(T_steps,5,N);
r = zeros(T_steps-1,5,N);

%Store control history
u_act = zeros(T_steps-1,2,N);

V_bnds = zeros(T_steps-1,N);

%Initialize
for k = 1:N
    x(1,:,k) = x_0(:,k)';
end
      
%% Simulate
disp('Ready to Simulate');
ode_options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

for k = 1:N
    state = x_0(:,k);
    for i = 1:T_steps-1
        
        xp = X_ref(i,:)'; 
        
        R_th = [cos(xp(3)), sin(xp(3));
               -sin(xp(3)), cos(xp(3))];
        r(i,:,k) = [R_th*(state(1:2)-xp(1:2));
                    wrapToPi(state(3)-xp(3));
                    state(4);
                    state(5)]';
        
        u_safety = u_sol(r(i,:,k)');
        V_bnds(i,k) = V_sol(c_sol(r(i,:,k)'));
        
        u = u_safety;
        
        u_act(i,:,k) = u';
        
        [d_t,d_x] = ode113(@(t,d_x)car5D_sim(t,d_x,u),[t_traj(i),t_traj(i+1)],state,ode_options);
        
        d_x(:,3) = wrapToPi(d_x(:,3));
        state = d_x(end,:)';
        x(i+1,:,k) = state';
    end
end

%% Plot

close all; 

figure(1)
hold on
set(gcf,'Color','w');
for i = 1:(0.5/dt):length(t_traj)
    proj_Ellipse([1:2],E_sol,1,[X_ref(i,1);X_ref(i,2)],30,'b',wrapToPi(X_ref(i,3)))
end
plot(X_ref(:,1),X_ref(:,2),'g-','linewidth',2); 
grid off
hold on
for k = 1:N
    plot(x(:,1,k),x(:,2,k),'r-','linewidth',1.0); 
end
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

figure()
subplot(2,1,1)
hold on
for k = 1:N
    plot(t_traj(1:end-1),u_act(:,1,k),'r-','linewidth',2);
end
line([t_traj(1),t_traj(end)],[-1,-1],'color','k','linewidth',2);
line([t_traj(1),t_traj(end)],[1,1],'color','k','linewidth',2);
grid on; xlabel('time [s]'); ylabel('a');

subplot(2,1,2)
hold on
for k = 1:N
    plot(t_traj(1:end-1),u_act(:,2,k),'r-','linewidth',2);
end
line([t_traj(1),t_traj(end)],[-3,-3],'color','k','linewidth',2);
line([t_traj(1),t_traj(end)],[3,3],'color','k','linewidth',2);
grid on; xlabel('time [s]'); ylabel('\alpha');


figure()
for k = 1:N
    plot(t_traj(1:end-1),V_bnds(:,k),'linewidth',2);
    hold all;
end
line([0,t_traj(end)],[rho_sol,rho_sol],'color','k','linewidth',2);
grid on
xlabel('Time [s]');
ylabel('V(r)');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
set(gcf,'Color','w');


