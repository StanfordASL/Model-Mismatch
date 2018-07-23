clear all; close all; clc;

ode_options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

%% Load SOS solution
load 'DUBINS_soln.mat';

%% Generate single-integrator trajectory

load 'SI_safe.mat';
%path in 10Hz
T = (size(path,2)-1)*0.1;
t_traj_0 = 0:0.1:T; t_traj_0 = t_traj_0';

%sim at 100Hz
dt = 0.01;
t_traj = 0:dt:T; t_traj = t_traj';
Xp = interp1(t_traj_0,path',t_traj);

Up = diff(Xp)/dt;
Up = [Up;Up(end,:)];

figure()
plot(Xp(:,1),Xp(:,2),'r-');
hold on
grid on

figure()
plot(t_traj,norms(Up,2,2));
grid on
xlabel('time [s]');
ylabel('$\|u_p\|$','interpreter','latex');

%% Simulate

x_0 = [Xp(1,1:2),0,atan2(Up(1,2),Up(1,1))]';

T_steps = length(t_traj);

%Store actual state
x = zeros(T_steps,4);
r = zeros(T_steps-1,3);

%Store control history
u_act = zeros(T_steps-1,2);

V_bnds = zeros(T_steps-1,1);

%Initialize
x(1,:) = x_0';
state = x_0;
      
%% Simulate
disp('Ready to Simulate');
kx = 1; ky = 1; kdx = 1; kdy = 1;

for i = 1:T_steps-1
        
    xp = Xp(i,:)'; 
    
    R = [cos(state(4)),sin(state(4));-sin(state(4)), cos(state(4))];
    r(i,:) = [R*(state(1:2)-xp);state(3)]';
    
    u = u_sol(r(i,:)');
    V_bnds(i) = V_sol(r(i,:)');
       
    u_act(i,:) = u';
   
    
    [d_t,d_x] = ode113(@(t,d_x)dubins_sim(t,d_x,u),[t_traj(i),t_traj(i+1)],state,ode_options);
    
    state = d_x(end,:)';
    state(4) = wrapToPi(state(4));
    x(i+1,:) = state';
end

%% Plot

figure(1)
hold on
plot(x(:,1),x(:,2),'b-'); 

figure()
subplot(2,1,1)
plot(t_traj(1:end-1),r(:,1:2)); 
grid on
xlabel('t'); ylabel('x_r,y_r');

subplot(2,1,2)
plot(t_traj,x(:,4));
grid on
xlabel('t'); ylabel('\theta');

figure()
plot(t_traj(1:end-1),u_act);
grid on
xlabel('t'); ylabel('u');


figure()
plot(t_traj(1:end-1),V_bnds);
hold on
line([0,t_traj(end)],[rho_sol,rho_sol],'color','k','linewidth',2);
grid on
legend('V','\rho');

