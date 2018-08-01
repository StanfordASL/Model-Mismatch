clear all; close all; clc;

%% Load SOS solution
load '../PVTOL_soln.mat'; %V_sol, u_sol, rho_sol, E_sol

%% Inertial prop

global m J l grav;

m = 0.486;
J = 0.00383;
l = 0.25;
grav = 9.81;

%% Generate double-integrator trajectory

p1 = 1.5*[[2.;2.], [4.;1.8], [4.1;4.4], [3.;5.], [2.;4.]]+repmat([0.6;0],1,5);
p2 = 1.5*[[6.;3.], [7.;.4], [8.; 3.], [7.5;5.5], [5.9;4.]];
p3 = 1.5*[[5.;6.], [8.6;6.5], [8.6;7.5], [4.;8.]]-repmat([0.5;0],1,4);
p4 = 1.5*[[3.;3.], [3.5;3.2], [1.4;9.7], [1.;8.]];

load 'DI_safe.mat';

Xp = x_y_xd_yd_xdd_ydd(1:4,:)';
Up = x_y_xd_yd_xdd_ydd(5:6,:)';

T = (size(Xp,1)-1)*0.1;
t_traj_0 = (0:0.1:T)';

dt = 0.005;
t_traj = (0:dt:T)';
Xp = interp1(t_traj_0,Xp,t_traj);

%% Simulate

N = 20;
start_rand_p = [-0.05;-0.005] + 2*[0.1;0.01].*rand(2,N);
% start_rand_p = zeros(2,N);
start_rand_v = zeros(2,N);
x_0 =   [kron(ones(1,N),Xp(1,1:2)')+0.5*start_rand_p;
         kron(ones(1,N),Xp(1,3:4)')+start_rand_v;
         zeros(2,N)]; 

T_steps = length(t_traj);

%Store actual state
x = zeros(T_steps,6,N);
r = zeros(T_steps-1,6,N);

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
        if mod(i,100)==0
            disp(i);
        end
        
        xp = Xp(i,:)'; 
        
        r(i,:,k) = [state(1:4)-xp;state(5:6)]';
        
        u_safety = u_sol([r(i,:,k)']);
        V_bnds(i,k) = V_sol(r(i,:,k)');
        
        u = u_safety;
        
        u_act(i,:,k) = u';
        
        [d_t,d_x] = ode113(@(t,d_x)pvtol_sim(t,d_x,u),[t_traj(i),t_traj(i+1)],state,ode_options);
        
        state = d_x(end,:)';
        x(i+1,:,k) = state';
        state(5) = wrapToPi(state(5));
    end
end

%% Plot

close all; 

figure(1)
hold on
patch(p1(1,:),p1(2,:),'k','FaceAlpha',1);
patch(p2(1,:),p2(2,:),'k','FaceAlpha',1);
patch(p3(1,:),p3(2,:),'k','FaceAlpha',1);
patch(p4(1,:),p4(2,:),'k','FaceAlpha',1);
set(gcf,'Color','w');
for i = 1:round(0.04/dt):length(t_traj)
    proj_Ellipse([1:2],E_sol,1,[Xp(i,1);Xp(i,2)],30,'k');
end
plot(Xp(:,1),Xp(:,2),'r-','linewidth',2); 
grid off
hold on
for k = 1:N
    plot(x(:,1,k),x(:,2,k),'b-','linewidth',1.0); 
end
set(gca,'gridlinestyle','--')
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

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

%%

thrust_low = 0.1*m*grav;
thrust_up = m*grav;
u1 = @(t1,t2) (t1+t2)/m;
u2 = @(t1,t2) (t1-t2)*(l/J);
thrust_V = [u1(thrust_low,thrust_low),u2(thrust_low,thrust_low);
            u1(thrust_low,thrust_up),u2(thrust_low,thrust_up);
            u1(thrust_up,thrust_low),u2(thrust_up,thrust_low);
            u1(thrust_up,thrust_up),u2(thrust_up,thrust_up)];
thrust_poly = Polyhedron('V',thrust_V).minHRep();

figure()
thrust_poly.plot('color','b','alpha',0.3); hold on
for k = 1:N
    plot(u_act(:,1,k),u_act(:,2,k),'r-','linewidth',2);
end
grid on



