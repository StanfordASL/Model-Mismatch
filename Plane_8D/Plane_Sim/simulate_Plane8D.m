disp('Ready to Simulate');
keyboard;

%% Data store

T_steps = length(t_ref);

%Store actual state
X = zeros(T_steps,8);
R = zeros(T_steps-1,8);

%Store control history
U = zeros(T_steps-1,3);
Up = zeros(T_steps-1,2);

V_bnds = zeros(T_steps-1,1);

X(1,:) = x_0';
x = x_0;

dVdt = zeros(T_steps-1,2);

%% Simulate

ode_options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

for i = 1:T_steps-1
    
    up_1 = wrapToPi(X_ref(i+1,4)-X_ref(i,4))/dT;
    up_2 = (X_ref(i+1,3)-X_ref(i,3))/dT;
    Up(i,:) = [up_1,up_2];
    
    %% Actual plane
    
    xp = X_ref(i+1,:)';
    
    rot = [cos(xp(4)), sin(xp(4));
        -sin(xp(4)), cos(xp(4))];
    
    r =   [rot*(x(1:2)-xp(1:2));
           x(3)-xp(3);
           wrapToPi(x(4)-xp(4));
           x(5);
           x(6:8)];
    
    R(i,:) = r';
    
    V_bnds(i) = V_sol(c_sol(r));
    
    u_fb = u_sol(r);
    F_drag = (air_den*area)*(r(5)^2)*(cd_0+4*(pi^2)*Kd*(r(8)^2));
    u_a = u_fb(1) + (F_drag/mass); %add in drag ffwd for accel
    u = [u_a; u_fb(2:3)];
    
    U(i,:) = u';
    
    %% Propagate true dynamics
    
    [~,d_x] = ode113(@(t,d_x)Plane8D_sim(t,d_x,u),[t_ref(i),t_ref(i+1)],x,ode_options);
    
    %% Check Error in relative dynamics
    
    dVdr = dVsol_dr(r);
    [r_dot,~] = RelP8D_sim(0,r,u_fb,Up(i,:));
    r_dot_approx = RelP8D_approx(0,r,u_fb,Up(i,:)');
    dVdt(i,1) = dVdr'*r_dot;
    dVdt(i,2) = dVdr'*r_dot_approx;
    
    %% Update
    
    x = d_x(end,:)';
    x([4;6;7;8]) = wrapToPi(x([4;6;7;8]));
    X(i+1,:) = x';
    
end

%% Plots

plot_Plane8D;