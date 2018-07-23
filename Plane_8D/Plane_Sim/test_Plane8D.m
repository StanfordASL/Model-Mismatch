%% dynamics params

global grav air_den area mass cd_0 Kd v_hat;

%% Setup smooth controller

smooth_cntrl = setup_Plane_cntrl();

%% Find some nominal trajectory

xp_f = [10,10,7,pi/2]';
[X_ref, ~, ~] = dubins_curve(xp_0([1;2;4]),xp_f([1;2;4]),min_turn_r,v_hat*dT,1);
T = (length(X_ref)-1)*dT;
t_ref = 0:dT:T;

X_ref(:,4) = wrapToPi(X_ref(:,3));

idx_mid = round(0.5*T/dT);
X_ref(1:idx_mid,3) = (linspace(xp_0(3),xp_0(3)+vz_p_lim*(T/2),idx_mid))';
z_mid = X_ref(idx_mid,3);
X_ref(idx_mid+1:end,3) = (linspace(z_mid,xp_0(3),length(t_ref)-idx_mid))';

fig_FMT = figure();

%% Simulate

simulate_Plane8D;




