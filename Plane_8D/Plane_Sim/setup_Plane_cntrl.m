function Prob  = setup_Plane_cntrl()

yalmip('clear');

%% variables

u = sdpvar(3,1);
dVdr_h = sdpvar(1);
dVdr_B = sdpvar(3,1);

%% Constraints

Constraints = [dVdr_h + dVdr_B'*u <= 0];

%% Define problem

Objective = u'*u;

Prob = optimizer(Constraints,Objective,sdpsettings('solver','mosek'),{dVdr_h,dVdr_B},u);


end