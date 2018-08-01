function [solved_1,L,u_sol,E_sol,gamma_sol] = solve_DUBINS_1(dynamics,r,rt,up,gp,gs,V,rho,E_w,toler)

delta_E = toler.delta_E;
n_c = length(dynamics.c);

%Initialize outputs
solved_1 = 0;
L = {};
u_sol = {};
E_sol = 0;
gamma_sol = 0;

%% Problems 1 & 2 combined

u_order = toler.all_deg.u_order;
ls_order = toler.all_deg.ls_order;
ly_order = toler.all_deg.ly_order; 
lp_order = toler.all_deg.lp_order;

[prog,gamma,Ly,Ls,u] = define_K_1(dynamics,r,rt,up,gp,gs,V,rho,u_order,ls_order,ly_order,lp_order);

%Strict feasibility
prog = prog.withPos(toler.slack_tol-gamma);

%Solver options
options = spot_sdp_default_options();
options.verbose = 0;

%Solve
fprintf('Solving Controller (1)...');
try
    SOS_soln = prog.minimize(gamma, @spot_mosek, options);
catch
    %failed
    fprintf('failed \n');
    return;
end

solved = SOS_soln.status.strcmp('STATUS_PRIMAL_AND_DUAL_FEASIBLE');

% Parse
if solved
	fprintf('solved \n');
    clean_eps = 1e-8;
    gamma_sol = clean(double(SOS_soln.eval(gamma)),clean_eps);
    Ly = clean(SOS_soln.eval(Ly),clean_eps);
    u_sol = clean(SOS_soln.eval(u),clean_eps);
	Ls = clean(SOS_soln.eval(Ls),clean_eps);
    
	clear prog;
else
	%failed
    fprintf('failed \n');
    return;
end

%% Problem 3

le_order = toler.all_deg.le_order;

[prog, obj, E, L_E] = define_K_2(r,rt,V,rho,dynamics.c,n_c,le_order,delta_E);

%Solver options
options = spot_sdp_default_options();
options.verbose = 0;

%Solve
fprintf('Solving Controller (2)...');
try
    SOS_soln = prog.minimize(-0.1*obj, @spot_mosek, options);
catch
    fprintf('failed \n');
    return;
end

solved = SOS_soln.status.strcmp('STATUS_PRIMAL_AND_DUAL_FEASIBLE');

% Parse
if solved
	fprintf('solved \n');
    clean_eps = 1e-8;
    L_E = clean(SOS_soln.eval(L_E),clean_eps);
	E_sol = clean(double(SOS_soln.eval(E)),clean_eps);
	clear SOS_soln prog;
else
    %failed
    fprintf('failed \n');
    return;
end

%% Assemble final multipliers

solved_1 = 1;
L = {Ly, Ls, L_E};
end