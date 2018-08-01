function [solved_1,L_sol,u_sol,E,gamma_sol] = solve_Plane8D_start_1(dynamics,r,rt,up,gp,gs,g,V_1,rho_1,toler)

delta_E = toler.delta_E;
n_c = length(dynamics.c);

%states: r = (x_r,y_r,th_r,v,om)

%Initialize outputs
solved_1 = 0;
L_sol = {};
u_sol = {};
E = 0;
gamma_sol = 0;

%% Problem 3

le_order = toler.all_deg.le_order; 
lg_order = toler.all_deg.lg_order;

[prog, obj, E, L_E] = define_K_2(r,rt,V_1,rho_1,dynamics.c,n_c,le_order,delta_E);
 
%State constraints
[prog, Lg] = prog.newSOSPoly(monomials(rt,0:lg_order),length(g));
for i = 1:length(g)
    prog = prog.withSOS( -g(i) + Lg(i)*(V_1-rho_1) );
end

%Solver options
options = spot_sdp_default_options();
options.verbose = 0;

%Solve
fprintf('Solving Controller (2)...');
try
    SOS_soln = prog.minimize(-obj, @spot_mosek, options);
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
    Lg = clean(SOS_soln.eval(Lg),clean_eps);
	E = clean(double(SOS_soln.eval(E)),clean_eps);
	clear SOS_soln prog;
else
    %failed
    fprintf('failed \n');
    return;
end

%% Problems 1 & 2 combined

u_order = toler.all_deg.u_order;
ls_order = toler.all_deg.ls_order;
ly_order = toler.all_deg.ly_order; 
lp_order = toler.all_deg.lp_order;

[prog,gamma,Ly,Ls,u,free_vars] = define_K_1_Plane(dynamics,r,rt,up,gp,gs,V_1,rho_1,u_order,ls_order,ly_order,lp_order);
% free_vars = [prog.coneVar; prog.freeVar];
len = length(free_vars);
[prog, reg] = prog.newPos(len);
prog = prog.withPos(reg-free_vars); %reg >= free_vars
prog = prog.withPos(free_vars + reg); %free_vars >= -reg

%Solver options
options = spot_sdp_default_options();
options.verbose = 2;

%Solve
fprintf('Solving Controller (1)...');
try
    SOS_soln = prog.minimize(5*gamma+(1e-2)*sum(reg), @spot_mosek, options);
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

%% Assemble final multipliers

solved_1 = 1;
L_sol = {Ly, Lg, Ls, L_E};

%% Save solution

gamma_1 = gamma_sol;
save('Plane1_solve.mat','V_1','rho_1','u_sol','L_sol','gamma_1');

end