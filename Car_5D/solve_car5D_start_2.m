function [solved_2,V_sol,rho,E,gamma_sol,u_s] = solve_car5D_start_2(dynamics,r,rt,up,gp,gs,g,u,Ly,Lg,Ls,L_E,gamma_1,toler)

delta_E = toler.delta_E;
 
n_r = 5;
n_c = length(dynamics.c);

%% Define Problem

V_order = toler.all_deg.V_order; 
lp_order = toler.all_deg.lp_order;

[prog, gamma, u_eps, E, V, rho] = define_V(dynamics,r,rt,up,gp,gs,u,Ly,Ls,L_E,V_order,lp_order,n_r,n_c,delta_E);

%slack improvement
prog = prog.withPos(gamma_1-gamma);

%state-space constraints
prog = prog.withSOS(-g + Lg*(V-rho));

%% Solve
options = spot_sdp_default_options();
options.verbose = 0;

fprintf('Solving V...');
try
    SOS_soln = prog.minimize(gamma+sum(u_eps), @spot_mosek, options);
catch
    %failed
    solved_2 = 0;
    V_sol = 0;
    rho = 0;
    E = 0;
    gamma_sol = 10;
    u_s = 10;
    return;
end

solved = SOS_soln.status.strcmp('STATUS_PRIMAL_AND_DUAL_FEASIBLE');

%% Parse
if solved
    fprintf('solved \n');
    clean_eps = 1e-9;
    solved_2 = 1;
    V_sol = clean(SOS_soln.eval(V),clean_eps);
    rho = double(SOS_soln.eval(rho));
    E = clean(double(SOS_soln.eval(E)),clean_eps);
    gamma_sol = double(SOS_soln.eval(gamma));
    u_s = max(double(SOS_soln.eval(u_eps)));
    
    clear SOS_soln prog;
else
    %failed
    fprintf('failed \n');
    solved_2 = 0;
    V_sol = 0;
    rho = 0;
    E = 0;
    gamma_sol = 10;
    u_s = 10;
    return;
end

end