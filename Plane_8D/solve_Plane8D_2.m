function [solved_2,E_sol,V_sol,rho_sol,gamma_sol,u_s] = solve_Plane8D_2(dynamics,r,rt,up,gp,gs,g,u,Ly,Lg,Ls,L_E,E_prev,alpha,E_w,toler)

delta_E = toler.delta_E;
lambda = toler.lambda;
alpha_l = toler.alpha_l;

n_r = length(r);
n_c = length(dynamics.c);
 
%% Define Problem

V_order = toler.all_deg.V_order; 
lp_order = toler.all_deg.lp_order;

[prog, gamma, u_eps, E, V, rho, free_vars] = define_V_Plane(dynamics,r,rt,up,gp,gs,u,Ly,Ls,L_E,V_order,lp_order,n_r,n_c,delta_E);
% free_vars = [prog.coneVar; prog.freeVar];
len = length(free_vars);
[prog, reg] = prog.newPos(len);
prog = prog.withPos(-free_vars + reg); %reg >= free_vars
prog = prog.withPos(free_vars + reg); %free_vars >= -reg

%state-space constraints
for i = 1:length(g)
    prog = prog.withSOS(-g(i) + Lg(i)*(V-rho));
end

% E trust region
prog = prog.withPSD(E - alpha_l*E_prev);
prog = prog.withPSD((1+alpha)*E_prev - E);

%Objective epigraph
[prog, obj] = maxdet(prog,E);

%% Solve
options = spot_sdp_default_options();
options.verbose = 0;

fprintf('Solving V...');
try
    SOS_soln = prog.minimize(3*(-obj+lambda*(gamma+sum(u_eps))) + (1e-2)*sum(reg), @spot_mosek, options);
catch
    %failed
    solved_2 = 0;
    V_sol = 0;
    rho_sol = 0;
    E_sol = 0;
    gamma_sol = 10;
    u_s = 10;
    return;
end

solved = SOS_soln.status.strcmp('STATUS_PRIMAL_AND_DUAL_FEASIBLE');

%% Parse
if solved
    fprintf('solved \n');
    clean_eps = 0;
    solved_2 = 1;
    V_sol = clean(SOS_soln.eval(V),clean_eps);
    rho_sol = double(SOS_soln.eval(rho));
    E_sol = clean(double(SOS_soln.eval(E)),clean_eps);
    gamma_sol = double(SOS_soln.eval(gamma));
    u_s = max(double(SOS_soln.eval(u_eps)));
    
    clear SOS_soln prog;
else
    %failed
    fprintf('failed \n');
    solved_2 = 0;
    V_sol = 0;
    rho_sol = 0;
    E_sol = 0;
    gamma_sol = 10;
    u_s = 10;
    return;
end

%% Save backup copy

u_sol = u;
save('Plane8D_raw_2.mat','V_sol','rho_sol','E_sol','u_sol');

%%

%lambda = 0.5;
% E trust region
% prog = prog.withPSD(E - alpha_l*E_prev);
% prog = prog.withPSD((1+alpha)*E_prev - E);

end