function E_sol = check_Plane8D(dynamics,r,rt,up,gp,gs,g,V,rho,u,toler)

global grav air_den area mass cd_0 Kd v_hat;

delta_E = toler.delta_E;
n_c = length(dynamics.c);
E_sol = zeros(n_c);

%% Ellipsoidal bound

le_order = toler.all_deg.le_order; 
lg_order = toler.all_deg.lg_order;

[prog, obj, E, ~] = define_C_3(r,rt,V,rho,dynamics.c,n_c,le_order,delta_E);
 
%State constraints
[prog, Lg] = prog.newSOSPoly(monomials(rt,0:lg_order),length(g));
for i = 1:length(g)
    prog = prog.withSOS(-g(i) + Lg(i)*(V-rho));
end

%Solver options
options = spot_sdp_default_options();
options.verbose = 0;

%Solve
fprintf('Solving Check (3)...');
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
    clean_eps = 1e-7;
	E_sol = clean(double(SOS_soln.eval(E)),clean_eps);
	clear SOS_soln prog;
else
    %failed
    fprintf('failed \n');
    return;
end

%% Stability constraint feasibility

%Initialize problem
prog = spotsosprog;
prog = prog.withIndeterminate(r);
prog = prog.withIndeterminate(up);

%Finite variables
[prog, gamma] = prog.newFree(1);

%Multiplier fncs
ly_order = toler.all_deg.ly_order;
[prog, Ly, ly_coeff] = prog.newFreePoly(monomials([r;up],0:ly_order),1);

lp_order = toler.all_deg.lp_order;
[prog, Lp, lp_coeff] = prog.newSOSPoly(monomials([r;up],0:lp_order),length(gp));

prog = prog.withPos(-0.05-gamma);

%Stabilizer
Vr = subs(V,rt,dynamics.c);
dV_dr = subs(diff(V,rt),rt,dynamics.c)*dynamics.Jc;
prog = prog.withSOS(-dV_dr*dynamics.h - dV_dr*dynamics.B*u + gamma + Ly*(Vr-rho) + Lp'*gp);

free_vars = [ly_coeff;lp_coeff];
len = length(free_vars);
[prog, reg] = prog.newPos(len);
prog = prog.withPos(-free_vars + reg); %reg >= free_vars
prog = prog.withPos(free_vars + reg); %free_vars >= -reg

%Solver options
options = spot_sdp_default_options();
options.verbose = 2;

%Solve
fprintf('Solving Check (2)...');
try
    SOS_soln = prog.minimize(gamma + (1e-2)*sum(reg), @spot_mosek, options);
catch
    %failed
    fprintf('failed \n');
    return;
end

solved = SOS_soln.status.strcmp('STATUS_PRIMAL_AND_DUAL_FEASIBLE');

% Parse
if solved
	fprintf('solved, gamma: %.6f \n',double(SOS_soln.eval(gamma)));
else
	%failed
    fprintf('failed \n');
    return;
end

%% Check control constraints

%Initialize problem
prog = spotsosprog;
prog = prog.withIndeterminate(r);
prog = prog.withIndeterminate(up);

%Finite variables
[prog, u_eps] = prog.newPos(gs.N_s);

%Multiplier fncs
ls_order = toler.all_deg.ls_order;
[prog, Ls] = prog.newSOSPoly(monomials(r,0:ls_order),gs.N_s);

%Control constraints
Vr = subs(V,rt,dynamics.c);
F_drag = (air_den*area)*(r(5)^2)*(cd_0+4*(pi^2)*Kd*(r(8)^2));
% sin_g = 0.50586*(r(6)/(pi/6));
u_a = u(1) + (F_drag/mass);%

for i = 1:gs.N_s
    prog = prog.withSOS(-(gs.A(i,:)*[u_a;u(2:3)] - gs.b(i) - u_eps(i)) + Ls(i)*(Vr-rho));
end

%Solver options
options = spot_sdp_default_options();
options.verbose = 0;

%Solve
fprintf('Solving Check (1)...');
try
    SOS_soln = prog.minimize(sum(u_eps), @spot_mosek, options);
catch
    %failed
    return;
end

solved = SOS_soln.status.strcmp('STATUS_PRIMAL_AND_DUAL_FEASIBLE');

% Parse
if solved
	fprintf('solved, u_s: %.6f \n',max(double(SOS_soln.eval(u_eps))));	
	clear SOS_soln prog;
	
else
    %failed
    fprintf('failed \n');
    return;
end

end