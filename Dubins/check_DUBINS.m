function E_sol = check_DUBINS(dynamics,r,rt,up,V,rho,u,gp,gs,toler)

delta_E = toler.delta_E;
n_c = length(dynamics.c);

%% check controller

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
for i = 1:gs.N_s
    prog = prog.withSOS(-(gs.A(i,:)*u - gs.b(i) - u_eps(i)) + Ls(i)*(Vr-rho));
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

%% Stability constraint

%Initialize problem
prog = spotsosprog;
prog = prog.withIndeterminate(r);
prog = prog.withIndeterminate(up);

%Finite variables
[prog, gamma] = prog.newFree(1);

%Multiplier fncs
ly_order = toler.all_deg.ly_order;
[prog, Ly] = prog.newFreePoly(monomials([r;up],0:ly_order),1);

lp_order = toler.all_deg.lp_order;
[prog, Lp] = prog.newSOSPoly(monomials([r;up],0:lp_order),length(gp));

%Stabilizer
Vr = subs(V,rt,dynamics.c);
dV_dr = subs(diff(V,rt),rt,dynamics.c)*dynamics.Jc;
prog = prog.withSOS(-dV_dr*dynamics.h - dV_dr*dynamics.B*u + gamma + Ly*(Vr-rho) + Lp'*gp);

%Strict feasibility
prog = prog.withPos(-gamma);

%Solver options
options = spot_sdp_default_options();
options.verbose = 0;

%Solve
fprintf('Solving Check (2)...');
try
    SOS_soln = prog.minimize(gamma , @spot_mosek, options);
catch
    %failed
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

%% Ellipsoidal bound

E_sol = zeros(n_c);

le_order = toler.all_deg.le_order; 

[prog, obj, E, ~] = define_C_3(r,rt,V,rho,dynamics.c,n_c,le_order,delta_E);
 
%Solver options
options = spot_sdp_default_options();
options.verbose = 0;

%Solve
fprintf('Solving Check (3)...');
try
    SOS_soln = prog.minimize(-0.1*obj, @spot_mosek, options);
catch
    return;
end

solved = SOS_soln.status.strcmp('STATUS_PRIMAL_AND_DUAL_FEASIBLE');

% Parse
if solved
	fprintf('solved \n');
    E_sol = double(SOS_soln.eval(E));
	clear SOS_soln prog;
else
    %failed
    fprintf('failed \n');
    return;
end

end