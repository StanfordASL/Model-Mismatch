function E_sol = check_PVTOL(dynamics,r,rt,up,gp,gs,g,V,rho,u,toler)

delta_E = toler.delta_E;
n_c = length(dynamics.c);

%% Most stabilizing controller

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
for i = 1:gs.N_s
    prog = prog.withSOS(-(gs.A(i,:)*u - gs.b(i) - u_eps(i)) + Ls(i)*(Vr-rho));
end

%Solver options
options = spot_sdp_default_options();
options.verbose = 2;

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

%% Stability constraint feasibility

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
options.verbose = 2;

%Solve
fprintf('Solving Check (2)...');
try
    SOS_soln = prog.minimize(gamma , @spot_mosek, options);
catch
    %failed
    fprintf('failed \n');
    return;
end

solved = SOS_soln.status.strcmp('STATUS_PRIMAL_AND_DUAL_FEASIBLE');

% Parse
if solved
	fprintf('solved, gamma: %.6f \n',double(SOS_soln.eval(gamma)));
    clear SOS_soln prog;
else
	%failed
    fprintf('failed \n');
    return;
end

%% Ellipsoidal bound

E_sol = zeros(n_c);

le_order = toler.all_deg.le_order; 
lg_order = toler.all_deg.lg_order;

[prog, obj, E, ~] = define_C_3(r,rt,V,rho,dynamics.c,n_c,le_order,delta_E);
 
%State constraints
[prog, Lg] = prog.newSOSPoly(monomials(rt,0:lg_order),1);
prog = prog.withSOS(-g + Lg*(V-rho));

%Solver options
options = spot_sdp_default_options();
options.verbose = 2;

%Solve
fprintf('Solving Check (3)...');
try
    SOS_soln = prog.minimize(-obj, @spot_mosek, options);
catch
    return;
end

solved = SOS_soln.status.strcmp('STATUS_PRIMAL_AND_DUAL_FEASIBLE');

% Parse
if solved
	fprintf('solved \n');
    E_sol = clean(double(SOS_soln.eval(E)),1e-7);
	clear SOS_soln prog;
else
    %failed
    fprintf('failed \n');
    return;
end

end