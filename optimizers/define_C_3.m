function [prog, obj, E, L_E] =  define_C_3(r,rt,V,rho,c,n_c,le_order,delta_E)

%% Initialize problem
prog = spotsosprog;
prog = prog.withIndeterminate(rt);

[prog, E] = prog.newSym(n_c);

%% Multiplier fncs

[prog, L_E] = prog.newSOSPoly(monomials(rt,0:le_order),1);

%% Constraints

% E PSD
prog = prog.withPSD(E - delta_E*eye(n_c));

%Bounding ellipsoid
prog = prog.withSOS( 1 - (rt'*E*rt) + L_E*(V-rho) );

%% Objective

[prog, obj] = maxdet(prog,E);

end