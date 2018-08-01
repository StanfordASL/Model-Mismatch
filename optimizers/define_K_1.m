function [prog,gamma,Ly,Ls,u] = define_K_1(dynamics,r,rt,up,gp,gs,V,rho,u_order,ls_order,ly_order,lp_order)

h = dynamics.h;
B = dynamics.B;
c = dynamics.c;
Jc = dynamics.Jc;

%% Initialize problem
prog = spotsosprog;
prog = prog.withIndeterminate(r);
prog = prog.withIndeterminate(up);

%% Finite variables
[prog, gamma] = prog.newFree(1);

%% Controller fnc
[prog, u] = prog.newFreePoly(monomials(r,0:u_order),size(gs.A,2));

%% Multiplier fncs
[prog, Ls] = prog.newSOSPoly(monomials(r,0:ls_order),gs.N_s);

[prog, Ly] = prog.newFreePoly(monomials([r;up],0:ly_order),1);

[prog, Lp] = prog.newSOSPoly(monomials([r;up],0:lp_order),length(gp));

%% Constraints

%Stability
dV_dr = subs(diff(V,rt),rt,c)*Jc;
Vr = subs(V,rt,c);

prog = prog.withSOS(-dV_dr*h - dV_dr*B*u + gamma + Ly*(Vr-rho) + Lp'*gp);

%Control constraints
for i = 1:gs.N_s
    prog = prog.withSOS(-(gs.A(i,:)*u - gs.b(i)) + Ls(i)*(Vr-rho));
end

end