function [prog,gamma,Ly,Ls,u,coeffs] = define_K_1_Plane(dynamics,r,rt,up,gp,gs,V,rho,u_order,ls_order,ly_order,lp_order)

global air_den area mass cd_0 Kd;

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
[prog, u, u_coeff] = prog.newFreePoly(monomials(r,0:u_order),size(gs.A,2));

%% Multiplier fncs
[prog, Ls, ls_coeff] = prog.newSOSPoly(monomials(r,0:ls_order),gs.N_s);

[prog, Ly, ly_coeff] = prog.newFreePoly(monomials([r;up],0:ly_order),1);

[prog, Lp, lp_coeff] = prog.newSOSPoly(monomials([r;up],0:lp_order),length(gp));

coeffs = [u_coeff;ls_coeff;ly_coeff;lp_coeff];

%% Constraints

%Stability
dV_dr = subs(diff(V,rt),rt,c)*Jc;
Vr = subs(V,rt,c);

prog = prog.withSOS( -dV_dr*h - dV_dr*B*u + gamma + Ly*(Vr-rho) + Lp'*gp );

%Control constraints
%add in drag term for acceleration
F_drag = (air_den*area)*(r(5)^2)*(cd_0+4*(pi^2)*Kd*(r(8)^2));
u_a = u(1) + (F_drag/mass);

for i = 1:gs.N_s
    prog = prog.withSOS( -(gs.A(i,:)*[u_a;u(2:3)] - gs.b(i)) + Ls(i)*(Vr-rho) );
end

end