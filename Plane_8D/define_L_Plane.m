function [prog, gamma, u_eps, E, V, rho, coeffs] = define_L_Plane(dynamics,r,rt,up,gp,gs,u,Ly,Ls,L_E,V_order,lp_order,n_r,n_c,delta_E,delta_rho)

global grav air_den area mass cd_0 Kd;

h = dynamics.h;
B = dynamics.B;
c = dynamics.c;
Jc = dynamics.Jc;

%% Initialize problem

prog = spotsosprog;
prog = prog.withIndeterminate(r);
prog = prog.withIndeterminate(rt);
prog = prog.withIndeterminate(up);

[prog, gamma] = prog.newFree(1);
[prog, u_eps] = prog.newPos(gs.N_s);
[prog, E] = prog.newSym(n_c);
[prog, rho] = prog.newPos(1);

%% Lyap. fnc

[prog, V, v_coeff] = prog.newFreePoly(monomials(rt,1:V_order),1);
Vr = subs(V,rt,c);
dV_dr = subs(diff(V,rt),rt,c)*Jc;

%% Define multipliers

[prog, Lp, lp_coeff] = prog.newSOSPoly(monomials([r;up],0:lp_order),length(gp));

coeffs = [lp_coeff];

%% Constraints

%Stability constraints
prog = prog.withSOS(-dV_dr*h - dV_dr*B*u + gamma + Ly*(Vr-rho) + Lp'*gp);

%Control constraints
%add in drag term for acceleration
F_drag = (air_den*area)*(r(5)^2)*(cd_0+4*(pi^2)*Kd*(r(8)^2));
% sin_g = 0.50586*(r(6)/(pi/6));
u_a = u(1) + (F_drag/mass);% + grav*sin_g;

for i = 1:gs.N_s
    prog = prog.withSOS(-(gs.A(i,:)*[u_a;u(2:3)] - gs.b(i) - u_eps(i)) + Ls(i)*(Vr-rho));
end

%Normalize V
prog = prog.withEqs(rho-1); %rho == 1

% E PSD
prog = prog.withPSD(E - delta_E*eye(n_c));

%Bounding ellipsoid
prog = prog.withSOS(1 - (rt'*E*rt) + L_E*(V-rho));

end