function dx_dt = Plane8D_sim(t,x,u)

global grav air_den area mass cd_0 Kd;

v = x(5);

psi = wrapToPi(x(4));
gamma = wrapToPi(x(6));
phi = wrapToPi(x(7));
alpha = wrapToPi(x(8));

F_lift = pi*air_den*area*alpha*(v^2);

F_drag = air_den*area*(v^2)*(cd_0 + 4*(pi^2)*Kd*(alpha^2) );


dx_dt = [v*cos(gamma)*cos(psi);
         v*cos(gamma)*sin(psi);
         v*sin(gamma);
         -(F_lift*sin(phi))/(mass*v*cos(gamma));
         -(F_drag/mass)-grav*sin(gamma)+u(1);
         (F_lift*cos(phi))/(mass*v)-(grav*cos(gamma)/v);
          u(2);
          u(3)];
      
end
         
   