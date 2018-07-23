function [dr_dt, h] = RelP8D_sim(t,r,u,up)

global grav air_den area mass cd_0 Kd v_hat;

%% True relative dynamics

v = r(5);

psi_r = wrapToPi(r(4));
gamma = wrapToPi(r(6));
phi = wrapToPi(r(7));
alpha = wrapToPi(r(8));

F_lift = pi*air_den*area*alpha*(v^2);

F_drag = air_den*area*(v^2)*(cd_0 + 4*(pi^2)*Kd*(alpha^2) );

h = [v*cos(gamma)*cos(psi_r)-v_hat+up(1)*r(2);
         v*cos(gamma)*sin(psi_r)-up(1)*r(1);
         v*sin(gamma)-up(2);
         -(F_lift*sin(phi))/(mass*v*cos(gamma))-up(1);
         -grav*sin(gamma);
         (F_lift*cos(phi))/(mass*v)-(grav*cos(gamma)/v);
          0;0];

dr_dt = [v*cos(gamma)*cos(psi_r)-v_hat+up(1)*r(2);
         v*cos(gamma)*sin(psi_r)-up(1)*r(1);
         v*sin(gamma)-up(2);
         -(F_lift*sin(phi))/(mass*v*cos(gamma))-up(1);
         -grav*sin(gamma)+u(1);
         (F_lift*cos(phi))/(mass*v)-(grav*cos(gamma)/v);
          u(2);
          u(3)];

end
         
   