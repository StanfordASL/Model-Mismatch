function dx_dt = car5D_sim(t,x,u)

theta = wrapToPi(x(3));
dx_dt = [x(4)*cos(theta);
         x(4)*sin(theta);
         x(5);
         u(1);
         u(2)];
end