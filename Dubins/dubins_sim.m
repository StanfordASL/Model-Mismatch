function dx_dt = dubins_sim(t,x,u)

th = wrapToPi(x(4));

dx_dt = [x(3)*cos(th);
         x(3)*sin(th);
         u(1);
         u(2)];
end