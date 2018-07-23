function dx_dt = pvtol_sim(t,x,u)

global grav;

dx_dt = [x(3);
         x(4);
        -u(1)*sin(x(5));
         u(1)*cos(x(5))-grav;
         x(6);
         u(2)];
end