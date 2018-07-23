function dxp_dt = DubinsP4D_sim(t,xp,up)

global v_hat;

psi_p = wrapToPi(xp(4));

dxp_dt = [v_hat*cos(psi_p);
          v_hat*sin(psi_p);
          up(2);
          up(1)];
      
end
         
   