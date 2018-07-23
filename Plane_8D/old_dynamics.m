dynamics.h = [r(5)*cos_g*cos_ps-v_hat+up(1)*r(2); %5
    r(5)*cos_g*sin_ps-up(1)*r(1); %4
    r(5)*sin_g-up(2); %2
    -C_con*r(8)*r(5)*sin_ph*sec_g-up(1); %5
    -(F_drag/mass)-grav*sin_g; %4
    C_con*r(5)*r(8)*cos_ph-grav*cos_g*v_inv; %4
    0;
    0];a

dynamics.B = [zeros(4,3);
    1,0,0;
    0,0,0;
    0,1,0;
    0,0,1];

dynamics.c = [r(1:4);
              r(5)-v_hat;
              r(6:7);
              r(8)-alpha_0];
dynamics.Jc = eye(8);

%%

%bounded state vector          
dynamics.c = [r(1:4);
              r(5)*cos_g*cos_ps-v_hat;
              r(5)*cos_g*sin_ps;
              r(5)*sin_g;
              r(6);
              r(7);
              r(8)-alpha_0];
          
%dynamics of rt_dot, so dV_dt = dV_drt * rt_dot

dynamics.h = [...
r(5)*cos_g*cos_ps-v_hat+up(1)*r(2); %xr_dot 
r(5)*cos_g*sin_ps-up(1)*r(1);  %yr_do
r(5)*sin_g-up(2); %zr_dot
-C_con*r(8)*r(5)*sin_ph*sec_g-up(1); %psi_r_dot
-C_con*r(8)*(r(5)^2)*(cos_ph*cos_ps+sin_ph*sin_ps)+r(5)*cos_g*sin_ps*up(1)-...
                                                   cos_g*cos_ps*(F_drag/mass); %v_par_dot
C_con*r(8)*(r(5)^2)*(sin_ph*cos_ps-cos_ph*sin_ps)-r(5)*cos_g*cos_ps*up(1)-...
                                                   cos_g*sin_ps*(F_drag/mass); %v_perp_dot
C_con*r(8)*(r(5)^2)*cos_g*cos_ph-sin_g*(F_drag/mass)-grav; %v_vert_dot
C_con*r(5)*r(8)*cos_ph-grav*cos_g*v_inv; %gamma_dot
0;
0];

dynamics.B = [...
zeros(4,3);
cos_g*cos_ps,0,0;
cos_g*sin_ps,0,0;
sin_g,0,0;
0,0,0;
0,1,0;
0,0,1];

g = [rt(4)^2-ps_max^2;
     rt(5)^2+rt(6)^2+rt(7)^2-bma_v^2;
     rt(8)^2-g_max^2;
     rt(9)^2-ph_max^2];