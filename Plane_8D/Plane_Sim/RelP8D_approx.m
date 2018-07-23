function dr_dt = RelP8D_approx(t,r,u,up)
global grav air_den area mass cd_0 Kd v_hat;

C_con = pi*air_den*area/mass;

%% Approximation made for SOS

sin_ps = 0.50586*(r(4)/(pi/6));
cos_ps = 0.9326 - 0.066986*(2*(r(4)/(pi/6))^2 -1);

sin_g = 0.50586*(r(6)/(pi/6));
cos_g = 0.9326 - 0.066986*(2*(r(6)/(pi/6))^2 -1);
sec_g = 1.07502 + 0.07728*(2*(r(6)/(pi/6))^2 -1);

sin_ph = 0.7264*(r(7)/(pi/4));
cos_ph = 0.8516 - 0.1464*(2*(r(7)/(pi/4))^2 -1);

%for v_hat ~ 6.5
v_min = 3; v_max = 10;
bma_v=0.5*(v_max-v_min);
bpa_v=0.5*(v_min+v_max);
v_inv = 0.182574 -0.1067*(r(5)-bpa_v)/bma_v + ...
        0.031181*(2*((r(5)-bpa_v)/bma_v)^2 -1) +...
       -0.00911*(4*((r(5)-bpa_v)/bma_v)^3 - 3*((r(5)-bpa_v)/bma_v));

% Dynamics

F_drag = (air_den*area)*(r(5)^2)*(cd_0+4*(pi^2)*Kd*(r(8)^2));

h = [r(5)*cos_g*cos_ps-v_hat+up(1)*r(2); %5
              r(5)*cos_g*sin_ps-up(1)*r(1); %4
              r(5)*sin_g-up(2); %2
              -C_con*r(8)*r(5)*sin_ph*sec_g-up(1); %5
              -grav*sin_g; %4
              C_con*r(5)*r(8)*cos_ph-grav*cos_g*v_inv; %4
              0;
              0];

B = [zeros(4,3);
              1,0,0;
              0,0,0;
              0,1,0;
              0,0,1];
      
dr_dt = h + B*u; 

end