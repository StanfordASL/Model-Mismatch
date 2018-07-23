function [path, T, C] = steer(v1,v2,return_path)

global v_hat dT min_turn_r vz_p_lim;

%% X-Y dubins

[dubins_param, C_xy] = dubins_cost(v1([1;2;4]),v2([1;2;4]),min_turn_r);
T = C_xy/v_hat;

Dz_abs = abs(v2(3)-v1(3));
%figure out if delta_z reachable:
if (Dz_abs <= vz_p_lim*T)
    %diagonal Dubins path
    dz = (v2(3)-v1(3))*dT/T;
    Dz_reach = 1;
    
    ds = v_hat*dT;
    N_steps = C_xy/ds;
    C = N_steps*sqrt(ds^2 + dz^2);    
% elseif (Dz_abs >= T*vz_p_lim+(2*pi*min_turn_r/v_hat)*vz_p_lim)
%     %diagonal Dubins path + helix
%     T_extra = (Dz_abs - vz_p_lim*T)/vz_p_lim;
%     dz = vz_p_lim*dT;
%     Dz_reach = 2;
%     
%     ds = v_hat*dT;
%     N_steps_1 = C_xy/ds;
%     C_1 = N_steps_1*sqrt(ds^2 + dz^2); %part 1
%     
%     om_helix = 2*pi/(T_extra);
%     N_steps_2 = T_extra/dT;
%     C_2 = N_steps_2*sqrt(ds^2 + dz^2); %part 2
%     C = C_1 + C_2;    
else
    C = inf;
    path = 0;
    return;
end

if (return_path && Dz_reach)
    %find path
    path_xy = get_dubins_curve(dubins_param,ds);
    N_steps = size(path_xy,1);
    path_z = (linspace(v1(3),v2(3),N_steps))';
    path = [path_xy(:,1:2),path_z,path_xy(:,3)];
else
    path = 0;
end

end