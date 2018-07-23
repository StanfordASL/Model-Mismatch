function [V_1, rho_1, E_1, u_sol] = initialize_PVTOL(dynamics,r,rt,up,gp,gs,g,toler)

%%
%states: r = (rx,ry,vx,vy,phi,phi_dot)

global grav;

%% Initial guess
% A = [0,0,1,0,0,0;
%      0,0,0,1,0,0;
%      0,0,0,0,-grav,0;
%      zeros(1,6);
%      zeros(1,5),1;
%      zeros(1,6)];
% 
% B = [zeros(3,2);
%      1,0;
%      0,0;
%      0,1];
% 
% [~,E_0,~] = lqr(A,B,diag([10,10,10,10,0.2,0.1]),1.0*eye(2)); 
% rho_1 = 1.0;
% V_1 = 5*rt'*E_0*rt;

load PVTOL_V.mat;

%% Loop
converged = 0;
itr = 1;
tol = toler.slack_tol;

fprintf('*******Starting initialization**********\n');
while (~converged)
    
    fprintf('Iteration: %d, \n',itr);
    
    [solved_1,L,u,E_1,gamma_1] = solve_PVTOL_start_1(dynamics,r,rt,up,gp,gs,g,V_1,rho_1,toler);
    fprintf('Step 1: gamma: %.6f \n',gamma_1);
    if (solved_1)
        L_sol = L;
        u_sol = u;
        if (gamma_1 <= tol)
            V_sol = V_1; rho_sol = rho_1; E_sol = E_1;
            converged = 1;
            continue;
        end
    else
        keyboard;
    end
    
    Ly = L_sol{1}; Lg = L_sol{2}; Ls = L_sol{3}; L_E = L_sol{4};
        
    fprintf('-------\n');
    [solved_2,V_2,rho_2,E_2,gamma_2,u_s] = solve_PVTOL_start_2(dynamics,r,rt,up,gp,gs,g,u_sol,Ly,Lg,Ls,L_E,gamma_1,toler);
    fprintf('Step 2: (gamma, u_s): (%6f:%6f),\n ',gamma_2,u_s);
    
    if (solved_2)
        if (max(gamma_2,u_s) <= tol)
            converged = 1;
            V_sol = V_2; rho_sol = rho_2; E_sol = E_2;
        else
            V_1 = V_2; rho_1 = rho_2;
            itr = itr + 1;
        end
    else
        keyboard;
    end
    
    fprintf('*******************\n');
end

V_1 = V_sol;
rho_1 = rho_sol;
E_1 = E_sol;


save('PVTOL_V.mat','V_1','rho_1','E_1','u_sol');

fprintf('***** Done initialization ****** \n');

end