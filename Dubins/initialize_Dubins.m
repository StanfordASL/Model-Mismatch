function [V_1,rho_1,E_1,u_sol] = initialize_Dubins(dynamics,r,rt,up,gp,gs,toler)

%%
%states: r = (rx,ry,v,th)

%% Initial guess
rho_s = 6*rand(1);
V_mon = [rt;rt.^2];
A = randn(length(V_mon)); A = (A'*A)/rho_s;
% rho_1 = 1.0;
% V_1 = V_mon'*A*V_mon;

load('Dubins_V.mat','V_1','rho_1');
V_1 = V_1 + V_mon'*A*V_mon;

%% Loop
converged = 0;
itr = 1;
tol = toler.slack_tol;

fprintf('*******Begin Initialization**********\n');
while (~converged)
    
    fprintf('Iteration: %d, \n',itr);
    
    [solved_1,L,u,E_1,gamma_1] = solve_DUBINS_start_1(dynamics,r,rt,up,gp,gs,V_1,rho_1,toler);
    fprintf('Step 1: gamma: %.6f \n',gamma_1);
    if (solved_1)
        L_sol = L;
        u_sol = u;
        if ( gamma_1 <= tol)
            V_sol = V_1; rho_sol = rho_1; E_sol = E_1;
            converged = 1;
            continue;
        end
    else
        keyboard;
    end
    
    Ly = L_sol{1}; Ls = L_sol{2}; L_E = L_sol{3};
        
    fprintf('-------\n');
    [solved_2,V_2,rho_2,E_2,gamma_2,u_s] = solve_DUBINS_start_2(dynamics,r,rt,up,gp,gs,u_sol,Ly,Ls,L_E,gamma_1,toler);
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
    fprintf('****************\n');
end

V_1 = V_sol;
rho_1 = rho_sol;
E_1 = E_sol;


save('Dubins_V.mat','V_1','rho_1','E_1','u_sol');

fprintf('***** Done initialization ****** \n');

end