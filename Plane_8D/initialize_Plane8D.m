function [V_1, rho_1, E_1, u_sol] = initialize_Plane8D(dynamics,r,rt,up,gp,gs,g,toler)

%%
%states: r = (rx,ry,rth,v,om)

%% Initial guess

rho_1 = 1;
V_mon = rt;

%         x y z ps vhx  vhy  vv  g  ph al
% A_diag = (1./ ([1;1;1;pi/6;
%                 2;1;1;
%                 pi/6;pi/4;pi/5]-0.1)).^2;

A_diag = (1./ ([1;1;1;pi/6;2;pi/6;pi/4;pi/5]-0.1)).^2;
A_rand = -0.1+0.2*rand(length(rt)); A_rand = A_rand'*A_rand;
A = diag(A_diag) + A_rand;

V_1 = V_mon'*A*V_mon;

solve_1 = 0; %set to zero if loading a u_sol, V_sol,gamma_1 from memory

%% Loop
converged = 0;
itr = 1;
tol = toler.slack_tol;

fprintf('*******Starting initialization**********\n');
while (~converged)
    
    fprintf('Iteration: %d, \n',itr);
    
    if (solve_1)
        [solved_1,L,u,E_1,gamma_1] = solve_Plane8D_start_1(dynamics,r,rt,up,gp,gs,g,V_1,rho_1,toler);
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
    else
        load Plane1_solve.mat; 
        solve_1 = 1;
    end
    
    Ly = L_sol{1}; Lg = L_sol{2}; Ls = L_sol{3}; L_E = L_sol{4};
        
    fprintf('-------\n');
    [solved_2,V_2,rho_2,E_2,gamma_2,u_s] = solve_Plane8D_start_2(dynamics,r,rt,up,gp,gs,g,u_sol,Ly,Lg,Ls,L_E,gamma_1,toler);
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


save('Plane8D_V.mat','V_1','rho_1','E_1','u_sol');

fprintf('***** Done initialization ****** \n');

end