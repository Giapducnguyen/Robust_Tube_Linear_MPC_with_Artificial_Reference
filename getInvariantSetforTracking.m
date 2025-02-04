% % Function helper: Calculation of invariant set for tracking
% % Reference:
% % "Robust tube-based MPC for tracking of constrained linear systems with
% %  additive disturbances" - Subsection 3.4       < D. Limon et al. >

function set_Xlambda_a = getInvariantSetforTracking(method, iter_max, ...
                           lambda, A, B, Kbar, M_theta, set_Xbar, set_Ubar)

set_Zbar = set_Xbar*set_Ubar;
[nx,nu] = size(B);
L = [-Kbar, eye(nu)]*M_theta;
A_a = [A+B*Kbar,        B*L;
       zeros(nu, nx),   eye(nu)];
if method == 1
    % % Method 1 - Direct computation with tightened constraints:
    
    % % Set constraint on the extended state x_a = (x, theta)
    A_Xlambda = [set_Zbar.A*[eye(nx), zeros(nx,nu); Kbar, L];
                 [zeros(size(set_Zbar.A,1),nx), set_Zbar.A*M_theta]];
    b_Xlambda = [set_Zbar.b; lambda.*set_Zbar.b];
    set_Xlambda = Polyhedron('A', A_Xlambda, 'b', b_Xlambda);
    set_Xlambda.minHRep();
    
    % % Computation of the maximal admissible invariant set X_{\lambda}^a
    found = false; % flag
    k = 0; set_Xk = set_Xlambda; % initialization
    
    figure; hold on; grid on; box on;
    title("Computation Evolution of $\Omega_{t,\bar{K}} = \mathrm{Proj}_x(\mathcal{X}_{\lambda}^{a})$","Interpreter","latex");
    plot(set_Xk.projection(1:nx), 'color', 'black', 'alpha', 0.01, 'edgecolor', 'black');
    
    fprintf("\nComputing the invariant set for tracking: \n");
    while (~found && k < iter_max)
        fprintf("Iteration: %d \n", k+1);
    
        set_Temp = Polyhedron('A', set_Xlambda.A * A_a^(k+1), 'b', set_Xlambda.b);
        set_Temp.minHRep();
        set_Xkp1 = intersect(set_Xk, set_Temp);
        set_Xkp1.minHRep();
        
        pause(0.10);
    
        if set_Xkp1 == set_Xk
            found = true;
            pltf = plot(set_Xkp1.projection(1:nx), 'color', 'green', 'alpha', 0.01, 'edgecolor', 'green');
        else
            k = k + 1;
            set_Xk = set_Xkp1;
            plot(set_Xkp1.projection(1:nx), 'color', 'black', 'alpha', 0.01, 'edgecolor', 'black');
        end
    end
    legend(pltf,"$\Omega_{t,\bar{K}} = \mathrm{Proj}_x(\mathcal{X}_{\lambda}^{a})$","Interpreter","latex");
    xlabel("$x_1$","Interpreter","latex");
    ylabel("$x_2$","Interpreter","latex");
    hold off;
    set_Xlambda_a = set_Xkp1;
    
else
    % % Method 2 - Indirect computation with transformation and additional constraints
    
    % % following the standard formulation of the paper:
    % % "Linear systems with state and control constraints: the theory 
    % %  and application of maximal output admissible sets." - E. G. Gilbert
    
    %
    % Step 1: System matrix transformation - Equations (2.10 & 4.1) - E. G. Gilbert
    U = [eye(nx), (A+B*Kbar - eye(nx))^-1*(-B*L); zeros(nu), eye(nu)];
    A_aU = U^-1*A_a*U;
    A_aS = A_aU(1:nx, 1:nx);
    C_aU = eye(nx+nu)*U;
    C_aS = C_aU(:, 1:nx);
    C_aL = C_aU(:, nx+1: end);
    
    % Step 2: Set constraint on the extended state (x, theta) - Equation (5.4) - E. G. Gilbert
    A_W1 = [set_Zbar.A*[eye(nx), zeros(nx,nu); Kbar, L];
            set_Zbar.A*[zeros(size(M_theta,1),nx), M_theta]];
    b_W1 = [set_Zbar.b; set_Zbar.b];
    set_W1 = Polyhedron('A', A_W1, 'b', b_W1);
    set_W1.minHRep();
    
    A_Xlambda = [set_Zbar.A*[eye(nx), zeros(nx,nu); Kbar, L];
                 set_Zbar.A*[zeros(size(M_theta,1),nx), M_theta]];
    b_Xlambda = [set_Zbar.b; lambda.*set_Zbar.b];
    set_Xlambda = Polyhedron('A', A_Xlambda, 'b', b_Xlambda);
    set_Xlambda.minHRep();
    
    % Step 3: Compute the first set for initialization - Equation (5.4) - E. G. Gilbert
    set_X0_A = [set_W1.A*[C_aS, C_aL]; 
                set_Xlambda.A*[zeros(size(C_aL,1),nx), C_aL]]; 
    set_X0_b = [set_W1.b; set_Xlambda.b];
    set_X0 = Polyhedron('A', set_X0_A, 'b', set_X0_b);
    set_X0.minHRep();
    
    % figure; hold on; grid on;
    % title("Computation Evolution of $\Omega_{t,\bar{K}} = \mathrm{Proj}_x(\mathcal{X}_{\lambda}^{a})$","Interpreter","latex");
    % plot(set_X0.projection(1:nx), 'color', 'black', 'alpha', 0.01, 'edgecolor', 'black');
    
    % Step 4: Computation of the maximal admissible invariant set 
    % - Equation (5.3) - E. G. Gilbert
    k = 0; set_Xk = set_X0;
    found = false;
    fprintf("\nComputing the invariant set for tracking: \n");
    while(~found && k < iter_max)
        fprintf('Iteration: %d \n',k+1);
    
        set_Temp_A = [set_W1.A*[C_aS*A_aS^(k+1), C_aL]; 
                       set_Xlambda.A*[zeros(size(C_aL,1),nx), C_aL]]; 
        set_Temp_b = [set_W1.b; set_Xlambda.b];
        set_Temp = Polyhedron('A', set_Temp_A, 'b', set_Temp_b);
        set_Temp.minHRep();
        set_Xkp1 = set_Xk.intersect(set_Temp);
        set_Xkp1.minHRep();
        pause(0.1);
    
        if set_Xkp1 == set_Xk
            found = true;
        else
            set_Xk = set_Xkp1;
            k = k + 1;
            % plot(set_Xk.projection(1:nx), 'color', 'black', 'alpha', 0.01, 'edgecolor', 'black');
        end
    end
    
    set_Xlambda_a = U*set_Xkp1; % Equation (2.10) - E. G. Gilbert
    set_Xlambda_a.minHRep();
    % pltf = plot(set_Xlambda_a.projection(1:nx), 'color', 'green', 'alpha', 0.01, 'edgecolor', 'green');
    % 
    % legend(pltf,"$\Omega_{t,\bar{K}} = \mathrm{Proj}_x(\mathcal{X}_{\lambda}^{a})$","Interpreter","latex");
    % xlabel("$x_1$","Interpreter","latex");
    % ylabel("$x_2$","Interpreter","latex");
    % hold off;
end