% % Simulation illustrating the paper:
% % "Robust tube-based MPC for tracking of constrained linear systems with
% %  additive disturbances" - D. Limon et al.

close all
clear all
clc

%% System info - Subsection 3.2
A = [1, 1; 0, 1];
B = [0, 0.5; 1, 0.5];
C = [1, 0];
D = [0, 0];

x_max = 5;
u_max = 0.3;
w_max = 0.1;

[nx,nu] = size(B);
ny = size(C,1);

%% State feedback computation - Subsection 3.2
Q = eye(nx);
R = 10*eye(nu);
[Kbar,P,~] = dlqr(A,B,Q,R);
Kbar = -Kbar;
T = 100.*P;

%% Minimal robust invariant set computation - Subsection 3.2
% % Error dynamics: e^{+} = A_K e + w;
set_X = Polyhedron('lb', -x_max*ones(nx,1), 'ub', x_max*ones(nx,1));
set_U = Polyhedron('lb', -u_max*ones(nu,1), 'ub', u_max*ones(nu,1));
set_W = Polyhedron('lb', -w_max*ones(nx,1), 'ub', w_max*ones(nx,1));
A_Kbar = A + B*Kbar;

epsilon = 1e-4;
set_PhiKbar = eRPIApprox(epsilon, A_Kbar, set_W);

% % Tightened sets with gain Kbar computed based on LQR
set_Xbar_LQR = set_X - set_PhiKbar;        set_Xbar_LQR.minHRep();
set_Ubar_LQR = set_U - (Kbar*set_PhiKbar); set_Ubar_LQR.minHRep();

%% Terminal region computed as the maximal invariant set - Subsection 3.2
r1 = -3;  % ref 1
r2 = 4.2; % ref 2

iter_max = 30;

set_Oinf_r1 = getSetOinf(iter_max, set_Xbar_LQR, set_Ubar_LQR, r1, A, B, Kbar);
set_Oinf_r2 = getSetOinf(iter_max, set_Xbar_LQR, set_Ubar_LQR, r2, A, B, Kbar);

% Get domain of attraction for N = 2
N = 2;
set_Xbar_r1 = getDomainOfAttraction(N, A, B, set_Oinf_r1, set_Xbar_LQR, set_Ubar_LQR);
set_Xbar_r1 =  set_Xbar_r1+set_PhiKbar; set_Xbar_r1.minHRep();
set_Xbar_r2 = getDomainOfAttraction(N, A, B, set_Oinf_r2, set_Xbar_LQR, set_Ubar_LQR);
set_Xbar_r2 =  set_Xbar_r2+set_PhiKbar; set_Xbar_r2.minHRep();

% Figure 2 - page 251
figure; grid on; hold on; box on;
plot(r1,0,'Marker','*','MarkerSize',6, 'Color','black');
plot(set_Oinf_r1,'color',[0.8,0.8,0.8],'alpha',0.01,'edgecolor','blue','Linestyle','--','LineWidth',1);
plot(set_Xbar_r1,'color',[0.8,0.8,0.8],'alpha',0.01,'edgecolor','black','Linestyle','--','LineWidth',1);
plot(r2,0,'Marker','o','MarkerSize',6, 'Color','black');
plot(set_Oinf_r2,'color',[0.8,0.8,0.8],'alpha',0.01,'edgecolor','blue','Linestyle','-','LineWidth',1);
plot(set_Xbar_r2,'color',[0.8,0.8,0.8],'alpha',0.01,'edgecolor','black','Linestyle','-','LineWidth',1);
legend({"$r_1$", "$O_{\infty}(r_1)$","$X_{2}(r_1)$",...
    "$r_2$", "$O_{\infty}(r_2)$","$X_{2}(r_2)$"}, 'Interpreter','latex','Fontsize',12);
xlabel("$x_1$","Interpreter","latex","FontSize",12);
ylabel("$x_2$","Interpreter","latex","FontSize",12);
axis([-5 5 -2.5 2.5]);
hold off;

%% Minimal robust invariant set computation - Subsection 7.2
lambdaK = 0.65; rho = 0.48;
[~, K] = getOptGainK(lambdaK, rho, A, B, set_X, set_U, set_W);
A_K = A + B*K;
set_PhiK = eRPIApprox(epsilon, A_K, set_W);

figure; grid on; hold on; box on;
plot(set_PhiKbar,'color',[0.8,0.8,0.8],'alpha',0.01,'edgecolor','blue','Linestyle','--','LineWidth',1);
plot(set_PhiK,'color',[0.8,0.8,0.8],'alpha',0.01,'edgecolor','red','Linestyle','--','LineWidth',1); 
legend({"$\phi_{\bar{K}}$ (based on LQR)","$\phi_{K}$ (optimized with LMI)"}, 'Interpreter','latex','Fontsize',12);
title("Disturbance invariant set comparison", 'Interpreter','latex','Fontsize',12);
xlabel("Error of $x_1$","Interpreter","latex","FontSize",12);
ylabel("Error of $x_2$","Interpreter","latex","FontSize",12);
hold off;

% % Tightened sets with gain K optimized using LMI
set_Xbar_opt = set_X - set_PhiK;             set_Xbar_opt.minHRep();
set_Ubar_opt = set_U - (K*set_PhiK);      set_Ubar_opt.minHRep();

figure; hold on; grid on; box on;
plot(set_U, 'color', [0.8,0.8,0.8], 'alpha', 0.01, 'edgecolor', 'black', 'linestyle', '-.', 'linewidth', 1);
plot(set_Ubar_LQR, 'color', [0.8,0.8,0.8], 'alpha', 0.01, 'edgecolor', 'blue', 'linestyle', '-', 'linewidth', 1);
plot(Kbar*set_PhiKbar, 'color', [0.8,0.8,0.8], 'alpha', 0.01, 'edgecolor', 'blue', 'linestyle', '--', 'linewidth', 1);
plot(set_Ubar_opt, 'color', [0.8,0.8,0.8], 'alpha', 0.01, 'edgecolor', 'red', 'linestyle', '-', 'linewidth', 1);
plot(K*set_PhiK, 'color', [0.8,0.8,0.8], 'alpha', 0.01, 'edgecolor', 'red', 'linestyle', '--', 'linewidth', 1);
legend({"$U$", "$\bar{U}_{\mathrm{LQR}}$","$\bar{K} \phi_{\bar{K}}$", ...
    "$\bar{U}_{\mathrm{op}}$","$K \phi_{K}$"}, 'Interpreter','latex','Fontsize',12);
title("Tightened control set comparison", 'Interpreter','latex','Fontsize',12);
xlabel("$u_1$","Interpreter","latex","FontSize",12);
ylabel("$u_2$","Interpreter","latex","FontSize",12);
hold off;

figure; hold on; grid on; box on;
plot(set_X, 'color', [0.8,0.8,0.8], 'alpha', 0.01, 'edgecolor', 'black', 'linestyle', '-.', 'linewidth', 1);
plot(set_Xbar_LQR, 'color', [0.8,0.8,0.8], 'alpha', 0.01, 'edgecolor', 'blue', 'linestyle', '-', 'linewidth', 1);
plot(set_PhiKbar, 'color', [0.8,0.8,0.8], 'alpha', 0.01, 'edgecolor', 'blue', 'linestyle', '--', 'linewidth', 1);
plot(set_Xbar_opt, 'color', [0.8,0.8,0.8], 'alpha', 0.01, 'edgecolor', 'red', 'linestyle', '-', 'linewidth', 1);
plot(set_PhiK, 'color', [0.8,0.8,0.8], 'alpha', 0.01, 'edgecolor', 'red', 'linestyle', '--', 'linewidth', 1);
legend({"$X$", "$\bar{X}_{\mathrm{LQR}}$","$\phi_{\bar{K}}$", ...
    "$\bar{X}_{\mathrm{op}}$","$\phi_{K}$"}, 'Interpreter','latex','Fontsize',12);
title("Tightened state set comparison", 'Interpreter','latex','Fontsize',12);
xlabel("$x_1$","Interpreter","latex","FontSize",12);
ylabel("$x_2$","Interpreter","latex","FontSize",12);
hold off;

%% Characterization of the steady states and inputs - Subsection 3.3
A_s = [A-eye(nx),   B,  zeros(nx,ny);
       C,           D,      -eye(ny)];
null_A = null(A_s,"rational");
M_theta = null_A(1:(nx+nu),:);
N_theta = null_A(nx+nu+1:end,:);

%% Calculation of invariant set for tracking - Subsection 3.4
lambda = 0.99;
% % Select 1 of the following computation method:
    % % Method 1 - Direct computation with tightened constraints
    % % Method 2 - Indirect computation with transformation and additional constraints
method = 2;
set_Xlambda_a = getInvariantSetforTracking(method, iter_max, lambda, A, B, Kbar, M_theta, set_Xbar_opt, set_Ubar_opt);
set_Omega_tKbar = set_Xlambda_a.projection(1:nx);

% Get domain of attraction for N = 2
N = 2;
set_XbarN = getDomainOfAttraction(N, A, B, set_Omega_tKbar, set_Xbar_opt, set_Ubar_opt);
set_Xbar_2 = set_XbarN+set_PhiK; set_Xbar_2.minHRep();
set_Xbar_r2 = getDomainOfAttraction(N, A, B, set_Oinf_r2, set_Xbar_LQR, set_Ubar_LQR);
set_Xbar_r2 =  set_Xbar_r2+set_PhiKbar; set_Xbar_r2.minHRep();

% Figure 3 - page 253
figure; grid on; hold on; box on;
plot(r2,0,'Marker','o','MarkerSize',6, 'Color','black');
plot(set_Oinf_r2,'color',[0.8,0.8,0.8],'alpha',0.01,'edgecolor','blue','Linestyle','-','LineWidth',1);
plot(set_Xbar_r2,'color',[0.8,0.8,0.8],'alpha',0.01,'edgecolor','blue','Linestyle','--','LineWidth',1);
plot(set_Omega_tKbar,'color',[0.8,0.8,0.8],'alpha',0.01,'edgecolor','red','Linestyle','-','LineWidth',1);
plot(set_Xbar_2,'color',[0.8,0.8,0.8],'alpha',0.01,'edgecolor','red','Linestyle','--','LineWidth',1);
legend({"$r_2$", "$O_{\infty}(r_2)$","$\bar{X}_{2}(r_2)$","$\Omega_{t,\bar{K}}$","$\bar{X}_2$"}, 'Interpreter','latex','Fontsize',12);
axis([-5 5 -2.5 2.5]);
hold off;

%% Closed-loop simulation settings
N = 3; % prediction horizon

% % Computed constraints
A_Xbar = set_Xbar_opt.A;        b_Xbar = set_Xbar_opt.b;
A_Ubar = set_Ubar_opt.A;        b_Ubar = set_Ubar_opt.b;
A_Xlambda_a = set_Xlambda_a.A;  b_Xlambda_a = set_Xlambda_a.b;
A_PhiK = set_PhiK.A;            b_PhiK = set_PhiK.b;

% Simulation parameters
Tsim = 150; % [s]
Ts = 1;
Nsim = Tsim/Ts; % Maximum simulation iterations

% Initialization of decision variable = [state; control; theta]
x0 = [4.29; -2]; % initial state

%% Case 1: Disturbance observer - OFF

% Initialization: decision variables
xbar_seq = zeros(nx, N+1); % state sequence
ubar_seq = zeros(nu, N); % control sequence
theta = zeros(nu, 1); % artificial parameter

% Initialization: actual state, control
xSim_seq = nan(nx, Nsim); xSim_seq(:,1) = x0;
uSim_seq = nan(nu, Nsim);

yref_seq = nan(ny, Nsim); % desired reference
yaRef_seq = nan(ny, Nsim); % artificial reference

% Initialization: nominal state and control
xbar0_seq = nan(nx, Nsim); % optimal x0
ubarSim_seq = nan(nu, Nsim); % optimal u0

% Initialization: actual disturbance
wSim_seq = nan(nx, Nsim);

% % fmincon options
options = optimoptions("fmincon","Algorithm","interior-point","MaxIterations",100,"Display","iter-detailed");

% Main simulation loop
for k = 1 : Nsim
    % Update current states:
    x0_1 = xSim_seq(:,k);

    % Update reference
    if k*Ts <= 50
        yt = 4.29;
    elseif k*Ts <= 100
        yt = -4.29;
    else
        yt = 0;
    end
    yref_seq(:,k) = yt;
    
    % Update disturbance
    cycle_idx = mod(floor(k*Ts / 12.5), 4);

    switch cycle_idx
        case 0
            w1 = w_max; w2 = w_max;
        case 1
            w1 = w_max; w2 = -w_max;
        case 2
            w1 = -w_max; w2 = -w_max;
        case 3
            w1 = -w_max; w2 = w_max;
    end
    wSim_seq(:,k) = [w1; w2];

    % Initialize decision variable for warm-start
    xi_seq0 = [reshape(xbar_seq, nx*(N+1), 1); 
              reshape(ubar_seq, nu*N, 1);
              theta];

    % Optimization
    [xi_seq_opt, fval, exitflag, output] = fmincon(...
        @(xi_seq) getCost(xi_seq, yt, N, Q, R, P, T, B, M_theta),... % fun
        xi_seq0,...                                                      % x0
        [], [], [], [], [], [], ...                                     % A, b, Aeq, beq, lb, ub
        @(xi_seq) getConstraints(xi_seq, x0_1, N, A, B, ...
        A_Xbar, b_Xbar, A_Ubar, b_Ubar, A_Xlambda_a, b_Xlambda_a, A_PhiK, b_PhiK),... % nonlcon
        options);                                                       % options


    % Regenerate state sequence from the optimal decision variable
    xbar_seqtemp = reshape(xi_seq_opt(1:nx*(N+1)), nx, N+1);
    
    % Regenerate control sequence from the optimal decision variable
    ubar_seq = reshape(xi_seq_opt(nx*(N+1)+1:nx*(N+1)+nu*N), nu, N);
    
    % Regenerate artificial reference from the optimal decision variable
    theta = xi_seq_opt(end-nu+1:end);
    z_s = M_theta*theta;
    x_s = z_s(1:nx);
    u_s = z_s(nx+1:end);
    y_s = N_theta*theta;
    
    % Store optimal nominal control
    ubarSim_seq(:,k) = ubar_seq(:,1);
    % Store actual control
    uSim_seq(:,k) = K*(xSim_seq(:,k)-xbar_seqtemp(:,1)) + ubar_seq(:,1);
    % Store optimal initial state
    xbar0_seq(:,k) = xbar_seqtemp(:,1);
    % Store artificial reference
    yaRef_seq(:,k) = y_s;

    % Shift solution 1 step for next iteration (warm-start)
    xbar_seq = [xbar_seqtemp(:,2:N+1), xbar_seqtemp(:,N+1)];
    ubar_seq = [ubar_seq(:,2:N), K*(xbar_seqtemp(:,N+1)-x_s) + u_s];

    % Propagate system
    xSim_seq(:,k+1) = A*xSim_seq(:,k) + B*uSim_seq(:,k) + wSim_seq(:,k);

end

%% Case 2: Disturbance observer - ON

% Initialization: decision variables
xbar_seq2 = zeros(nx, N+1); % state sequence
ubar_seq2 = zeros(nu, N); % control sequence
theta2 = zeros(nu, 1); % artificial parameter

% Initialization: actual state, control
xSim_seq2 = nan(nx, Nsim); xSim_seq2(:,1) = x0;
uSim_seq2 = nan(nu, Nsim);

yref_seq2 = nan(ny, Nsim); % desired reference
yaRef_seq2 = nan(ny, Nsim); % artificial reference

% Initialization: nominal state and control
xbar0_seq2 = nan(nx, Nsim); % optimal x0
ubarSim_seq2 = nan(nu, Nsim); % optimal u0

% Initialization: actual disturbance
wSim_seq2 = nan(nx, Nsim);

% Initialization: estimated disturbance - Section 6
what_seq = nan(nx, Nsim);
b = 0.95;
H = (C+D*K)*(eye(nx)-A_K)^-1;
yhatt_seq = nan(ny, Nsim); % modified desired output

% % fmincon options
options = optimoptions("fmincon","Algorithm","interior-point","MaxIterations",100,"Display","iter-detailed");

% Main simulation loop
for k = 1 : Nsim
    % Update current states:
    x0_2 = xSim_seq2(:,k);

    % Update reference
    if k*Ts <= 50
        yt = 4.29;
    elseif k*Ts <= 100
        yt = -4.29;
    else
        yt = 0;
    end
    yref_seq2(:,k) = yt;
    
    % Update disturbance
    cycle_idx = mod(floor(k*Ts / 12.5), 4);

    switch cycle_idx
        case 0
            w1 = w_max; w2 = w_max;
        case 1
            w1 = w_max; w2 = -w_max;
        case 2
            w1 = -w_max; w2 = -w_max;
        case 3
            w1 = -w_max; w2 = w_max;
    end
    wSim_seq2(:,k) = [w1; w2];

    % Initialize decision variable for warm-start
    xi_seq0 = [reshape(xbar_seq2, nx*(N+1), 1); 
               reshape(ubar_seq2, nu*N, 1);
               theta2];
    
    % Estimate disturbance
    if k == 1
        what_seq(:,k) = zeros(nx,1);
    else
        what_seq(:, k) = b*(xSim_seq2(:,k) - (A*xSim_seq2(:,k-1) + B*uSim_seq2(:,k-1))) + (1-b)*what_seq(:,k-1);
    end
    
    % Modify desired output
    yhat_t = yt - H*what_seq(:,k);
    yhatt_seq(:, k) = yhat_t;

    % Optimization
    [xi_seq_opt2, fval2, exitflag2, output2] = fmincon(...
        @(xi_seq) getCost(xi_seq, yhat_t, N, Q, R, P, T, B, M_theta),... % fun
        xi_seq0,...                                                      % x0
        [], [], [], [], [], [], ...                                     % A, b, Aeq, beq, lb, ub
        @(xi_seq) getConstraints(xi_seq, x0_2, N, A, B, ...
        A_Xbar, b_Xbar, A_Ubar, b_Ubar, A_Xlambda_a, b_Xlambda_a, A_PhiK, b_PhiK),... % nonlcon
        options);                                                       % options

    % Regenerate state sequence from the optimal decision variable
    xbar_seqtemp2 = reshape(xi_seq_opt2(1:nx*(N+1)), nx, N+1);
    
    % Regenerate control sequence from the optimal decision variable
    ubar_seq2 = reshape(xi_seq_opt2(nx*(N+1)+1:nx*(N+1)+nu*N), nu, N);
    
    % Regenerate artificial reference from the optimal decision variable
    theta2 = xi_seq_opt2(end-nu+1:end);
    z_s = M_theta*theta2;
    x_s = z_s(1:nx);
    u_s = z_s(nx+1:end);
    y_s = N_theta*theta2;
    
    % Store optimal nominal control
    ubarSim_seq2(:,k) = ubar_seq2(:,1);
    % Store actual control
    uSim_seq2(:,k) = K*(xSim_seq2(:,k)-xbar_seqtemp2(:,1)) + ubar_seq2(:,1);
    % Store optimal initial state
    xbar0_seq2(:,k) = xbar_seqtemp2(:,1);
    % Store artificial reference
    yaRef_seq2(:,k) = y_s;

    % Shift solution 1 step for next iteration (warm-start)
    xbar_seq2 = [xbar_seqtemp2(:,2:N+1), xbar_seqtemp2(:,N+1)];
    ubar_seq2 = [ubar_seq2(:,2:N), K*(xbar_seqtemp2(:,N+1)-x_s) + u_s];

    % Propagate system
    xSim_seq2(:,k+1) = A*xSim_seq2(:,k) + B*uSim_seq2(:,k) + wSim_seq2(:,k);

end

%% Plots - Figure 4 - page 254
figure; tiledlayout(3,1, "TileSpacing","compact", "Padding","compact");
Tseq = (0:Ts:Tsim);
nexttile; hold on; grid on; box on;
plot(Tseq, [wSim_seq2(1,:), nan], 'LineStyle','-','Color','red','LineWidth',1);
plot(Tseq, [wSim_seq2(2,:), nan], 'LineStyle','--','Color','blue','LineWidth',1);
legend({"$w_1$","$w_2$"},"Interpreter","latex","FontSize",12);
yticks([-w_max 0 w_max]);
hold off;

nexttile; hold on; grid on; box on;
plot(Tseq, xSim_seq(1,:), 'LineStyle','-','Color','red','LineWidth',1);
plot(Tseq, [yref_seq, nan], 'LineStyle','--','Color','black','LineWidth',1);
plot(Tseq, [yaRef_seq, nan], 'LineStyle','-.','Color','blue','LineWidth',1);
legend({"$y$","$y_t$","$\bar{y}_s$"},"Interpreter","latex","FontSize",12);
title("Disturbance Observer - OFF","Interpreter","latex","FontSize",12);
yticks([-5 0 5]);
hold off;

nexttile; hold on; grid on; box on;
plot(Tseq, xSim_seq2(1,:), 'LineStyle','-','Color','red','LineWidth',1);
plot(Tseq, [yref_seq2, nan], 'LineStyle','--','Color','black','LineWidth',1);
plot(Tseq, [yhatt_seq, nan], 'LineStyle',':','Color','magenta','LineWidth',2);
plot(Tseq, [yaRef_seq2, nan], 'LineStyle','-.','Color','blue','LineWidth',1);
legend({"$y$","$y_t$","$\hat{y}_t$","$\bar{y}_s$"},"Interpreter","latex","FontSize",12);
xlabel("Time","Interpreter","latex", "FontSize", 12);
title("Disturbance Observer - ON","Interpreter","latex","FontSize",12);
yticks([-5 0 5]);
hold off;