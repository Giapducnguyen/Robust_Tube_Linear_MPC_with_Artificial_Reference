%% Constraint function
function [c_ineq, c_eq] = getConstraints(xi_seq, x0, N, A, B, ...
        A_Xbar, b_Xbar, A_Ubar, b_Ubar, A_Xlambda_a, b_Xlambda_a, A_PhiK, b_PhiK)
% Dimensions
[nx, nu] = size(B);
% Extract nominal state sequence from the decision variable
xbar_seq = reshape(xi_seq(1:nx*(N+1)), nx, N+1);
% Extract nominal control sequence from the decision variable
ubar_seq = reshape(xi_seq(nx*(N+1)+1:nx*(N+1)+nu*N), nu, N);
% Extract artificial reference from the decision variable
theta = xi_seq(end-nu+1:end);

% State continuity constraints
c_eq_sc = nan(nx*N,1); 
for k = 1 : N
    xkp1 = A*xbar_seq(:,k) + B*ubar_seq(:,k);
    c_eq_sc(nx*(k-1)+1:nx*(k-1)+nx, 1) = xbar_seq(:,k+1)-xkp1;
end

% Initial constraint:
c_ineq_x0 = -A_PhiK*xbar_seq(:,1) - b_PhiK + A_PhiK*x0;

% State hard constraints from 0 to N-1
size1 = size(b_Xbar, 1);
c_ineq_x = nan(size1*N,1);
for k = 1 : N
    c_ineq_x(size1*(k-1)+1 : size1*(k-1)+size1) = A_Xbar*xbar_seq(:,k)-b_Xbar;
end

% Control hard constraints from 0 to N-1
size2 = size(b_Ubar, 1);
c_ineq_u = nan(size2*N, 1);
for k = 1 : N
    c_ineq_u(size2*(k-1)+1 : size2*(k-1)+size2) = A_Ubar*ubar_seq(:,k)-b_Ubar;
end

% Terminal constraint at N: 
% Terminal constraint
c_ineq_xf = A_Xlambda_a * [xbar_seq(:,N+1); theta] - b_Xlambda_a;

% % Concatenate constraints
% Inequality constraints
c_ineq = [c_ineq_x0; % initial state constraint
          c_ineq_x; % hard constraint on state(from 0 to N-1)
          c_ineq_u; % hard constraint on control(from 0 to N-1)
          c_ineq_xf]; % Terminal constraint on state and artificial reference
% Equality constraint
c_eq = c_eq_sc; % Continuity constraints (Multiple shooting method)

end