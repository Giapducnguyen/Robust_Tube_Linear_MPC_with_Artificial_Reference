%% Cost function
function V_N = getCost(xi_seq, yt, N, Q, R, P, T, B, M_theta)

% Dimensions
[nx, nu] = size(B);

% Extract nominal state sequence from the decision variable
xbar_seq = reshape(xi_seq(1:nx*(N+1)), nx, N+1);

% Extract nominal control sequence from the decision variable
ubar_seq = reshape(xi_seq(nx*(N+1)+1:nx*(N+1)+nu*N), nu, N);

% Extract artificial reference from the decision variable
theta = xi_seq(end-nu+1:end);
z_s = M_theta*theta;
x_s = z_s(1:nx);
u_s = z_s(nx+1:end);

% Initialize cost value
V_N = 0;
% Compute cost over prediction horizon
for k = 1 : N
    % Term 1: state tracking
    term1 = (xbar_seq(:,k)-x_s)'*Q*(xbar_seq(:,k)-x_s);
    % Term 2: control tracking
    term2 = (ubar_seq(:,k)-u_s)'*R*(ubar_seq(:,k)-u_s);
    V_N = V_N + term1 + term2;
end
% Term 3: terminal state tracking
term3 = (xbar_seq(:,N+1)-x_s)'*P*(xbar_seq(:,N+1)-x_s);
% Term 4: artificial steady state deviation from target reference
term4 = (x_s-[yt;0])'*T*(x_s-[yt;0]);
% Total cost
V_N = V_N + term3 + term4;
end
