%% Function helper: Find gain K to minimize disturbance invariant set
% % Reference:
% % "Robust tube-based MPC for tracking of constrained linear systems with
% %  additive disturbances" - Subsection 7.2       < D. Limon et al. >

function [set_E, K] = getOptGainK(lambda, rho, A, B, set_X, set_U, set_W)
% Transform constraint sets 
% 1 - State constraint: from A x <= b  to  H x <= 1, i = 1, ..., nrx
nrx = size(set_X.A, 1); nrx = 2*nrx;
H = set_X.A ./ set_X.b; H = [H; -H];
% 2 - State constraint: from A u <= b  to  L u <= 1, i = 1, ..., nru
nru = size(set_X.A, 1); nru = 2*nru;
L = set_U.A ./ set_U.b; L = [L; -L];

% Optimized variables
[nx,nu] = size(B);
W = sdpvar(nx, nx);
Y = sdpvar(nu, nx);
gamma = sdpvar;

% Vertices of disturbance set
V_W = set_W.V;
m = size(V_W,1);

% % Constraint 1
Cons = [];
for k = 1:m
    w = V_W(k,:)';
    Cons = [Cons, ...
                [lambda*W,    zeros(1,nx)', (A*W+B*Y)';
                 zeros(1,nx), 1-lambda,     w';
                 (A*W+B*Y),   w,            W          ] >= 0];
end
clearvars k;
% % Constraint 2
for k = 1:nru
    Cons = [Cons, ...
                  [rho^2,        (Y'*L(k,:)')';
                   (Y'*L(k,:)'), W             ] >= 0];
end
clearvars k;
% % Constraint 3
for k = 1: nrx
    Cons = [Cons, ...
                [gamma,      (W*H(k,:)')';
                (W*H(k,:)'), W            ] >= 0];
end
clearvars k;

% % Options
Opts = sdpsettings('solver', 'sedumi', 'sedumi.eps', 1e-9, 'verbose', 1);

% % Optimization
optimize(Cons, gamma, Opts);

% % Extract for disturbance invariant set and gain K
set_E = inv(value(W));
K = value(Y)*inv(value(W));

end
