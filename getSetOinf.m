%% Function helper: Find terminal region
% % Reference:
% % "Robust tube-based MPC for tracking of constrained linear systems with
% %  additive disturbances" - Subsection 3.2       < D. Limon et al. >

% System:                       x^+ = A x + B u
% Steady state and control:     (xs, us)
% Control law:                  u = K(x - xs) + us
% Error dynamics:               e^+ = (A + B K) e = A_K e

% Step 1: find invariant set for the error dynamics
% Step 2: find corresponding invariant set for the state

function set_Oinf = getSetOinf(iter_max, set_X, set_U, r, A, B, K)
A_K = A + B*K;
[nx,~] = size(B);
xs = [r; 0];    us = B\((eye(nx)-A)*xs); % steady state and control

k = 0; found = false;
set_Es = Polyhedron('A', [set_X.A; set_U.A*K], ...
            'b', [set_X.b-set_X.A*xs; set_U.b-set_U.A*us]);
set_Es.minHRep();
set_Es_k = set_Es;

fprintf("\nComputing the maximal invariant set: \n");
% figure; grid on; hold on; box on;
while (k < iter_max && ~found)
    fprintf("Iteration: %d \n", k+1);

    set_temp = invAffineMap(set_Es_k, A_K);

    set_temp.minHRep();
    set_Es_kp1 = set_temp.intersect(set_Es_k);
    set_Es_kp1.minHRep();
    pause(0.1);
    if set_Es_kp1 == set_Es_k
        found = true;
        set_Oinf = set_Es_k+xs; set_Oinf.minHRep();
        % plot(set_Xf,'color',[0.8,0.8,0.8],'alpha',0.01,'edgecolor','red','Linestyle','-','LineWidth',1);
    else
        set_Es_k = set_Es_kp1;
        % plot(set_Es_k+xs,'color',[0.8,0.8,0.8],'alpha',0.01,'edgecolor',[0.8,0.8,0.8],'Linestyle','--','LineWidth',1);
        k = k + 1;
    end
end
% hold off;
end
