% % Domain of attraction with prediction horizon N


function set_XN = getDomainOfAttraction(N, A, B, set_Xf, set_X, set_U)
set_Xk = set_Xf; % initialization -> calculation backward

% figure; hold on;
% title("Computation Evolution of $X_N$","Interpreter","latex","FontSize",12);
% pltb = plot(set_Xf, 'color', 'red', 'alpha', 0.05, 'edgecolor', 'red','Linestyle','-','LineWidth',1);

fprintf('\n Computing the domain of attraction for N = %d \n', N);
for k = 1 : N
    fprintf('Iteration: %d \n',k);

    set_temp = set_Xk + (-B*set_U); set_temp.minHRep();
    % set_Xk = intersect(set_X, Polyhedron('A', set_temp.A*A, 'b', set_temp.b));
    set_temp2 = set_temp.invAffineMap(A); set_temp2.minHRep();
    set_Xk = intersect(set_X, set_temp2);
    set_Xk.minHRep();
    pause(0.1);
    % plot(set_Xk, 'color', [0.8, 0.8, 0.8], 'alpha', 0.01, 'edgecolor', [0.8, 0.8, 0.8],'Linestyle','--');
end
set_XN = set_Xk;
% pltf1 = plot(set_XN, 'color', [0.8, 0.8, 0.8], 'alpha', 0.01, 'edgecolor', 'black','Linestyle','-');
% legend([pltb, pltf1],["$X_{f}$","$X_N$"],"Interpreter","latex","FontSize",12);
% xlabel("$x_1$","Interpreter","latex","FontSize",12);
% ylabel("$x_2$","Interpreter","latex","FontSize",12);
% hold off;
end