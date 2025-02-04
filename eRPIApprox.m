%% Function Helper 1: epsilon- outer approximation of RPI set
%{
Reference:
Rakovic, Sasa V., et al. "Invariant approximations of the minimal robust 
positively invariant set." IEEE Transactions on automatic control 50.3 (2005): 406-410.
%}
function Fs_alpha = eRPIApprox(epsilon,A_K,set_W)
[nx,~] = size(A_K);
Ms = 1000;
s = 0;
alpha = 1000;

mss = zeros(2*nx,1);

while(alpha > epsilon/(epsilon+Ms))
    s = s+1;
    alpha = max(set_W.support(A_K^s*set_W.A')./set_W.b);
    mss = mss + set_W.support([A_K^s, -A_K^s]);
    Ms = max(mss);
end

% figure; hold on;
Fs = set_W;
for i = 1:s-1
    Fs = Fs + A_K^i*set_W; Fs.minHRep(); Fs.minVRep();
    % % Fs.plot('Color','red','alpha',0.1); 
    pause(0.1);
    % plot(Fs,'color',[0.8,0.8,0.8],'alpha',0.01,'edgecolor',[0.8,0.8,0.8],'Linestyle','--','LineWidth',1);
end
Fs_alpha = ((1-alpha)^-1)*Fs;
Fs_alpha.minHRep();
% % Fs_alpha.plot('Color','yellow','alpha',0.35);
% plt = plot(Fs_alpha,'color',[0.8,0.8,0.8],'alpha',0.01,'edgecolor','blue','Linestyle','--','LineWidth',1);
% legend(plt, "$\phi_{\mathrm{K}}$", 'Interpreter','latex','Fontsize',12);
% title('Disturbance invariant set (minimal RPI set)','interpreter','latex','Fontsize',12);
% xlabel('Error of state $x_1$','interpreter','latex','Fontsize',12);
% ylabel('Error of state $x_2$','interpreter','latex','Fontsize',12);
% hold off;

end