LW = 'linewidth';

% fig5a-1
load data_GE_60.mat
figure
plot(mid([y(:,1),y(:,2)])',mid([y(:,6),y(:,6)])','Color',[0 0.4470 0.7410],LW,1.6)
% stairs(mid([y(:,1)]),mid(y(:,6)),LW,1.6)
xlabel('$t$','interpreter','latex')
ylabel('$\|\bar{a}\|_X$','interpreter', 'latex')
axis([0,0.21,0,1.05*max(mid(y(:,6)))])
title('The values of $\bar{a}$ in $X$-norm ($\theta=\pi/3$)','interpreter', 'latex')
% % fig5a-2
figure
semilogy(mid([y(:,1);y(end,2)]),mid([y(:,8);err_at_endpoint]),'.-',LW,1.6)
xlabel('$t$','interpreter','latex')
ylabel('$\varepsilon_i$','interpreter', 'latex')
yticks([1e-14,1e-12,1e-10,1e-8])
% xticklabels({
%   '0','10^{-16}','10^{-14}','10^{-12}','10^{-10}',...
%   '10^{-8}','10^{-6}','10^{-4}','10^{-2}','1'
% })
axis([0,0.21,1e-14,2*max(mid([y(:,8);err_at_endpoint]))])
title('The point-wise error estimate ($\theta=\pi/3$)','interpreter', 'latex')

% fig5b-1
load data_GE_45.mat
figure
plot(mid([y(:,1),y(:,2)])',mid([y(:,6),y(:,6)])','Color',[0 0.4470 0.7410],LW,1.6)
% stairs(mid([y(:,1)]),mid(y(:,6)),LW,1.6)
xlabel('$t$','interpreter','latex')
ylabel('$\|\bar{a}\|_X$','interpreter', 'latex')
axis([0,0.155,0,1.05*max(mid(y(:,6)))])
title('The values of $\bar{a}$ in $X$-norm ($\theta=\pi/4$)','interpreter', 'latex')
% fig5b-2
figure
semilogy(mid([y(:,1);y(end,2)]),mid([y(:,8);err_at_endpoint]),'.-',LW,1.6)
xlabel('$t$','interpreter','latex')
ylabel('$\varepsilon_i$','interpreter', 'latex')
yticks([1e-14,1e-12,1e-10,1e-8,1e-6])
axis([0,0.155,1e-14,2*max(mid([y(:,8);err_at_endpoint]))])
title('The point-wise error estimate ($\theta=\pi/4$)','interpreter', 'latex')

% fig5c-1
load data_GE_30.mat
figure
plot(mid([y(:,1),y(:,2)])',mid([y(:,6),y(:,6)])','Color',[0 0.4470 0.7410],LW,1.6)
xlabel('$t$','interpreter','latex')
ylabel('$\|\bar{a}\|_X$','interpreter', 'latex')
axis([0,0.132,0,1.05*max(mid(y(:,6)))])
title('The values of $\bar{a}$ in $X$-norm ($\theta=\pi/6$)','interpreter', 'latex')
% fig5c-2
figure
semilogy(mid([y(:,1);y(end,2)]),mid([y(:,8);err_at_endpoint]),'.-',LW,1.6)
xlabel('$t$','interpreter','latex')
ylabel('$\varepsilon_i$','interpreter', 'latex')
yticks([1e-14,1e-12,1e-10,1e-8,1e-6,1e-4])
axis([0,0.132,1e-14,2*max(mid([y(:,8);err_at_endpoint]))])
title('The point-wise error estimate ($\theta=\pi/6$)','interpreter', 'latex')

% fig5d-1
load data_GE_15.mat
figure
plot(mid([y(:,1),y(:,2)])',mid([y(:,6),y(:,6)])','Color',[0 0.4470 0.7410],LW,1.6)
xlabel('$t$','interpreter','latex')
ylabel('$\|\bar{a}\|_X$','interpreter', 'latex')
axis([0,0.125,0,1.05*max(mid(y(:,6)))])
title('The values of $\bar{a}$ in $X$-norm ($\theta=\pi/12$)','interpreter', 'latex')
% fig5d-2
figure
semilogy(mid([y(:,1);y(end,2)]),mid([y(:,8);err_at_endpoint]),'.-',LW,1.6)
xlabel('$t$','interpreter','latex')
ylabel('$\varepsilon_i$','interpreter', 'latex')
yticks([1e-12,1e-8,1e-4])
axis([0,0.125,1e-14,2*max(mid([y(:,8);err_at_endpoint]))])
title('The point-wise error estimate ($\theta=\pi/12$)','interpreter', 'latex')
