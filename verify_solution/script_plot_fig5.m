load data_upper_bound.mat

figure

LW = 'linewidth'; MS = 'markersize';
semilogy([0:128],mid([y(:,8);eps_all]),'.-',LW,1.6,MS,12)

xlabel('$i$','interpreter','latex')
ylabel('$\varepsilon_i$','interpreter', 'latex')
% axis([0,33,1e-15,2e-4])

figure

semilogy(mid(y(:,10)),'.-',LW,1.6,MS,12)
xlabel('$i$','interpreter','latex')
ylabel('$\rho_i$','interpreter', 'latex')
% axis([0,33,5e-11,2e-3])