% Demo code for complex blowup
% 
%% Complex Ginzburg-Landau equation
% JJIAM 2016 Fig.3
figure
hold on
for i=[15,30,45,60]
%   theat1d_cheb_fft(20,150,0.08,i);
  theat1d_cheb(20,150,0.08,i,3);
end
xlabel('$\rho$','interpreter','latex'), ylabel('$\|u\|_{\infty}$','interpreter', 'latex')
title('Fig. 3')
legend({'$\gamma=15$','$\gamma=30$','$\gamma=45$','$\gamma=60$'},'Interpreter','Latex')

%% Fig.4
% tic
% theat1d_cheb(20,150,0.08,45,1);
% toc
% tic
theat1d_cheb1(20,150,0.08,45,1);
% toc

%% Fourier modes
theat1d_cheb(20,150,0.08,45,2);

%% Nonlinear Schoedinger equation
% Fig.7
nschoe_cheb_test(20,100,0.15,1);

%% Fig.8(a)
figure
nschoe_cheb_test(20,100,1,3);
xlabel('$s$','interpreter','latex'), ylabel('$\|u\|_{\infty}$','interpreter', 'latex')
title('Fig. 8(a)')

%% Fourier modes
nschoe_cheb_test(20,100,0.15,2);

%% Fig.11
nschoe_cheb_test2(20,100,0.25,3);
xlabel('$s$','interpreter','latex'), ylabel('$\|u\|_{\infty}$','interpreter', 'latex')
title('Fig. 11')
