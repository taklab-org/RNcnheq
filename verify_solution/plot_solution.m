function plot_solution(a,tspan,index)
% if isintval(a)
%   a = mid(a);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Input is the one-sided in Chebyshev %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chebfunpref.setDefaults('factory');
% index = 1;
[n,M] = size(a); %M=2*N+1
N = (M-1)/2;
dx = 1/(2*N-1);
x = dx*(0:2*N);
a0 = [a(1,:);2*a(2:end,:)];
if length(tspan)<2
  tspan=chebpts(n,[0,tspan]);
else
  tspan=chebpts(n,tspan);
end

if index==1
% Plot profile:
% figure
subplot(1,2,1);
mesh(x,tspan,(2*N+1)*real(ifft(ifftshift(chebcoeffs2chebvals(a0),2).')).')
hold on
xlabel('$x$','interpreter','latex'), ylabel('$t$','interpreter', 'latex'), zlabel('$\mathrm{Re}(u)$','interpreter', 'latex')
title('Real part')
% figure
subplot(1,2,2);
mesh(x,tspan,(2*N+1)*imag(ifft(ifftshift(chebcoeffs2chebvals(a0),2).')).')
hold on
xlabel('$x$','interpreter','latex'), ylabel('$t$','interpreter', 'latex'), zlabel('$\mathrm{Im}(u)$','interpreter', 'latex')
title('Imaginary part')

% 
elseif index==2
% figure
% Plot fourier modes:
k = (-N:N)';
subplot(1,2,1);
mesh(k,tspan,abs(real(chebcoeffs2chebvals(a0))))
% hold on
xlabel('$k$','interpreter','latex'), ylabel('$t$','interpreter', 'latex'), zlabel('$\mathrm{Re}(u)$','interpreter', 'latex')
title('Real part')
set(gca, 'ZScale', 'log')
% figure
subplot(1,2,2);
mesh(k,tspan,abs(imag(chebcoeffs2chebvals(a0))))
% hold on
xlabel('$k$','interpreter','latex'), ylabel('$t$','interpreter', 'latex'), zlabel('$\mathrm{Im}(u)$','interpreter', 'latex')
title('Imaginary part')
set(gca, 'ZScale', 'log')
elseif index==3
% norm plot
figure
LW = 'linewidth'; lw = 1.6;
y=ifft(ifftshift(chebcoeffs2chebvals(a0),2).');
plot(tspan,abs(max(y)),LW,lw);
xlabel('$t$','interpreter','latex'), ylabel('$\|u(t)\|_{\infty}$','interpreter', 'latex')
title('The maximum norm of the solution.')
elseif index==4
% norm plot
figure
LW = 'linewidth'; lw = 1.6;
y=a0;
semilogy(abs(y),LW,lw);
xlabel('$\ell$','interpreter','latex'), ylabel('$|\tilde{a}_{\ell,k}|$','interpreter', 'latex')
title('Chebyshev coefficients.')
end