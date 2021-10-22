function plot_profile(a,color)
% if isintval(a)
%   a = mid(a);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Input is the one-sided in Chebyshev %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chebfunpref.setDefaults('factory');
% index = 1;
[~,M] = size(a); %M=2*N+1
N_pad = 150;
M = M + 2*N_pad;

N = (M-1)/2;
dx = 1/(2*N-1);
x = dx*(0:2*N);
a0 = [zeros(1,N_pad),a,zeros(1,N_pad)];
% a0 = [a(1,:);2*a(2:end,:)];

% Plot profile:
% figure
subplot(1,2,1);
plot(x,(2*N+1)*real(ifft(ifftshift(chebcoeffs2chebvals(a0),2).')).',color,'linewidth',1.6)
hold on
xlabel('$x$','interpreter','latex'), ylabel('$\mathrm{Re}(\bar{u})$','interpreter', 'latex')%, zlabel('$\mathrm{Re}(u)$','interpreter', 'latex')
% title('Real part')
% figure
subplot(1,2,2);
plot(x,(2*N+1)*imag(ifft(ifftshift(chebcoeffs2chebvals(a0),2).')).',color,'linewidth',1.6)
hold on
xlabel('$x$','interpreter','latex'), ylabel('$\mathrm{Im}(\bar{u})$','interpreter', 'latex')%, zlabel('$\mathrm{Im}(u)$','interpreter', 'latex')
% title('Imaginary part')