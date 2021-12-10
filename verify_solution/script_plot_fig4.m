%% preliminary
clf
% addpath('../../toolbox/chebfun-master/')
addpath('../verify_defect/')

N = 25; % # of Fourier coefficients
n = 13; % # of Chebyshev coefficients
% stepsize = 0.08/2^4; % length of time step
% tspan = [0,stepsize];
rigorous = 0;

% angle = 45;
% blowuptime = 0.0119;
% ontop = blowuptime/cos(angle/180*pi);
% timedivide = 2^4;
% stepsize = ontop/timedivide;
% tspan = [0,stepsize];
angle = 60;
upper_bound_of_blowuptime = 0.0145; % a bit improved!
timedivide = 2^4;
stepsize = upper_bound_of_blowuptime/timedivide;
tspan = [0,stepsize];


a0 = zeros(1,2*N+1); % Initial sequence
a0(N+1) = 50; a0(N) = -25; a0(N+2) = -25;

y_local = zeros(1,10);
y = []; % Data container

%% getting approximate solution and residual bounds
% timestep = 1;
for timestep=1:2*timedivide
%   timestep
%   if mod(timestep-1,2)==0 && angle < 0
%   if mod(timestep-1,2)==0 && angle > 0 || timestep == timedivide+1
  if mod(timestep-1,2)==0
%     y = [y,timestep];
    if timestep == timedivide+1 || timestep == 1
      iro = 'r';
    else
      iro = 'k';
    end
    plot_profile(a0,iro)
  end
  if timestep == timedivide + 1
    angle = -angle;
    subplot(1,2,1)
    text(0.68,90,'$O$','interpreter','latex','Color','red')
    text(0.45,32,'$z_A$','interpreter','latex','Color','red')
    subplot(1,2,2)
    text(0.51,5,'$O$','interpreter','latex','Color','red')
    text(0.85,35,'$z_A$','interpreter','latex','Color','red')
    figure
    plot_profile(a0,'r')
  end
  [a, d_N, d_infty] = getting_the_solution_timestepping(N,n,tspan,a0,angle,rigorous);% Output is one-sided Chebyshev!
  disp(['delta_N = ',num2str(d_N)])
  disp(['delta_tail = ',num2str(d_infty)])
  d_all = d_N+d_infty;
  h = stepsize;

  a0 = sum(a)+sum(a(2:end,:));% initial sequence of next step
  tspan = tspan + stepsize;% next time step
end
plot_profile(a0,'r')
subplot(1,2,1)
text(0.85,-30,'$z_A$','interpreter','latex','Color','red')
text(0.7,-420,'$z_C$','interpreter','latex','Color','red')
subplot(1,2,2)
text(0.65,30,'$z_A$','interpreter','latex','Color','red')
text(0.8,500,'$z_C$','interpreter','latex','Color','red')