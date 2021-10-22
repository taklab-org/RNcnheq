%% preliminary
clear
% addpath('../../toolbox/chebfun-master/')
% addpath('../verify_defect/')
% addpath('../variational_problem/')
addpath('../chebfun-master/')
addpath('../verify_defect/')
addpath('../variational_problem/')
%
N = 25; % # of Fourier coefficients
n = 13; % # of Chebyshev coefficients
%
rigorous = 1;% switch rigosous computation (for test)
%
angle = 60;
% blowuptime = 0.0119;
% ontop = blowuptime/cos(angle/180*pi);
upper_bound_of_blowuptime = 0.0145; % a bit improved!
timedivide = 2^6;
stepsize = upper_bound_of_blowuptime/timedivide;
tspan = [0,stepsize];
core = 0; % core = 0 is best!
%
% Initial sequence
a0 = zeros(1,2*N+1);
a0(N+1) = 50; a0(N) = -25; a0(N+2) = -25;
%
if rigorous>0
  y_local = intval(zeros(1,10));
else
  y_local = zeros(1,10);
end
y = []; % Data container
%
eps_all = 0;
%% getting approximate solution and residual bounds
for timestep=1:2*timedivide
  timestep
  if timestep == timedivide + 1
    angle = -angle;
  end
  [a, d_N, d_infty] = getting_the_solution_timestepping(N,n,tspan,a0,angle,rigorous);% Output is one-sided Chebyshev!
  disp(['delta_N = ',num2str(d_N)])
  disp(['delta_tail = ',num2str(d_infty)])
  h = stepsize;

  if rigorous>0
    if sup(compute_eps(intval(a),intval(a0)))/norm(a0,1) > 1e-14
      error('approxiate solution is rough')
    end
    eps_all = sup(eps_all + compute_eps(intval(a),intval(a0)));
    d_all = sup(d_N+d_infty);
  else
    eps_all = compute_eps(a,a0);
    d_all = d_N+d_infty;
  end
  
  a0 = sum(a)+sum(a(2:end,:));% initial sequence of next step (two-sided)
  
  %   plot_solution(a,tspan,1)
  
  
  %% start to solve variational problem
  [C,r_minus] = solve_variational_equation(a,angle,h,N,n,core,rigorous);
  % r_minus
  if rigorous>0
    C0 = reshape(C,n,2*core+1,2*core+1);
    C0k = reshape(sum(abs(intval(C0(1,:,:))),2),2*core+1,1);
    M_phi = max(2*(sum(abs(intval(C)))'+r_minus)-C0k);
  else
    C0 = reshape(C,n,2*core+1,2*core+1);
    C0k = reshape(sum(abs(C0(1,:,:)),2),2*core+1,1);
    M_phi = max(2*(sum(abs(C))'+r_minus)-C0k);
  end
  
  [C_backward,r_minus_backward] = solve_variational_equation_backward(a,angle,h,N,n,core,rigorous);
  if rigorous>0
    C0_backward = reshape(C_backward,n,2*core+1,2*core+1);
    C0k_backward = reshape(sum(abs(intval(C0_backward(1,:,:))),2),2*core+1,1);
    M_psi = max(2*(sum(abs(intval(C_backward)))'+r_minus_backward)-C0k_backward);
  else
    C0_backward = reshape(C_backward,n,2*core+1,2*core+1);
    C0k_backward = reshape(sum(abs(C0_backward(1,:,:)),2),2*core+1,1);
    M_psi = max(2*(sum(abs(C_backward))'+r_minus_backward)-C0k_backward);
  end
  
  M0 = M_phi*M_psi;
  
  C0(2:end,:,:) = 2*C0(2:end,:,:);
  if rigorous>0
    phi_at_endpoint = norm(reshape(sum(intval(C0),1),2*core+1,2*core+1),1);
  else
    phi_at_endpoint = norm(reshape(sum(C0,1),2*core+1,2*core+1),1);
  end
  
  %% start to estimate the evolution operator
  if rigorous>0
    ipi = intval('pi');
    alp_N = (core+1)^2*(2*ipi)^2*cos(angle*ipi/180);
    a_X = a_norm(intval(a));
    a_infty = a; a_infty(:,N+1) = 0;
    a_inf_X = a_norm(intval(a_infty));
  else
    alp_N = (core+1)^2*(2*pi)^2*cos(angle*pi/180);
    a_X = a_norm(a);
    a_infty = a; a_infty(:,N+1) = 0;
    a_inf_X = a_norm(a_infty);
  end
  
  if alp_N>2*a_X
    M_infinity = (1-exp(-(alp_N-2*a_X)*h))/(alp_N-2*a_X);
    M_bar_infty = (h-M_infinity)/(alp_N-2*a_X);
    kappa = 1 - 4*M0*a_inf_X^2*M_bar_infty;
  else
    M_infinity = (exp((2*a_X-alp_N)*h)-1)/(2*a_X-alp_N);
    M_bar_infty = (M_infinity-h)/(2*a_X-alp_N);
    kappa = 1 - 4*M0*a_inf_X^2*M_bar_infty;
  end
  M_infinite_sup = max(1,exp((2*a_X-alp_N)*h));
  M_infinite_at_endpoint = exp(-(alp_N-2*a_X)*h);
  
  if kappa>0
    U_matrix = [M0/kappa, 2*M0*M_infinity*a_inf_X/kappa;...
      2*M0*M_infinity*a_inf_X/kappa, M_infinite_sup+4*M0*M_infinity^2*a_inf_X^2/kappa];
    M = sup(norm(U_matrix,1));
    U_matrix_at_endpoint1 = [phi_at_endpoint*(1+2*M_psi*a_inf_X*h*U_matrix(2,1)),...
      2*phi_at_endpoint*M_psi*a_inf_X*h*U_matrix(2,2);...
      2*M_infinity*a_inf_X*U_matrix(1,1),...
      exp((2*a_X-alp_N)*h)+2*M_infinity*a_inf_X*U_matrix(1,2)];
    U_matrix_at_endpoint2 = [phi_at_endpoint*M_psi*(1+2*a_inf_X*h*U_matrix(2,1)),...
      2*phi_at_endpoint*M_psi*a_inf_X*h*U_matrix(2,2);...
      2*M_infinity*a_inf_X*U_matrix(1,1),...
      M_infinite_sup+2*M_infinity*a_inf_X*U_matrix(1,2)];
    M_at_endpoint = min(M,sup(norm(U_matrix_at_endpoint1,1)));
    Ms = min(M,sup(norm(U_matrix_at_endpoint2,1)));
  else
    error('Linearized problem is not solved (kappa<0)')
  end
  
  
  %% verify the contraction mapping: require INTLAB for using "verifynlssall"
  F = @(x) M*(eps_all+stepsize*(2*x.^2+d_all))-x;
  
  [xx,xx_cand,data] = verifynlssall(F,infsup(0,2));
  if ~any(xx>0)
    error('contraction mapping is not verified!')
  end
  while(1)
    if isempty(xx_cand)
      err = sup(xx(:,all(F(sup(xx))<0,1)));
      break
    else
      [xx,xx_cand,data] = verifynlssall(data);
    end
  end
  
  err_at_endpoint = intval(M_at_endpoint)*eps_all+intval(Ms)*stepsize*(2*intval(err)^2+intval(d_all));
  err_at_endpoint = min(err,sup(err_at_endpoint));
  
  %% Data
  y_local(1) = tspan(1);
  y_local(2) = tspan(2);
  y_local(3) = M0;
  y_local(4) = M_infinite_sup;
  y_local(5) = M;
  y_local(6) = a_X;
  y_local(7) = kappa;
  y_local(8) = eps_all;
  y_local(9) = d_N+d_infty;
  y_local(10) = err;
  y = [y;y_local];
  
  %% Update the initial error and time interval
  eps_all = err_at_endpoint;
  tspan = tspan + stepsize;% next time step
  
end

%% Proof of non-negativity of u(t,x)
a0 = intval(a0);
a0_flip = fliplr(a0);
if ~in(0, 0.5*norm(imag(a0) + imag(a0_flip) - 1i*(real(a0)-real(a0_flip)),1) - eps_all)
    disp('Existence of branching singularity is proved.')
end
% printresult_timestepping
save data_upper_bound.mat