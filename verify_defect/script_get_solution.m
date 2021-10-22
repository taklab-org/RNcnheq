N = 50; % # of Fourier coefficients
n = 141; % # of Chebyshev coefficients
tmax = 0.08;
angle = 45;

rigorous = 0;

[a, d_N, d_infty] = getting_the_solution(N,n,tmax,angle,rigorous);
disp(d_N)
disp(d_infty)
h = tmax;
% save approximate_solution_v2 a angle h N n
% save(['approximate_solution_',datestr(now,'yyyymmdd_HHMMSS'),'.mat'],'a','angle', 'h', 'N', 'n')