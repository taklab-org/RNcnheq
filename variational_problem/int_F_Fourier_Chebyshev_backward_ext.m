function [F] = int_F_Fourier_Chebyshev_backward_ext(c,ta,b,angle,h,N,n,m)

%%% INPUTS
% c = (c_{ell,k})_{k=-N,...,N , ell = 0,...,n-1} : the Fourier-Chebyshev
%                    space-time coefficients of c (in k: Fourier and in ell: Chebyshev)
% ta = (ta_{ell,k})_{k=-N,...,N , ell = 0,...,n-1} : the Fourier-Chebyshev
%                    space-time coefficients of the approximate solution ta
% b = (b_k)_{k=-N,...,N} : row vector of Fourier coeff. of the initial condition
% angle : the angle parameter
% h : step size in time
% N : determines the # of Fourier coefficients which is 2*N+1
% n : # of Chebyshev coefficients

ipi = intval('pi');
theta = angle/180*ipi;

omega = 2*ipi;

c = reshape(c,n,2*m+1); c = [c;zeros(1,2*m+1)];
ta = reshape(ta,n,2*N+1); ta = [ta;zeros(1,2*N+1)];

c = intval(c); ta = intval(ta); b = intval(b); h = intval(h);

%%% Beginning of the computation of the convolution

c_ext = [flipud(c(2:end,:));c];
ta_ext = [flipud(ta(2:end,:));ta];

[~,ta_c_full] = convtensor(ta_ext,[zeros(2*n+1,N-m),c_ext,zeros(2*n+1,N-m)]);

ta_c = [ta_c_full(2*n+1:end,2*N+1-m:2*N+1+m);zeros(1,2*m+1)];

%%% End of the computation of the convolution

F = intval(zeros(2*n-1,2*m+1));

%%% Initial condition

N = m;

for k = -N:N
    F(1,k+N+1) = c(1,k+N+1) + 2*sum( ((-1).^(1:n-1)').*c(2:n,k+N+1)) - b(k+N+1);
end

%%% The rest of the equations

c = [c;zeros(n-1,2*N+1)];

for k = -N:N
    lambda_k = -(h/2)*exp(1i*theta)*omega^2*k^2;
    ind = (1:2*n-2)';
    weight=ones(2*n-2,1); weight(2*n-2)=0;
    F(ind+1,k+N+1) = 2*ind.*c(ind+1,k+N+1) - (-lambda_k*c(ind,k+N+1)+lambda_k*c(ind+2,k+N+1)+h*exp(1i*theta)*weight.*(ta_c(ind+2,k+N+1)-ta_c(ind,k+N+1)));
end

F = reshape(F,(2*N+1)*(2*n-1),1);

end