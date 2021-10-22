
clear
h = .01; N = 30; n = 30; angle = 45; ta = rand((2*N+1)*n,1); c = rand((2*N+1)*n,1); b = rand(1,2*N+1);

DF = DF_Fourier_Chebyshev(ta,angle,h,N,n);
DF1 = finite_diff_F_Fourier_Chebyshev(c,ta,b,angle,h,N,n);

norm(DF-DF1,1)
