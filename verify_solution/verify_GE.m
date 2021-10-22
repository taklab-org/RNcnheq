function success = verify_GE(rc,rs,mu,theta,ac)
% rc: absolute value of center direction s.t. |a_0|<rc
% rs: absolute value of stable direction s.t. |a_s|<rs
% mu: cos(theta)

% getting the rho value
D = ((mu-4*rc-2*rs)/(4*rs))^2 - 8*rs;
if D>=0
%   rho = (mu-4*rc-2*rs)/(8*rs)-0.5*sqrt(D);
  rho = 4*rs/((mu-4*rc-2*rs)/(4*rs)+sqrt(D));
  if rho<0 || isinf(rho) || rho>0.02
    success = 0;
    disp('cannot obtain rho')
    return
  end
else
  success = 0;
  disp('cannot obtain rho')
  return
end

% check the global existence
delta1 = 2*rc+(1+2*rho)*rs; if delta1 >= mu, success = 0; return, end
delta2 = 2*rc+2*(1+rho)*rs; if delta2 >= mu, success = 0; return, end
% delta3 = 2*(rho*(rc+rho*rs)+rs);
delta4 = 2*(rc+2*rho*rs+rs); if delta4 >= mu, success = 0; return, end
lambda = 4*(rho*(rc+rho*rs)+rs)*rs/((mu-delta1)*(mu-delta4)) + 2*(rc+rho*rs)/(mu-delta1);
% 
if mu>delta2+2*rc && lambda<1 && rho*rs < min(0.02*rs, -real(exp(1i*theta)*ac)) % it is assumed that rc := |a_c| + 0.02*|a_s| and r_s := |a_s|
  success = 1;
  disp('success for GE')
else
  success = 0;
  disp('fail to obtain GE')
end