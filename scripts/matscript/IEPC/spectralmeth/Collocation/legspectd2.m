%
% This function returns the second derivative spectral matrix
%

N = 3;

Np=N+1;
T = zeros([Np-2 Np]);
for m = 0:N-2
  p = m+2:2:N;
  T(m+1,p+1) = (m+0.5)* (p.*(p+1) - m*(m+1))
  pause;
end

