%
% This function returns the second derivative spectral matrix
%

N = 8;

Np=N+1;
T = zeros([Np-2 Np]);
for m = 0:N-2
  p = m+2:2:N;
  T(m+1,p+1) = p.*(p.^2 - m^2)
  pause
end
T(1,:)=0.5*T(1,:)

