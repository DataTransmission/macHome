function T = spectd2(N,basisflag);
%
% This function returns the second derivative spectral matrix
% Usage is T = spectd2(N,basisflag)
% where N         = order of spectral truncation
%       basisflag = type of basis (Chebyshev=1, Legendre=2)
%

Np=N+1;
T = zeros([Np Np]);
switch (basisflag)
 case 1
   for m = 0:N-2
     p = m+2:2:N;
     T(m+1,p+1) = p.*(p.^2 - m^2);
   end
   T(1,:)=0.5*T(1,:);
 case 2
   for m = 0:N-2
     p = m+2:2:N;
     T(m+1,p+1) = (m+0.5)* (p.*(p+1) - m*(m+1));
   end
end
