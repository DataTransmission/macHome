function Ds = spectd1(N,basisflag);
%
% This function returns the first derivative spectral matrix
% Usage is Ds = spectd2(N,basisflag)
% where N         = order of spectral truncation
%       basisflag = type of basis (Chebyshev=1, Legendre=2)
%

Np=N+1;
Ds = zeros([Np Np]);
switch (basisflag)
 case 1
   for m = 0:N-1
     p = m+1:2:N;
     Ds(m+1,p+1) = 2*p;
   end
   Ds(1,:)=0.5*Ds(1,:);
 case 2
   for m = 0:N-1
     p = m+1:2:N;
     Ds(m+1,p+1) = (2*m+1);
   end
end
